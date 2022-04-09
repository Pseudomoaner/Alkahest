function [gradXYZ,gradTheta] = calcPotentialGradsGPU(inField)
%CALCPORENTIALGRADSGPU uses the GPU-compatible version of the potential gradient
%calculation functions, if necessary hardware is available.
%
%   INPUTS:
%       -inField: The input WensinkField object
%
%   OUTPUTS:
%       -gradXYZ: The potential gradient, x, y and z components.
%       -gradTheta: The potential gradient, for reorientations in the
%       xy-plane
%       -gradPhi: The potential gradient, for reorientations in the polar
%       angle (wrt the xy-plane).
%
% Note that performance tests suggest that this code is actually slower
% than the .mex based versions of the potential gradient calculations, at
% least for systems of ~2000 rods. I have left it here in case greater
% efficiencies can be squeezed out of it in the future (or if we start
% running larger simulations, where the relative performance may increase)
%
%   Author: Oliver J. Meacock, (c) 2020

boundX = inField.xWidth/2;
boundY = inField.yHeight/2;
Height = inField.yHeight;
Width = inField.xWidth;
lam = inField.lam;
U0 = inField.U0;

includeMat = inField.cellDists < inField.distThresh;
includeMat(logical(diag(ones(length(includeMat),1)))) = 0;
includeMat = includeMat(1:length(inField.nCells),:); %Alpha rods should include all cells, beta rods all cells and all barrier rods

%To parallelize, need to feed data to the GPU in an efficient manner. In
%this section, I will construct Nx1 vectors encoding all the data needed to
%calculate the interactions between rods
noInts = sum(sum(includeMat));
xBets = zeros(noInts,1); xAlphs = zeros(noInts,1);
yBets = zeros(noInts,1); yAlphs = zeros(noInts,1);
nBets = zeros(noInts,1); nAlphs = zeros(noInts,1);
lBets = zeros(noInts,1); lAlphs = zeros(noInts,1);
thetBets = zeros(noInts,1); thetAlphs = zeros(noInts,1);
alphInds = zeros(noInts,1);

currInd = 1;
for i = 1:length(inField.nCells) %The index of the cell alpha.
    %Get indices of cells that this cell (alpha) interacts with
    betInds = find(includeMat(i,:));

    %Store data associated with these cells (beta)
    xs = inField.xCells; xBets(currInd:currInd + numel(betInds) - 1) = xs(betInds);
    ys = inField.yCells; yBets(currInd:currInd + numel(betInds) - 1) = ys(betInds);    
    ns = inField.nCells; nBets(currInd:currInd + numel(betInds) - 1) = ns(betInds);
    ls = inField.lCells; lBets(currInd:currInd + numel(betInds) - 1) = ls(betInds);
    thets = inField.thetCells; thetBets(currInd:currInd + numel(betInds) - 1) = thets(betInds);
    
    %Store data associated with alpha
    xAlphs(currInd:currInd + numel(betInds) - 1) = inField.xCells(i);
    yAlphs(currInd:currInd + numel(betInds) - 1) = inField.yCells(i);
    nAlphs(currInd:currInd + numel(betInds) - 1) = inField.nCells(i);
    lAlphs(currInd:currInd + numel(betInds) - 1) = inField.lCells(i);
    thetAlphs(currInd:currInd + numel(betInds) - 1) = inField.thetCells(i);
    alphInds(currInd:currInd + numel(betInds) - 1) = i;

    %Increment storage index
    currInd = currInd + numel(betInds);
end

%Create versions of these vectors on the GPU
xAlphsG = gpuArray(xAlphs);
yAlphsG = gpuArray(yAlphs);
nAlphsG = gpuArray(nAlphs);
lAlphsG = gpuArray(lAlphs);
thetAlphsG = gpuArray(thetAlphs);
alphIndsG = gpuArray(alphInds);

xBetsG = gpuArray(xBets);
yBetsG = gpuArray(yBets);
nBetsG = gpuArray(nBets);
lBetsG = gpuArray(lBets);
thetBetsG = gpuArray(thetBets);

%And versions of constants on the GPU
U0G = gpuArray(inField.U0);
lamG = gpuArray(inField.lam);
boundXG = gpuArray(boundX);
boundYG = gpuArray(boundY);
WidthG = gpuArray(Width);
HeightG = gpuArray(Height);

%Now use arrayfun to apply necessary calculations to each element of these
%vectors. Note we are using nested functions here, rather than hiving off
%functionality to separate scripts (as in other versions of the potential
%interaction code).
if strcmp(inField.boundConds,'none')
    [dUdxG,dUdyG,dUdthetG] = arrayfun(@gpuCalcEnergyGradients,xBetsG,yBetsG,lBetsG,nBetsG,thetBetsG,xAlphsG,yAlphsG,lAlphsG,nAlphsG,thetAlphsG,U0G,lamG);
elseif strcmp(inField.boundConds,'periodic')
    [dUdxG,dUdyG,dUdthetG] = arrayfun(@gpuCalcEnergyGradientsPeriodic,xBetsG,yBetsG,lBetsG,nBetsG,thetBetsG,xAlphsG,yAlphsG,lAlphsG,nAlphsG,thetAlphsG,U0G,lamG,boundXG,boundYG,WidthG,HeightG);
end

%Finally, transfer data back to the CPU and add up potential contributions to each rod alpha
gradXYg = -[accumarray(alphIndsG,dUdxG),accumarray(alphIndsG,dUdyG)]/2;
gradThetaG = -accumarray(alphIndsG,dUdthetG)/2;

gradXYZ = gather(gradXYg);
gradTheta = gather(gradThetaG);
end

function [dUdx,dUdy,dUdthet] = gpuCalcEnergyGradientsPeriodic(xBet,yBet,lBet,nBet,thetBet,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height)
    preFac = U0/(nAlph*nBet);
    
    dUdx = 0;
    dUdy = 0;
    
    dUdthet = 0;
    
    %Pairwise comparison of segments - note vectorization doesn't work on
    %GPUs
    for i = 1:nAlph
        for j = 1:nBet
            alphPos = i - ((nAlph+1)/2);
            betPos = j - ((nBet+1)/2);
            
            xiAlph = xAlph + (lAlph * alphPos * cos(thetAlph));
            xjBet = xBet + (lBet * betPos * cos(thetBet));
            yiAlph = yAlph + (lAlph * alphPos * sin(thetAlph));
            yjBet = yBet + (lBet * betPos * sin(thetBet));
            
            x = xiAlph - xjBet;
            y = yiAlph - yjBet;
            
            periodX = abs(x) > boundX; %These segements are closer in the wrap-around x direction.
            periodY = abs(y) > boundY; %Likewise for y.
            
            %Find the distance between these segments in the wrap-around direction
            if periodX
                tmpX = x;
                absX = abs(tmpX);
                sgnX = tmpX./absX;
                x = -sgnX*(Width - absX);
            end
            if periodY
                tmpY = y;
                absY = abs(tmpY);
                sgnY = tmpY./absY;
                y = -sgnY*(Height - absY);
            end
                        
            r = (x^2 + y^2)^0.5;
            invR = 1/r;
            
            drdx = invR*x;
            drdy = invR*y;
            
            drdthet = lAlph*alphPos*invR*(cos(thetAlph)*y - sin(thetAlph)*x);
            
            postFac = (exp(-r/lam) * (lam + r)) / (lam * r^2);
            
            dUdx = dUdx + preFac * drdx * postFac;
            dUdy = dUdy + preFac * drdy * postFac;
            
            dUdthet = dUdthet + preFac * drdthet * postFac;
        end
    end
end

function [dUdx,dUdy,dUdthet] = gpuCalcEnergyGradients(xBet,yBet,lBet,nBet,thetBet,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam)
    preFac = U0/(nAlph*nBet);
    
    dUdx = 0;
    dUdy = 0;
    
    dUdthet = 0;
    
    %Pairwise comparison of segments - note vectorization doesn't work on
    %GPUs
    for i = 1:nAlph
        for j = 1:nBet
            alphPos = i - ((nAlph+1)/2);
            betPos = j - ((nBet+1)/2);
            
            xiAlph = xAlph + (lAlph * alphPos * cos(thetAlph));
            xjBet = xBet + (lBet * betPos * cos(thetBet));
            yiAlph = yAlph + (lAlph * alphPos * sin(thetAlph));
            yjBet = yBet + (lBet * betPos * sin(thetBet));

            x = xiAlph - xjBet;
            y = yiAlph - yjBet;
    
            r = (x^2 + y^2)^0.5;
            invR = 1/r;
            
            drdx = invR*x;
            drdy = invR*y;
            
            drdthet = lAlph*alphPos*invR*(cos(thetAlph)*y - sin(thetAlph)*x);
            
            postFac = (exp(-r/lam) * (lam + r)) / (lam * r^2);
            
            dUdx = dUdx + preFac * drdx * postFac;
            dUdy = dUdy + preFac * drdy * postFac;
            
            dUdthet = dUdthet + preFac * drdthet * postFac;
        end
    end
end