function [gradXYZ,gradTheta,gradPhi] = calcPotentialGradsCompiled(inField)
%CALCPORENTIALGRADSCOMPILED uses the .mex version of the potential gradient
%calculation functions, if available.
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

gradXYZ = zeros(length(inField.nCells),2);
gradTheta = zeros(length(inField.nCells),1);
for i = 1:length(inField.nCells) %The index of the cell alpha.
    %Get indices of cells that this cell (alpha) interacts with
    betInds = find(includeMat(i,:));
    xs = inField.xCells; xBets = xs(betInds);
    ys = inField.yCells; yBets = ys(betInds);
    ns = inField.nCells; nBets = ns(betInds);
    ls = inField.lCells; lBets = ls(betInds);
    thets = inField.thetCells; thetBets = thets(betInds);
    
    %Get dynamics
    
    if length(inField.aCells) >= 1
        xAlph = inField.xCells(i);
        yAlph = inField.yCells(i);
        lAlph = inField.lCells(i);
        nAlph = inField.nCells(i);
        thetAlph = inField.thetCells(i);
    end
    
    if strcmp(inField.boundConds,'none')
        [dUdx,dUdy,dUdthet] = mexCalcEnergyGradients(xBets,yBets,lBets,nBets,thetBets,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height);
    elseif strcmp(inField.boundConds,'periodic')
        [dUdx,dUdy,dUdthet] = mexCalcEnergyGradientsPeriodic(xBets,yBets,lBets,nBets,thetBets,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height);
    end
    
    gradXYZ(i,:) = -[sum(dUdx)/2,sum(dUdy)/2];
    gradTheta(i) = -sum(dUdthet)/2;
end