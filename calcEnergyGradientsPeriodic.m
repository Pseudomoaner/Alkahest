function [dUdx,dUdy,dUdthet] = calcEnergyGradientsPeriodic(xBets,yBets,lBets,nBets,thetBets,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height)
%  CALCENERGYGRADIENTSPERIODIC is a native MATLAB-compatable implementation
%  of mexCalcEnergyGradientsPeriodic.c. Best to compile that if you can,
%  but this will work. This function calculates the interaction between a
%  rod alpha and each of its neighbours, beta.
%  
%  For full details, see:
% 
%  Wensink, H. H., & Löwen, H. (2012). Emergent states in dense systems of 
%  active rods: From swarming to turbulence. Journal of Physics Condensed 
%  Matter, 24(46). https://doi.org/10.1088/0953-8984/24/46/464130
% 
%  Author: Oliver J. Meacock

dUdx = zeros(size(xBets,1),1);
dUdy = zeros(size(xBets,1),1);
dUdthet = zeros(size(xBets,1),1);

for bet = 1:size(xBets,1)
    xBet = xBets(bet);
    yBet = yBets(bet);
    lBet = lBets(bet);
    nBet = nBets(bet);
    thetBet = thetBets(bet);
    
    preFac = U0/(nAlph*nBet);
    
    %Pairwise comparison of segments
    alphPos = repmat((1:nAlph) - ((nAlph+1)/2),[nBet,1]); %Rows are the alpha index, columns the beta index.
    betPos = repmat((1:nBet)' - ((nBet+1)/2),[1,nAlph]);
    
    xiAlph = xAlph + lAlph*alphPos*cos(thetAlph);
    xjBet = xBet + lBet*betPos*cos(thetBet);
    yiAlph = yAlph + lAlph*alphPos*sin(thetAlph);
    yjBet = yBet + lBet*betPos*sin(thetBet);
    
    x = xiAlph - xjBet;
    y = yiAlph - yjBet;
    
    periodX = abs(x) > boundX; %These segements are closer in the wrap-around x direction.
    periodY = abs(y) > boundY; %Likewise for y.
    
    %Find the distance between these segments in the wrap-around direction
    tmpX = x(periodX);
    tmpY = y(periodY);
    
    absX = abs(tmpX);
    sgnX = tmpX./absX;
    absY = abs(tmpY);
    sgnY = tmpY./absY;
    
    x(periodX) = -sgnX.*(Width - absX);
    y(periodY) = -sgnY.*(Height - absY);
    
    r = (x.^2 + y.^2).^0.5;
    
    drdx = (1./r).*x;
    drdy = (1./r).*y;
    
    yTerm = lAlph*cos(thetAlph)*alphPos.*y;
    xTerm = lAlph*sin(thetAlph)*alphPos.*x;
    drdthet = (1./r).*(yTerm - xTerm);
    postFac = (exp(-r/lam) .* (lam + r)) ./ (lam * r.^2);
    
    dUdx(bet) = sum(sum(preFac * drdx .* postFac));
    dUdy(bet) = sum(sum(preFac * drdy .* postFac));
    
    dUdthet(bet) = sum(sum(preFac * drdthet .* postFac));
end
