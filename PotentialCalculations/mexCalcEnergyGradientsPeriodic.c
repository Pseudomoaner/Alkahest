#include <mex.h>
#include <math.h>
#include <matrix.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //xBets,yBets,lBets,nBets,thetBets,xAlph,yAlph,lAlph,nAlph,thetAlph,U0,lam,boundX,boundY,Width,Height are the input variables
    double *xBets, *yBets, *lBets, *thetBets, *betInds, *nBets, *dUdx, *dUdy, *dUdthet;
    double xAlph, yAlph, lAlph, thetAlph, U0, lam, boundX, boundY, Width, Height;
    int  nAlph, noBets;
    
    //Things declared for the energy calculation itself
    int betInd, bet, nBet, alphSeg, betSeg;
    double xBet, yBet, lBet, thetBet, preFac, postFac, drdx, drdy, drdthet, x, y, r, rInv, xiAlph, xjBet, yiAlph, yjBet, alphPos, betPos, absX, absY, sgnX, sgnY; 
    
    if (nrhs != 16) {
        mexErrMsgTxt("Need 16 (!) input arguments!");
    }
    
    //Unload the input variables
    xBets = mxGetPr(prhs[0]);
    yBets = mxGetPr(prhs[1]);
    lBets = mxGetPr(prhs[2]);
    nBets = mxGetPr(prhs[3]);
    thetBets = mxGetPr(prhs[4]);
    xAlph = mxGetScalar(prhs[5]);
    yAlph = mxGetScalar(prhs[6]);
    lAlph = mxGetScalar(prhs[7]);
    nAlph = (int)mxGetScalar(prhs[8]);
    thetAlph = mxGetScalar(prhs[9]);
    U0 = mxGetScalar(prhs[10]);
    lam = mxGetScalar(prhs[11]);
    boundX = mxGetScalar(prhs[12]);
    boundY = mxGetScalar(prhs[13]);
    Width = mxGetScalar(prhs[14]);
    Height = mxGetScalar(prhs[15]);
    
    noBets = mxGetM(prhs[1]);
    
    //Create output matrices:
    plhs[0] = mxCreateDoubleMatrix(1,noBets,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,noBets,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,noBets,mxREAL);
    dUdx = mxGetPr(plhs[0]);
    dUdy = mxGetPr(plhs[1]);
    dUdthet = mxGetPr(plhs[2]);
    
    //Use the main energy determination function
    for (bet = 0; bet < noBets; bet++) { //For each other cell beta
        
        xBet = xBets[bet];
        yBet = yBets[bet];
        lBet = lBets[bet];
        nBet = nBets[bet];
        thetBet = thetBets[bet];
        
        preFac = U0/(nAlph * nBet);
        
        //Pairwise comparison of segments
        for (alphSeg = 0; alphSeg < nAlph; alphSeg++) {
            alphPos = (double)(alphSeg) - ((double)(nAlph-1)/2);
            for (betSeg = 0; betSeg < nBet; betSeg++) {
                betPos = (double)(betSeg) - ((double)(nBet-1)/2);
                
                xiAlph = xAlph + (lAlph * alphPos * cos(thetAlph));
                xjBet = xBet + (lBet * betPos * cos(thetBet));
                yiAlph = yAlph + (lAlph * alphPos * sin(thetAlph));
                yjBet = yBet + (lBet * betPos * sin(thetBet));
                
                x = xiAlph - xjBet;
                y = yiAlph - yjBet;
                
                //Deal with the periodic boundaries
                absX = fabs(x);
                absY = fabs(y);
                
                if (absX > boundX) { //These segements are closer in the wrap-around x direction.
                    sgnX = x / absX;
                    x = -sgnX * (Width - absX);
                }
                if (absY > boundY) { //Likewise for y.
                    sgnY = y / absY;
                    y = -sgnY * (Height - absY);
                }
                
                //Distance between two segmenets (in periodic space)
                r = sqrt(x * x + y * y);
                rInv = 1/r;
                
                drdx = rInv * x;
                drdy = rInv * y;
                
                drdthet = rInv * lAlph * alphPos * (cos(thetAlph) * y - sin(thetAlph) * x);
                
                postFac = (exp(-r / lam) * (lam + r)) / (lam * r * r);
                
                dUdx[bet] = dUdx[bet] + (preFac * drdx * postFac);
                dUdy[bet] = dUdy[bet] + (preFac * drdy * postFac);
                
                dUdthet[bet] = dUdthet[bet] + (preFac * drdthet * postFac);
            }
        }
    }

    return;
}