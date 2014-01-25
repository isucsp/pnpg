/*
 * ====================================================================
 *
 *       Filename:  parFwdPrj.c
 *
 *    Description:  mexFunction to calculate parallel fwd projection
 *
 *        Version:  1.0
 *        Created:  09/12/2012 09:49:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Renliang Gu (), renliang@iastate.edu
 *   Organization:  Iowa State University
 *
 * ====================================================================
 */
#include "mex.h"
#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "parPrj.h"

/*  The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *inpr, *outpr;
    double *maskIdx;
    char* cmd;

    /*  check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }

    cmd=mxArrayToString(prhs[3]);
    inpr=mxGetPr(prhs[0]);
    maskIdx=mxGetPr(prhs[1]);
    size_t totPixels=mxGetNumberOfElements(prhs[1]);

    parConf config;
    config.bw=mxGetScalar(mxGetField(prhs[2],0,"bw"));
    config.nc=mxGetScalar(mxGetField(prhs[2],0,"nc"));
    config.nr=mxGetScalar(mxGetField(prhs[2],0,"nr"));
    config.np=mxGetNumberOfElements(mxGetField(prhs[2],0,"theta"));
    config.prjWidth=mxGetScalar(mxGetField(prhs[2],0,"prjWidth"));
    config.theta=mxGetPr(mxGetField(prhs[2],0,"theta"));

    if(!strcmp(cmd,"forward")){config.fwd=1;}
    else config.fwd=0;
    mxFree(cmd);
    /*printf("prjWidth=%d\tnp=%d\n",config.prjWidth, config.np);*/

    /*  create the output matrix */
    if(config.fwd){      /* forward projection */
        plhs[0] = mxCreateDoubleMatrix(config.prjWidth*config.np,1,mxREAL);
        outpr = mxGetPr(plhs[0]);
        parPrjConf(inpr, outpr, maskIdx, &config, totPixels);
    }else{
        plhs[0] = mxCreateDoubleMatrix(totPixels,1,mxREAL);
        outpr = mxGetPr(plhs[0]);
        parPrjConf(outpr, inpr, maskIdx, &config, totPixels);
    }
#if MULTHREAD
    parPrjThread();
#else
    parPrjRun();
#endif
    return;
}

