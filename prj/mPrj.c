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
#include "prj.h"

#ifdef __cplusplus
extern "C" 
#endif
extern struct prjConf* pConf;

/*  The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *imgt, *sinot;
    //double *maskIdx;
    char* cmd;
    if(nrhs!=3){
        mexPrintf("number of input arguments are: %d\n",nrhs);
        return;
    }
    cmd=mxArrayToString(prhs[2]);

    if(!strcmp(cmd,"config")){
        struct prjConf config;
        config.n=mxGetScalar(mxGetField(prhs[1],0,"n"));
        config.prjWidth=mxGetScalar(mxGetField(prhs[1],0,"prjWidth"));
        config.np=mxGetScalar(mxGetField(prhs[1],0,"np"));
        config.prjFull=mxGetScalar(mxGetField(prhs[1],0,"prjFull"));
        config.dSize=(float)mxGetScalar(mxGetField(prhs[1],0,"dSize"));
        config.effectiveRate=(float)mxGetScalar(mxGetField(prhs[1],0,"effectiveRate"));
        config.d=(float)mxGetScalar(mxGetField(prhs[1],0,"d"));
#if DEBUG
        mexPrintf("\nConfiguring the operator...\n");
        mexPrintf("config.n=%d\n",config.n);
        mexPrintf("config.prjWidth=%d\n",config.prjWidth);
        mexPrintf("config.np=%d\n",config.np);
        mexPrintf("config.prjFull=%d\n",config.prjFull);
        mexPrintf("config.dSize=%g\n",config.dSize);
        mexPrintf("config.effectiveRate=%g\n",config.effectiveRate);
        mexPrintf("config.d=%g\n",config.d);
#endif

        setup(config.n, config.prjWidth, config.np, config.prjFull, config.dSize, config.effectiveRate, config.d);
    }else{
        ft* img=(ft*)malloc(pConf->imgSize*sizeof(ft));
        ft* sino=(ft*)malloc(pConf->sinoSize*sizeof(ft));
        if(!strcmp(cmd,"forward")){      /* forward projection */
            plhs[0] = mxCreateNumericMatrix(pConf->prjWidth*pConf->np,1,mxDOUBLE_CLASS,mxREAL);
            imgt = mxGetPr(prhs[0]);
            sinot = mxGetPr(plhs[0]);
            for(int i=0; i<pConf->n; i++)
                for(int j=0; j<pConf->n; j++)
                    img[i*pConf->n+j]=(ft)imgt[i+j*pConf->n];
#if GPU
            gpuPrj(img, sino, FWD_BIT);
#else
            cpuPrj(img, sino, FWD_BIT);
#endif
            for(int i=0; i<pConf->sinoSize; i++)
                sinot[i]=sino[i];
        }else if(!strcmp(cmd,"backward")){
            plhs[0] = mxCreateNumericMatrix(pConf->n*pConf->n,1,mxDOUBLE_CLASS,mxREAL);
            imgt = mxGetPr(plhs[0]);
            sinot = mxGetPr(prhs[0]);
            for(int i=0; i<pConf->sinoSize; i++)
                sino[i]=(ft)sinot[i];
#if GPU
                gpuPrj(img, sino, BWD_BIT );
#else
                cpuPrj(img, sino, BWD_BIT );
#endif
            for(int i=0; i<pConf->n; i++)
                for(int j=0; j<pConf->n; j++)
                    imgt[i*pConf->n+j]=img[i+j*pConf->n];
        }else if(!strcmp(cmd,"FBP")){
            plhs[0] = mxCreateNumericMatrix(pConf->n*pConf->n,1,mxDOUBLE_CLASS,mxREAL);
            imgt = mxGetPr(plhs[0]);
            sinot = mxGetPr(prhs[0]);
            for(int i=0; i<pConf->sinoSize; i++)
                sino[i]=(ft)sinot[i];
#if GPU
                gpuPrj(img, sino, FBP_BIT );
#else
                cpuPrj(img, sino, FBP_BIT );
#endif
            for(int i=0; i<pConf->n; i++)
                for(int j=0; j<pConf->n; j++)
                    imgt[i*pConf->n+j]=img[i+j*pConf->n];
        }else{
            showSetup();
            mexPrintf("\nPrinting the current configuration ...\n");
            mexPrintf("config.n=%d\n",pConf->n);
            mexPrintf("config.prjWidth=%d\n",pConf->prjWidth);
            mexPrintf("config.np=%d\n",pConf->np);
            mexPrintf("config.prjFull=%d\n",pConf->prjFull);
            mexPrintf("config.dSize=%g\n",pConf->dSize);
            mexPrintf("config.effectiveRate=%g\n",pConf->effectiveRate);
            mexPrintf("config.d=%g\n",pConf->d);
        }
    }
    return;
}

