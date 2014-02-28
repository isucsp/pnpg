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

extern struct prjConf* pConf;

/*  The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *imgt, *sinot;
    //double *maskIdx;
    char* cmd;
    if(nrhs!=3){
        printf("number of input arguments are: %d\n",nrhs);
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
        printf("\nConfiguring the operator...\n");
        printf("config.n=%d\n",config.n);
        printf("config.prjWidth=%d\n",config.prjWidth);
        printf("config.np=%d\n",config.np);
        printf("config.prjFull=%d\n",config.prjFull);
        printf("config.dSize=%g\n",config.dSize);
        printf("config.effectiveRate=%g\n",config.effectiveRate);
        printf("config.d=%g\n",config.d);

        setup(config.n, config.prjWidth, config.np, config.prjFull, config.dSize, config.effectiveRate, config.d);
    }else{
        ft img[pConf->imgSize];
        ft sino[pConf->sinoSize];
        if(!strcmp(cmd,"forward")){      /* forward projection */
            plhs[0] = mxCreateNumericMatrix(pConf->prjWidth*pConf->np,1,mxDOUBLE_CLASS,mxREAL);
            imgt = mxGetPr(prhs[0]);
            sinot = mxGetPr(plhs[0]);
            for(int i=0; i<pConf->n; i++)
                for(int j=0; j<pConf->n; j++)
                    img[i*pConf->n+j]=(ft)imgt[i+j*pConf->n];
#if GPU
            gpuPrj(img, sino, RENEW_MEM | FWD_BIT);
#else
            cpuPrj(img, sino, RENEW_MEM | FWD_BIT);
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
                gpuPrj(img, sino, RENEW_MEM );
#else
                cpuPrj(img, sino, RENEW_MEM );
#endif
            for(int i=0; i<pConf->n; i++)
                for(int j=0; j<pConf->n; j++)
                    imgt[i*pConf->n+j]=img[i+j*pConf->n];
        }else{
            showSetup();
            printf("\nPrinting the current configuration ...\n");
            printf("config.n=%d\n",pConf->n);
            printf("config.prjWidth=%d\n",pConf->prjWidth);
            printf("config.np=%d\n",pConf->np);
            printf("config.prjFull=%d\n",pConf->prjFull);
            printf("config.dSize=%g\n",pConf->dSize);
            printf("config.effectiveRate=%g\n",pConf->effectiveRate);
            printf("config.d=%g\n",pConf->d);
        }
    }
    return;
}

