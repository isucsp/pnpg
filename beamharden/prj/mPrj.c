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
    ft *inpr, *outpr;
    //double *maskIdx;
    char* cmd;

    if(nrhs!=3){
        printf("number of input arguments are: %d\n",nrhs);
        return;
    }
    cmd=mxArrayToString(prhs[2]);

    /*  create the output matrix */
    if(!strcmp(cmd,"forward")){      /* forward projection */
        inpr=(ft*)mxGetPr(prhs[0]);
        if(sizeof(ft)==sizeof(float))
            plhs[0] = mxCreateNumericMatrix(pConf->prjWidth*pConf->np,1,mxSINGLE_CLASS,mxREAL);
        else
            plhs[0] = mxCreateNumericMatrix(pConf->prjWidth*pConf->np,1,mxDOUBLE_CLASS,mxREAL);
        outpr = (ft *)mxGetPr(plhs[0]);
#if GPU
        gpuPrj(inpr, outpr, RENEW_MEM | FWD_BIT);
#else
        cpuPrj(inpr, outpr, RENEW_MEM | FWD_BIT);
#endif
    }else if(!strcmp(cmd,"backward")){
        inpr=(ft*)mxGetPr(prhs[0]);
        if(sizeof(ft)==sizeof(float)){
            plhs[0] = mxCreateNumericMatrix(pConf->n*pConf->n,1,mxSINGLE_CLASS,mxREAL);
            outpr = (ft*)mxGetData(plhs[0]);
        }else{
            plhs[0] = mxCreateNumericMatrix(pConf->n*pConf->n,1,mxDOUBLE_CLASS,mxREAL);
            outpr = (ft*)mxGetPr(plhs[0]);
        }
#if GPU
        gpuPrj(outpr, inpr, RENEW_MEM );
#else
        cpuPrj(outpr, inpr, RENEW_MEM );
#endif
    }else if(!strcmp(cmd,"config")){
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
    return;
}

