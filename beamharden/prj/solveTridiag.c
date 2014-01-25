/*
 * ====================================================================
 *
 *       Filename:  solveTridiag.c
 *
 *    Description:  Solve a general tridiagonal system
 *
 *        Version:  1.0
 *        Created:  09/29/2012 11:07:06 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Organization:  Iowa State University
 *
 * ====================================================================
 */

#include "mex.h"
#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#define DEBUG 0

void solveMatrix(int n, double *a, double *b, double *c, double *vv, double *x){
    /**
     * n - number of equations
     * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
     * b - the main diagonal
     * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
     * v - right part
     * x - the answer
     */
    int i;
    double *v=vv;

#if 0
    v=calloc(n,sizeof(double));
    memcpy(v,vv,n*sizeof(double));
#endif
    for (i = 1; i < n; i++){
        double m = a[i]/b[i-1];
        b[i] = b[i] - m * c[i - 1];
        v[i] = v[i] - m*v[i-1];
    }

    x[n-1] = v[n-1]/b[n-1];

    for (i = n - 2; i >= 0; --i)
        x[i] = (v[i] - c[i] * x[i+1]) / b[i];
#if 0
    free(v);
#endif
}

int forward(int n, double *a, double *b, double *c, double *v, double *x){
    int i;
    v[0]=b[0]*x[0]+c[0]*x[1];
    for(i=1; i<n-1; i++)
        v[i]=a[i]*x[i-1]+b[i]*x[i]+c[i]*x[i+1];
    v[n-1]=a[n-1]*x[n-2]+b[n-1]*x[n-1];
}

/*  The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *a, *b, *c, *v, *outpr;

    /*  check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }

    a=mxGetPr(prhs[0]);
    b=mxGetPr(prhs[1]);
    size_t N=mxGetNumberOfElements(prhs[1]);
    c=mxGetPr(prhs[2]);
    v=mxGetPr(prhs[3]);
    if(N!=mxGetNumberOfElements(prhs[3]))
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","The demension doesn't match");
    char *cmd=mxArrayToString(prhs[4]);

    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    outpr = mxGetPr(plhs[0]);
    if(!strcmp(cmd,"inv")){
#if DEBUG
        printf("solve the inverse problem\n");
#endif
        solveMatrix(N,a,b,c,v,outpr);
    }else
        forward(N,a,b,c,outpr,v);
    return;
}
