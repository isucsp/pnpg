/*
 * ====================================================================
 *
 *       Filename:  parPrj.c
 *
 *    Description:  calculate the parallel projection
 *
 *        Version:  1.0
 *        Created:  09/12/2012 06:09:58 PM
 *       Revision:  1.0
 *       Compiler:  gcc
 *
 *         Author:  Renliang Gu (), renliang@iastate.edu
 *   Organization:  Iowa State University
 *
 * ====================================================================
 */

/*#include <unistd.h>*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "common/thread.h"
#include "parPrj.h"

const double pi=3.14159265;
int nr=1024, nc=1024; /* The origin is top-left corner*/
int xC, yC, pC, totPixels;
double beamWidth=1;
double hbeamW;
double *img, *prj, *maskIdx;
parConf *opt;
int blockSz, angleSz;

int parPrjConf(double *t_img, double *t_prj, double *t_maskIdx, parConf *t_opt,\
        size_t t_totPixels){
    img=t_img; prj=t_prj; maskIdx=t_maskIdx; opt=t_opt; totPixels=t_totPixels;
    beamWidth=opt->bw;
    nr=opt->nr; nc=opt->nc;
    hbeamW=beamWidth/2;
    xC=nc/2; yC=nr/2; pC=xC>yC? xC : yC;
}

int parPrjThread(){
    int ncpu=32;
    int res;
    CUTThread *a_thread;
    void *thread_result;
    size_t i=0; 

    a_thread=calloc(ncpu,sizeof(pthread_t));
    if(opt->fwd){
        angleSz=(opt->np+ncpu-1)/(ncpu);
        /*printf("angleSz=%d given ncpu=%d\n",angleSz,ncpu);*/
        for(i=0; i<ncpu; i++)
            a_thread[i]=start_thread(parPrjFor, (void *)i);
        /*printf("Waiting for thread to finish...\n");*/
    }else{
        blockSz=(totPixels+ncpu-1)/(ncpu);
        /*printf("blockSz=%d given ncpu=%d\n",blockSz,ncpu);*/
        for(i=0; i<ncpu; i++){
            a_thread[i]=start_thread(parPrjBack, (void *)i);
        }
    }
    /*printf("Waiting for thread to finish...\n");*/
    wait_for_threads(a_thread,ncpu);
    free(a_thread);
    if(opt->fwd)
        for(i=0; i<opt->np*opt->prjWidth; i++)
            prj[i]/=beamWidth;
    else
        for(i=0; i<totPixels; i++)
            img[i]/=beamWidth;
    return 0;
}

int parPrjRun(){
    int i=0; 
    for(i=0; i<opt->np; i++) parPrj(0,totPixels,i);
    return 0;
}

CUT_THREADPROC parPrjFor(void *arg){
    int i=0; 
    int istart=((size_t)arg)*angleSz;
    int iend=istart+angleSz;
    if(iend>opt->np) iend=opt->np;
    for(i=istart; i<iend; i++) parPrj(0,totPixels,i);
    CUT_THREADEND;
}

CUT_THREADPROC parPrjBack(void *arg){
    int i=0; 
    int istart=((size_t)arg)*blockSz;
    int iend=istart+blockSz;
    if(iend>totPixels) iend=totPixels;
    for(i=0; i<opt->np; i++) parPrj(istart,iend,i);
    CUT_THREADEND;
}

void *parPrj(int istart, int iend, int itheta){
    double *localPrj=prj+opt->prjWidth*itheta;
    double theta=opt->theta[itheta];
#if DEBUG
    printf("prject theta=%g\n",theta);
#endif
    double temp=abs(fmod(theta,90)-45)*pi/180; /* test value of temp */
    double height=1/sin(temp+pi/4);
    double bias1=sqrt(2)/2*sin(temp);
    double bias2=sqrt(2)/2*cos(temp);

    double areaT=0.5*height*(bias2-bias1);
    double areaR=height*bias1*2;
    double weight, dist;
    double ub, lb;
    int x,y,k,i;

    /*printf("#dist\tweight\n");*/
    for(i=istart; i<iend; i++){
        /* row index: *maskIdx%nr
           col index: *maskIdx/nr */
        x=(((int)maskIdx[i])/nr)-xC; y=(((int)maskIdx[i])%nr)-yC;
        dist=x*cos(theta*pi/180)+y*sin(theta*pi/180);
        ub=floor(dist+bias2+beamWidth/2-EPS);
        lb=ceil(dist-bias2-beamWidth/2+EPS);
#if DEBUG
        printf("dist=%g\tub=%g\tlb=%g\n",dist,ub,lb);
#endif
        for(k=(int)lb; k<=ub; k++){
            if(k-hbeamW<=dist-bias2)
                if(k+hbeamW<=dist-bias1)
                    weight=areaT*pow((k+hbeamW-(dist-bias2))/(bias2-bias1),2);
                else if(k+hbeamW<=dist+bias1)
                    weight=areaT+height*(k+hbeamW-(dist-bias1));
                else
                    weight=areaT*2+areaR-\
                           areaT*pow((dist+bias2-k-hbeamW)/(bias2-bias1),2);
            else if(k-hbeamW<=dist-bias1)
                if(k+hbeamW<=dist-bias1)
                    weight=height/(bias2-bias1)*(k-dist+bias2)*beamWidth;
                else if(k+hbeamW<=dist+bias1)
                    weight=areaT*(1-pow((k-hbeamW-dist+bias2)/(bias2-bias1),2))\
                           +height*(k+hbeamW-(dist-bias1));
                else if(k+hbeamW<=dist+bias2)
                    weight=areaT*(1-pow((k-hbeamW-dist+bias2)/(bias2-bias1),2))\
                           +areaR+\
                           areaT*(1-pow((dist+bias2-k-hbeamW)/(bias2-bias1),2));
                else
                    weight=areaT*(1-pow((k-hbeamW-dist+bias2)/(bias2-bias1),2))\
                           +areaR+areaT;
            else if(k-hbeamW<=dist+bias1)
                if(k+hbeamW<=dist+bias1)
                    weight=height*beamWidth;
                else if(k+hbeamW<=dist+bias2)
                    weight=height*(dist+bias1-k+hbeamW)+\
                           areaT*(1-pow((dist+bias2-k-hbeamW)/(bias2-bias1),2));
                else
                    weight=height*(dist+bias1-k+hbeamW)+areaT;
            else{
                if(k+hbeamW<=dist+bias2)
                    weight=(height)/(bias2-bias1)*(dist+bias2-k)*beamWidth;
                else
                    weight=0.5*height/(bias2-bias1)*pow(dist+bias2-k+hbeamW,2);
            }
            /*printf("%d\t%e\n",k+pC,weight);*/
            if((k<-pC)||((k+pC)>=nc)){
                /* printf("Mask is not proper: %d\n",k+pC); */
                continue;
            }
            if(weight!=0){
                if(opt->fwd)         /* Forward projection */
                    localPrj[k+pC] += weight * img[i];
                else            /* Backward projection */
                    img[i] += localPrj[k+pC] * weight;
            }
        }
    }
    return 0;
}
