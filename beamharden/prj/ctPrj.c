/*
 * ====================================================================
 *
 *       Filename:  ctPrj.cu
 *
 *    Description:  use GPU to calculate parallel beam, fan beam projection 
 *    and back projection
 *
 *        Version:  1.0
 *        Created:  09/12/2012 06:09:58 PM
 *       Revision:  1.0
 *       Compiler:  nvcc
 *
 *         Author:  Renliang Gu (), renliang@iastate.edu
 *   Organization:  Iowa State University
 *
 * ====================================================================
 */

/*#include <unistd.h>*/
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <pthread.h>
#include "ctPrj.h"
#include "common/cpu_bitmap.h"
//#include "mex.h"

#define DIM   1024
#define XC(N) (N/2)
#define YC(N) (N/2)

typedef float ft;

const ft PI = 3.14159265;
const ft SQRT2 = sqrt(2);
int nr=1024, nc=1024; /* The origin is top-left corner*/
int xC, yC, pC, totPixels;
ft beamWidth=1;
ft hbeamW;
ft *img, *prj, *maskIdx, *sino;
cpuPrjConf conf;
int blockSz, angleSz;

/*
   */
ft getWeight(ft dist, ft beamWidth, ft cosR, ft sinR){
    // range of gamma is in [-alpha,pi/4+alpha], where alpha is small
    //ft temp=abs(fmod(theta,90)-45)*PI/180; /* test value of temp */
    ft height=1/MAX(cosR,sinR);
    ft bias1, bias2;             //bias2 should be larger.
    ft temp1=fabs(cosR+sinR)/2;
    ft temp2=fabs(cosR-sinR)/2;
    //if(cosR==1 && sinR==0) printf("temp1=%f, temp2=%f\n",temp1,temp2);
    if(temp1>temp2){
        bias2=temp1; bias1=temp2;
    }else{
        bias1=temp1; bias2=temp2;
    }
    //if(cosR==1 && sinR==0) printf("height=%f, bias1=%f, bias2=%f\n",height,bias1,bias2);

    ft areaT=0.5*height*(bias2-bias1);
    ft areaR=height*bias1*2;
    ft area=height*(bias1+bias2);
    ft weight=0;

    // if pixel is to the right of the ray, dist > 0, otherwise dist<=0
    ft lfoot=dist-bias2, rfoot=dist+bias2;
    ft lshoulder = dist-bias1, rshoulder=dist+bias1;
    ft hbeamW = beamWidth/2;
    if(-hbeamW<=lfoot){
        if (hbeamW<=lfoot) weight=0;
        if(hbeamW<=lshoulder)
            weight=areaT*pow((hbeamW-(dist-bias2))/(bias2-bias1),2);
        else if(hbeamW<=rshoulder)
            weight=areaT+height*(hbeamW-(dist-bias1));
        else if(hbeamW<=rfoot)
            weight=areaT*2+areaR-\
                   areaT*pow((dist+bias2-hbeamW)/(bias2-bias1),2);
        else weight=area;
    }else if(-hbeamW<=lshoulder){
        if(hbeamW<=lshoulder)
            weight=height/(bias2-bias1)*(-dist+bias2)*beamWidth;
        else if(hbeamW<=rshoulder)
            weight=areaT*(1-pow((-hbeamW-dist+bias2)/(bias2-bias1),2))\
                   +height*(hbeamW-(dist-bias1));
        else if(hbeamW<=rfoot)
            weight=areaT*(1-pow((-hbeamW-dist+bias2)/(bias2-bias1),2))\
                   +areaR+\
                   areaT*(1-pow((dist+bias2-hbeamW)/(bias2-bias1),2));
        else
            weight=areaT*(1-pow((-hbeamW-dist+bias2)/(bias2-bias1),2))\
                   +areaR+areaT;
    }else if(-hbeamW<=rshoulder){
        if(hbeamW<=rshoulder)
            weight=height*beamWidth;
        else if(hbeamW<=rfoot)
            weight=height*(dist+bias1+hbeamW)+\
                   areaT*(1-pow((dist+bias2-hbeamW)/(bias2-bias1),2));
        else
            weight=height*(dist+bias1+hbeamW)+areaT;
    }else{
        if(-hbeamW<=rfoot)
            if(hbeamW<=rfoot)
                weight=(height)/(bias2-bias1)*(dist+bias2)*beamWidth;
            else
                weight=0.5*height/(bias2-bias1)*pow(dist+bias2+hbeamW,2);
    }
    return weight;
}

void rayDrive(int thetaIdx, int tIdx){
    //printf("entering rayDrive\n");
    // detector array is of odd size with the center at the middle
    //int tIdx, thetaIdx; // index counts from 0
    int sinoIdx;
    //printf("tIdx=%d, thetaIdx=%d\n",tIdx, thetaIdx);

    if(tIdx==0 && conf.prjWidth%2==0){
        //printf("thetaIdx=%d\n",thetaIdx);
        //printf("(x,y)=(%d,%d)\n",gridDim.x,gridDim.y);
        sinoIdx=thetaIdx*conf.prjWidth+tIdx;
        sino[sinoIdx]=0;
        return; //make sure the size is odd
    }

    ft theta;       // the current angle, range in [0 45]
    ft t;           // position of current detector
    ft d;           // the distance from rotating center to X-ray source in pixel
    ft dSize;       // size of each detector in pixel

    int N;   // N is even// size of image is NxN, and centers at (N/2,N/2)
    int pC = conf.prjWidth/2;

    if(conf.nc==conf.nr) N=conf.nc;
    else
        N=0;
    theta = thetaIdx*2*PI/conf.prjFull;
    dSize = conf.dSize;
    t = (tIdx-pC)*dSize;
    d = conf.d;

    ft cosT, sinT;  // cosine and sine of theta
    cosT=cos(theta); sinT=sin(theta);

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft bq=sqrt(d*d+t*t);
    ft cosB=d/bq, sinB=t/bq; //, tanB=t/d;
    ft cosR=cosB*cosT-sinB*sinT;
    ft sinR=sinB*cosT+cosB*sinT;
    ft beamWidth = dSize*conf.effectiveRate;
    ft beamwidthCosB=beamWidth*cosB;
    ft beamwidthCosBOverD=beamwidthCosB/d;
    ft sinTBeamwidthCosBOverD=sinT*beamwidthCosBOverD;

    //if(tIdx==pC && thetaIdx==30){
    //    printf("theta=%f,cosT=%f, sinT=%f\n",theta,cosT,sinT);
    //    printf("cosB=%f, sinB=%f\n",cosB,sinB);
    //}

    ft slopeXYc=(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft slopeXYl=((t-beamWidth/2)*cosT+d*sinT)/((t-beamWidth/2)*sinT-d*cosT);
    ft slopeXYr=((t+beamWidth/2)*cosT+d*sinT)/((t+beamWidth/2)*sinT-d*cosT);
    ft xc = -t*d/(t*sinT-d*cosT);
    ft xl = -(t-beamWidth/2)*d/((t-beamWidth/2)*sinT-d*cosT);
    ft xr = -(t+beamWidth/2)*d/((t+beamWidth/2)*sinT-d*cosT);
    int x,y=N/2-1;
    xl = xl + (y-1/2)*slopeXYl;
    xc = xc + y*slopeXYc;
    xr = xr + (y+1/2)*slopeXYr;
    int dxl=round(MIN(xl,xr)), dxr=round(MAX(xl,xr));
    dxl=MAX(dxl,-(N-1)/2); dxr=MIN(dxr,(N-1)/2);

    ft bwPre, bw, dist, weight;
    int temp, imgIdx;
    while(y>-N/2){
        for(x=dxl,
                bwPre = (x*sinT-y*cosT),
                bw = beamwidthCosB+bwPre*beamwidthCosBOverD;
                x<=dxr;
                x++,
                bw+=sinTBeamwidthCosBOverD){

            dist=x*cosR+y*sinR-d*t/bq;
            //if(thetaIdx==0 && tIdx==pC) printf("dist=%f, bw=%f, cosR=%f, sinR=%f\n",dist,bw,cosR,sinR);
            weight=getWeight(dist,bw,cosR,sinR);
            if(tIdx==pC && thetaIdx==0){
                //printf("(%d,%d)\n",x,y);
                //printf("%f,%f;\n",dist,weight);
            }

            sinoIdx=thetaIdx*conf.prjWidth+tIdx;
            imgIdx = (y+N/2)*N+x+N/2;
            sino[sinoIdx]+=weight*img[imgIdx];
            
            temp=thetaIdx+conf.prjFull/4;
            if(temp<conf.np){
                sinoIdx=temp*conf.prjWidth+tIdx;
                imgIdx = (x+N/2)*N-y+N/2;
                sino[sinoIdx]+=weight*img[imgIdx];
            }
            temp+=conf.prjFull/4;
            if(temp<conf.np){
                sinoIdx=temp*conf.prjWidth+tIdx;
                imgIdx = (-y+N/2)*N-x+N/2;
                sino[sinoIdx]+=weight*img[imgIdx];
            }
            temp+=conf.prjFull/4;
            if(temp<conf.np){
                sinoIdx=temp*conf.prjWidth+tIdx;
                imgIdx = (-x+N/2)*N+y+N/2;
                sino[sinoIdx]+=weight*img[imgIdx];
            }

            if(thetaIdx!=0 && thetaIdx!=conf.prjFull/8){
                temp=conf.prjFull/2-thetaIdx;
                if(temp<conf.np){
                    sinoIdx=temp*conf.prjWidth-tIdx+2*pC;
                    imgIdx = (-y+N/2)*N+x+N/2;
                    sino[sinoIdx]+=weight*img[imgIdx];
                }
                temp=3*conf.prjFull/4-thetaIdx;
                if(temp<conf.np){
                    sinoIdx=temp*conf.prjWidth-tIdx+2*pC;
                    imgIdx = (x+N/2)*N+y+N/2;
                    sino[sinoIdx]+=weight*img[imgIdx];
                }
                temp=conf.prjFull-thetaIdx;
                if(temp<conf.np){
                    sinoIdx=temp*conf.prjWidth-tIdx+2*pC;
                    imgIdx = (y+N/2)*N-x+N/2;
                    sino[sinoIdx]+=weight*img[imgIdx];
                }
                temp=conf.prjFull/4-thetaIdx;
                if(temp<conf.np){
                    sinoIdx=temp*conf.prjWidth-tIdx+2*pC;
                    imgIdx = (-x+N/2)*N-y+N/2;
                    sino[sinoIdx]+=weight*img[imgIdx];
                }
            }
        }
        xl -= slopeXYl; xc -= slopeXYc; xr -= slopeXYr;
        dxl=round(MIN(xl,xr)), dxr=round(MAX(xl,xr));
        dxl=MAX(dxl,-(N-1)/2); dxr=MIN(dxr,(N-1)/2);
        y--;
    }
}

void *parPrjFor(void *arg){
    //printf("test\n");
    int i=0; 
    int angleIdx=((size_t)arg);
    for(int tIdx = 0; tIdx < conf.prjWidth; tIdx++){
        //printf("angleIdx=%d, tIdx=%d\n",angleIdx, tIdx);
        rayDrive(angleIdx,tIdx);
    }
    return 0;
}

int multiThread(){
    conf.sec = MIN(conf.np, conf.prjFull/8+1);
    //conf.nthreadx = ((conf.sec+conf.angleStep-1)/conf.angleStep);
    //conf.nthready = ((conf.prjWidth+conf.prjStep-1)/conf.prjStep);
    int ncpu = conf.sec;
    int res;
    pthread_t *a_thread;
    void *thread_result;
    size_t i=0; 


    a_thread=(pthread_t*)calloc(ncpu,sizeof(pthread_t));
    if(conf.fwd){
        printf("forward\n");
        if(conf.rayDrive){
            printf("angleSz=%d given ncpu=%d\n",angleSz,ncpu);
            for(i=0; i<ncpu; i++){
                res = pthread_create(&a_thread[i], NULL, parPrjFor, (void *)i);
                if (res != 0) {
                    perror("Thread creation failed");
                    exit(EXIT_FAILURE);
                }
            }
            /*printf("Waiting for thread to finish...\n");*/
        }
    }else{
        printf("backward\n");
        //blockSz=(totPixels+ncpu-1)/(ncpu);
        ///*printf("blockSz=%d given ncpu=%d\n",blockSz,ncpu);*/
        //for(i=0; i<ncpu; i++){
        //    res = pthread_create(&a_thread[i], NULL, parPrjBack, (void *)i);
        //    if (res != 0) {
        //        perror("Thread creation failed");
        //        exit(EXIT_FAILURE);
        //    }
        //}
    }
    /*printf("Waiting for thread to finish...\n");*/
    for(i=0; i<ncpu; i++){
        res = pthread_join(a_thread[i], &thread_result);
        if (res != 0) {
            perror("Thread join failed");
            exit(EXIT_FAILURE);
        }
        /*printf("Thread joined, it returned %s\n", (char *)thread_result);*/
    }
    free(a_thread);
    //if(conf.fwd)
    //    for(i=0; i<conf.np*conf.prjWidth; i++)
    //        prj[i]/=beamWidth;
    //else
    //    for(i=0; i<totPixels; i++)
    //        img[i]/=beamWidth;
    return 0;
}

int main( void ) {
    int N = 1024/1;
    conf.nc=N;
    conf.nr=N;
    conf.prjWidth=N;
    conf.np=360;
    conf.prjFull=360;
    conf.dSize = 1;
    conf.effectiveRate=0.9;
    conf.d=80000;
    conf.fwd=1;
    conf.angleStep=1;
    conf.prjStep=N;
    conf.rayDrive=1;

    conf.imgSize=conf.nc*conf.nr;
    conf.sinoSize=conf.prjWidth*conf.np;

    CPUBitmap image( conf.nr, conf.nc );
    unsigned char *ptr = image.get_ptr();

    img = (ft*) malloc(conf.imgSize*sizeof(ft));
    int offset;
    int YC = conf.nr/2, XC = conf.nc/2;
    unsigned char value;
    for(int i=0; i < conf.nr; i++)
        for(int j=0; j < conf.nc; j++){
            offset = i*conf.nc+j;
            if( (i-YC-0)*(i-YC-0)+(j-XC)*(j-XC) < conf.nr*conf.nc/16){
                img[offset]=1;
            }else
                img[offset]=0;
//            if(i<5 && j < 5) img[i][j]=1;
            value = (unsigned char)(0xff * img[offset]);
            ptr[(offset<<2)+0] = value;
            ptr[(offset<<2)+1] = value;
            ptr[(offset<<2)+2] = value;
            ptr[(offset<<2)+3] = 0xff;
        }
    //image.display_and_exit();

    double elapsedTime;

    // start the timers
    timeval start;
    gettimeofday(&start,NULL);

    sino = (ft *) malloc(conf.sinoSize*sizeof(ft));

    memset(sino,0,conf.sinoSize*sizeof(ft));

    multiThread();

    timeval end;
    gettimeofday(&end,NULL);
    elapsedTime = (end.tv_sec-start.tv_sec)*1e3+(double(end.tv_usec-start.tv_usec))/1000;

    printf( "Time taken:  %f ms\n", elapsedTime );

    CPUBitmap sinogram( conf.prjWidth, conf.np );
    unsigned char *sinoPtr = sinogram.get_ptr();

    ft temp=0;
    for(int i=0; i < conf.np; i++)
        for(int j=0; j < conf.prjWidth; j++){
            offset = i*conf.prjWidth+j;
            if( sino[offset]>temp){
                //printf("%f",temp);
                temp=sino[offset];
                //printf(" -> %f\n",temp);
            }
        }
    for(int i=0; i < conf.np; i++)
        for(int j=0; j < conf.prjWidth; j++){
            offset = i*conf.prjWidth+j;
            value = (unsigned char)(255 * sino[offset]/temp);
            sinoPtr[(offset<<2)+0] = value;
            sinoPtr[(offset<<2)+1] = value;
            sinoPtr[(offset<<2)+2] = value;
            sinoPtr[(offset<<2)+3] = 0xff;
        }

    sinogram.display_and_exit();
    return 0;
}

