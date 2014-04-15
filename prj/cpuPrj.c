/*
 *   Description: use CPU to calculate parallel beam, fan beam projection 
 *   and their corresponding back projection.
 *
 *   Reference:
 *   Author: Renliang Gu (renliang@iastate.edu)
 *   $Revision: 0.1 $ $Date: Sun 10 Nov 2013 01:23:32 AM CST
 *
 *   v_0.1:     first draft
 */
#define CPU 1

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "./common/kiss_fft.h"

#ifdef __cplusplus
extern "C"{
#endif
#include "prj.h"
#ifdef __cplusplus
}
#endif

#if SHOWIMG
#include "./common/cpu_bitmap.h"
#endif

struct prjConf config;
struct prjConf* pConf = &config;

static unsigned int nthread=64;
static int fSize, bSize;
static ft *pImg, *pSino;
void (*rayDrive)(ft*, ft*, int);
void (*pixelDrive)(ft*, ft*, int);

void pixelDrivePar(ft* img, ft* sino, int threadIdx){
    //printf("entering pixelDrive\n");
    // detector array is of odd size with the center at the middle
    ft theta;       // the current angle, range in [0 45]
    int thetaIdx; // index counts from 0
    struct prjConf* conf = pConf;

    int pC = conf->prjWidth/2;
    int N=conf->n;

    int x=0, y=threadIdx;
    while(y>x){ x++; y-=x; }
    // make sure the border is set to be zero when size of image is odd
    if(x==N/2 && N%2==0){
        for(int i=0; i<N; i++){
            img[i]=0; img[i*N]=0;
        }
        return;
    }

    ft cosT, sinT;  // cosine and sine of theta

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft oc;
    ft beamWidth = conf->dSize*conf->effectiveRate;

    ft dist, weight;
    int temp, imgIdx;
    ft imgt[8];
    ft t, tl, tr;
    int dt, dtl, dtr;
    for(int i=0; i<8; i++) imgt[i]=0;
    for(thetaIdx=0; thetaIdx<conf->prjFull/2;thetaIdx++){

        theta  = thetaIdx *2*PI/conf->prjFull;
        cosT = cos(theta ); sinT = sin(theta );

        // up left
        oc = (x-0.5)*cosT + (y-0.5)*sinT;
        tl = oc; tr = tl;

        // up right
        //qe = (x+0.5)*sinT - (y-0.5)*cosT + d;
        //oc = (x+0.5)*cosT + (y-0.5)*sinT;
        oc = oc + cosT; t = oc;
        tl = MIN(tl, t); tr=MAX(tr,t);

        // bottom right
        //qe = (x+0.5)*sinT - (y+0.5)*cosT + d;
        //oc = (x+0.5)*cosT + (y+0.5)*sinT;
        oc=oc+sinT; t = oc;
        tl = MIN(tl, t); tr=MAX(tr,t);

        // bottom left
        //qe = (x-0.5)*sinT - (y+0.5)*cosT + d;
        //oc = (x-0.5)*cosT + (y+0.5)*sinT;
        oc = oc-cosT; t = oc;
        tl = MIN(tl, t); tr=MAX(tr,t);

        //qe = d+x*sinT-y*cosT; // positive direction is qo
        dtl = MAX((int)round((tl)/conf->dSize),-(conf->prjWidth-1)/2);
        dtr = MIN((int)round((tr)/conf->dSize), (conf->prjWidth-1)/2);

        for(dt=dtl; dt<=dtr; dt++){
            t = dt*conf->dSize;
            dist=x*cosT+y*sinT-t;
            weight=getWeight(dist,beamWidth,cosT,sinT);

            if(thetaIdx<conf->np){
                imgt[0] += weight*sino[thetaIdx*conf->prjWidth+dt+pC];
                imgt[2] += weight*sino[thetaIdx*conf->prjWidth-dt+pC];
            }

            temp=(thetaIdx+conf->prjFull/4)%conf->prjFull;
            if(temp<conf->np){
                imgt[1] += weight*sino[temp*conf->prjWidth+dt+pC];
                imgt[3] += weight*sino[temp*conf->prjWidth-dt+pC];
            }

            temp=(thetaIdx+conf->prjFull/2)%conf->prjFull;
            if(temp<conf->np){
                imgt[2] += weight*sino[temp*conf->prjWidth+dt+pC];
                imgt[0] += weight*sino[temp*conf->prjWidth-dt+pC];
            }

            temp=(thetaIdx+3*conf->prjFull/4)%conf->prjFull;
            if(temp<conf->np){
                imgt[3] += weight*sino[temp*conf->prjWidth+dt+pC];
                imgt[1] += weight*sino[temp*conf->prjWidth-dt+pC];
            }

            temp=(3*conf->prjFull/2-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[4] += weight*sino[temp*conf->prjWidth-dt+pC];
                imgt[6] += weight*sino[temp*conf->prjWidth+dt+pC];
            }

            temp=(7*conf->prjFull/4-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[5] += weight*sino[temp*conf->prjWidth-dt+pC];
                imgt[7] += weight*sino[temp*conf->prjWidth+dt+pC];
            }

            temp=(conf->prjFull-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[6] += weight*sino[temp*conf->prjWidth-dt+pC];
                imgt[4] += weight*sino[temp*conf->prjWidth+dt+pC];
            }

            temp=(5*conf->prjFull/4-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[7] += weight*sino[temp*conf->prjWidth-dt+pC];
                imgt[5] += weight*sino[temp*conf->prjWidth+dt+pC];
            }
        }
    }
    if(conf->cmd & FBP_BIT) for(int i=0; i<8; i++) imgt[i]=imgt[i]*PI/conf->np;
    if(x>(N-1)/2 || y>(N-1)/2) return;
    imgIdx = ( y+N/2)*N+x+N/2; img[imgIdx] = imgt[0]/conf->effectiveRate;
    imgIdx = ( x+N/2)*N-y+N/2; img[imgIdx] = imgt[1]/conf->effectiveRate;
    imgIdx = (-y+N/2)*N-x+N/2; img[imgIdx] = imgt[2]/conf->effectiveRate;
    imgIdx = (-x+N/2)*N+y+N/2; img[imgIdx] = imgt[3]/conf->effectiveRate;
    if(y==0 || y>=x) return;
    imgIdx = (-y+N/2)*N+x+N/2; img[imgIdx] = imgt[4]/conf->effectiveRate;
    imgIdx = ( x+N/2)*N+y+N/2; img[imgIdx] = imgt[5]/conf->effectiveRate;
    imgIdx = ( y+N/2)*N-x+N/2; img[imgIdx] = imgt[6]/conf->effectiveRate;
    imgIdx = (-x+N/2)*N-y+N/2; img[imgIdx] = imgt[7]/conf->effectiveRate;
}

void pixelDriveFan(ft* img, ft* sino, int threadIdx){
    //printf("entering pixelDrive\n");
    // detector array is of odd size with the center at the middle
    ft theta;       // the current angle, range in [0 45]
    int thetaIdx; // index counts from 0
    struct prjConf* conf = pConf;

    int pC = conf->prjWidth/2;
    ft d=conf->d; // the distance from rotating center to source in pixel
    int N=conf->n;

    int x=0, y=threadIdx;
    while(y>x){ x++; y-=x; }

    ft cosT, sinT;  // cosine and sine of theta

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft qe, oc, qa;
    ft bq;
    ft cosB, sinB; //, tanB=t/d;
    ft cosR, sinR;
    ft beamWidth = conf->dSize*conf->effectiveRate;
    ft bw;

    ft dist, weight;
    int temp, imgIdx;
    ft imgt[8];
    ft t, tl, tr;
    int dt, dtl, dtr;
    for(int i=0; i<8; i++) imgt[i]=0;
    for(thetaIdx=0; thetaIdx<conf->prjFull;thetaIdx++){

        theta  = thetaIdx *2*PI/conf->prjFull;
        cosT = cos(theta ); sinT = sin(theta );

        // up letf
        qe = (x-0.5)*sinT - (y-0.5)*cosT + d;
        oc = (x-0.5)*cosT + (y-0.5)*sinT;
        tl = d*oc/qe; tr = tl;

        // up right
        //qe = (x+0.5)*sinT - (y-0.5)*cosT + d;
        //oc = (x+0.5)*cosT + (y-0.5)*sinT;
        qe = qe+sinT; oc = oc + cosT; t = d*oc/qe;
        tl = MIN(tl, t); tr=MAX(tr,t);

        // bottom right
        //qe = (x+0.5)*sinT - (y+0.5)*cosT + d;
        //oc = (x+0.5)*cosT + (y+0.5)*sinT;
        qe = qe-cosT; oc=oc+sinT; t = d*oc/qe;
        tl = MIN(tl, t); tr=MAX(tr,t);

        // bottom left
        //qe = (x-0.5)*sinT - (y+0.5)*cosT + d;
        //oc = (x-0.5)*cosT + (y+0.5)*sinT;
        qe = qe-sinT; oc = oc-cosT; t = d*oc/qe;
        tl = MIN(tl, t); tr=MAX(tr,t);

        //qe = d+x*sinT-y*cosT; // positive direction is qo
        qe = qe + sinT/2 +cosT/2;
        dtl = MAX((int)round((tl)/conf->dSize),-(conf->prjWidth-1)/2);
        dtr = MIN((int)round((tr)/conf->dSize), (conf->prjWidth-1)/2);
        qa = sqrt(2*d*qe-d*d+x*x+y*y);

        for(dt=dtl; dt<=dtr; dt++){
            t = dt*conf->dSize;
            bq=sqrt(d*d+t*t);
            cosB=d/bq; sinB=t/bq; //, tanB=t/d;
            cosR=cosB*cosT-sinB*sinT;
            sinR=sinB*cosT+cosB*sinT;
            dist=x*cosR+y*sinR-d*t/bq;
            bw = qe*beamWidth*cosB/d;
            weight=getWeight(dist,bw,cosR,sinR);

            // method provide by the books
            //if(conf->cmd & FBP_BIT) weight = weight*d*d/qe/qe;
            // The one I think should be
            if(conf->cmd & FBP_BIT) weight = weight*d/qa;
            //if(conf->cmd & FBP_BIT) weight = weight*d/qe;

            if(thetaIdx<conf->np)
                imgt[0] += weight*sino[thetaIdx*conf->prjWidth+dt+pC];

            temp=(thetaIdx+conf->prjFull/4)%conf->prjFull;
            if(temp<conf->np){
                imgt[1] += weight*sino[temp*conf->prjWidth+dt+pC];
            }

            temp=(thetaIdx+conf->prjFull/2)%conf->prjFull;
            if(temp<conf->np){
                imgt[2] += weight*sino[temp*conf->prjWidth+dt+pC];
            }

            temp=(thetaIdx+3*conf->prjFull/4)%conf->prjFull;
            if(temp<conf->np){
                imgt[3] += weight*sino[temp*conf->prjWidth+dt+pC];
            }

            temp=(3*conf->prjFull/2-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[4] += weight*sino[temp*conf->prjWidth-dt+pC];
            }

            temp=(7*conf->prjFull/4-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[5] += weight*sino[temp*conf->prjWidth-dt+pC];
            }

            temp=(conf->prjFull-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[6] += weight*sino[temp*conf->prjWidth-dt+pC];
            }

            temp=(5*conf->prjFull/4-thetaIdx)%conf->prjFull;
            if(temp<conf->np){
                imgt[7] += weight*sino[temp*conf->prjWidth-dt+pC];
            }
        }
    }
    if(conf->cmd & FBP_BIT) for(int i=0; i<8; i++) imgt[i]=imgt[i]*PI/conf->np;
    if(x>(N-1)/2 || y>(N-1)/2) return;
    imgIdx = ( y+N/2)*N+x+N/2; img[imgIdx] = imgt[0]/conf->effectiveRate;
    imgIdx = ( x+N/2)*N-y+N/2; img[imgIdx] = imgt[1]/conf->effectiveRate;
    imgIdx = (-y+N/2)*N-x+N/2; img[imgIdx] = imgt[2]/conf->effectiveRate;
    imgIdx = (-x+N/2)*N+y+N/2; img[imgIdx] = imgt[3]/conf->effectiveRate;
    if(y==0 || y>=x) return;
    imgIdx = (-y+N/2)*N+x+N/2; img[imgIdx] = imgt[4]/conf->effectiveRate;
    imgIdx = ( x+N/2)*N+y+N/2; img[imgIdx] = imgt[5]/conf->effectiveRate;
    imgIdx = ( y+N/2)*N-x+N/2; img[imgIdx] = imgt[6]/conf->effectiveRate;
    imgIdx = (-x+N/2)*N-y+N/2; img[imgIdx] = imgt[7]/conf->effectiveRate;
}

void rayDrivePar(ft* img, ft* sino, int threadIdx){
    //printf("entering rayDrive\n");
    // detector array is of odd size with the center at the middle
    int temp;
    const struct prjConf * conf = pConf;
    int sinoIdx;
    int thetaIdx, tIdx;
    if(conf->prjWidth%2==0){
        if(threadIdx==fSize-1){
            for(int i=0; i<conf->np; i++)
                sino[i*conf->prjWidth]=0;
            return;
        }
        thetaIdx= threadIdx/(pConf->prjWidth/2);
        tIdx = threadIdx%(pConf->prjWidth/2)+1;
    }else{
        thetaIdx = threadIdx/((pConf->prjWidth+1)/2);
        tIdx = threadIdx%((pConf->prjWidth+1)/2);
    }

    ft theta;       // the current angle, range in [0 45]
    ft t;           // position of current detector

    int N =conf->n;// N is of size NxN centering at (N/2,N/2)
    int pC = conf->prjWidth/2;
    ft beamWidth = conf->dSize*conf->effectiveRate;
    ft hbeamW = beamWidth/2;

    theta = thetaIdx*2*PI/conf->prjFull;
    t = (tIdx-pC)*conf->dSize;
    ft tl = t-hbeamW,
       tr = t+hbeamW;

    ft cosT, sinT;  // cosine and sine of theta
    cosT=cos(theta); sinT=sin(theta);

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);

    ft   xl, xr;
    int dxl,dxr;

    int x,y=(N-1)/2;
    //xc = xc + y*slopeXYc;

    // beamw is based on the position of this pixel
    ft dist, weight;
    ft sinot[8];
    for(int i=0; i<8; i++) sinot[i]=0;
    for(y=(N-1)/2; y>=-(N-1)/2; y--){
        xl = (tl-sinT *(y+0.5f))/cosT;
        xr = (tr-sinT *(y-0.5f))/cosT;

        dxl =MAX((int)round(xl ),-(N-1)/2);
        dxr =MIN((int)round(xr ), (N-1)/2);

        //if(thetaIdx==45 && blockIdx.y==1 && threadIdx.x==0){
        //    printf("\nthetaIdx=%d, t=%f, y=%d, %d->%d\n \t",
        //            thetaIdx,t,y,dxll,dxrr);
        //    for(temp=0; temp<=dxrr-dxll; temp++){
        //        if(shared[1][temp]>0)
        //            printf("%d: %d: %f",temp,temp+dxll,shared[1][temp]);
        //    }
        //}
        for(x=dxl; x<=dxr; x++){
            dist=x*cosT+y*sinT-t;
            weight=getWeight(dist,beamWidth,cosT,sinT);

            sinot[0]+=weight*img[( y+N/2)*N+x+N/2];
            //temp=thetaIdx+conf->prjFull/4;
            sinot[1]+=weight*img[( x+N/2)*N-y+N/2]; //img[imgIdx];
            //temp+=conf->prjFull/4;
            sinot[2]+=weight*img[(-y+N/2)*N-x+N/2]; //img[imgIdx];
            //temp+=conf->prjFull/4;
            sinot[3]+=weight*img[(-x+N/2)*N+y+N/2]; //img[imgIdx];
            //temp=conf->prjFull/2-thetaIdx;
            sinot[4]+=weight*img[(-y+N/2)*N+x+N/2]; //img[imgIdx];
            //temp=3*conf->prjFull/4-thetaIdx;
            sinot[5]+=weight*img[( x+N/2)*N+y+N/2]; //img[imgIdx];
            //temp=conf->prjFull-thetaIdx;
            sinot[6]+=weight*img[( y+N/2)*N-x+N/2]; //img[imgIdx];
            //temp=conf->prjFull/4-thetaIdx;
            sinot[7]+=weight*img[(-x+N/2)*N-y+N/2]; //img[imgIdx];
        }
    }

    if(tIdx>pC) return;

    if(thetaIdx<conf->np){
        sinoIdx=thetaIdx*conf->prjWidth;
        sino[sinoIdx+tIdx]=sinot[0]/conf->effectiveRate;
        sino[sinoIdx+2*pC-tIdx]=sinot[2]/conf->effectiveRate;
    }

    temp = thetaIdx+conf->prjFull/4;
    if(temp<conf->np){
        sinoIdx = temp*conf->prjWidth;
        sino[sinoIdx+tIdx]=sinot[1]/conf->effectiveRate;
        sino[sinoIdx+2*pC-tIdx]=sinot[3]/conf->effectiveRate;
    }

    temp = thetaIdx+conf->prjFull/2;
    if(temp<conf->np){
        sinoIdx = temp*conf->prjWidth;
        sino[sinoIdx+tIdx]=sinot[2]/conf->effectiveRate;
        sino[sinoIdx+2*pC-tIdx]=sinot[0]/conf->effectiveRate;
    }

    temp = thetaIdx+3*conf->prjFull/4;
    if(temp<conf->np){
        sinoIdx = temp*conf->prjWidth;
        sino[sinoIdx+tIdx]=sinot[3]/conf->effectiveRate;
        sino[sinoIdx+2*pC-tIdx]=sinot[1]/conf->effectiveRate;
    }

    if(thetaIdx>0 && thetaIdx<conf->prjFull*0.125f){
        tIdx = 2*pC-tIdx;

        temp = conf->prjFull/2-thetaIdx;
        if(temp<conf->np){
            sinoIdx = temp*conf->prjWidth;
            sino[sinoIdx+tIdx]=sinot[4]/conf->effectiveRate;
            sino[sinoIdx+2*pC-tIdx]=sinot[6]/conf->effectiveRate;
        }

        temp = 3*conf->prjFull/4-thetaIdx;
        if(temp<conf->np){
            sinoIdx = temp*conf->prjWidth;
            sino[sinoIdx+tIdx]=sinot[5]/conf->effectiveRate;
            sino[sinoIdx+2*pC-tIdx]=sinot[7]/conf->effectiveRate;
        }

        temp = conf->prjFull-thetaIdx;
        if(temp<conf->np){
            sinoIdx = temp*conf->prjWidth;
            sino[sinoIdx+tIdx]=sinot[6]/conf->effectiveRate;
            sino[sinoIdx+2*pC-tIdx]=sinot[4]/conf->effectiveRate;
        }

        temp = conf->prjFull/4-thetaIdx;
        if(temp<conf->np){
            sinoIdx = temp*conf->prjWidth;
            sino[sinoIdx+tIdx]=sinot[7]/conf->effectiveRate;
            sino[sinoIdx+2*pC-tIdx]=sinot[5]/conf->effectiveRate;
        }
    }
}

void rayDriveFan(ft* img, ft* sino, int threadIdx){
    //printf("entering rayDrive\n");
    // detector array is of odd size with the center at the middle
    const struct prjConf * conf = pConf;
    int sinoIdx;
    int thetaIdx, tIdx;
    if(conf->prjWidth%2==0){
        thetaIdx= threadIdx/(pConf->prjWidth-1);
        tIdx = threadIdx%(pConf->prjWidth-1)+1;
    }else{
        thetaIdx= threadIdx/pConf->prjWidth;
        tIdx = threadIdx%pConf->prjWidth;
    }
    //printf("tIdx=%d, thetaIdx=%d\n",tIdx, thetaIdx);

    //if(blockIdx.x==0 && blockIdx.y==0)
    //    printf("gridDim=(%d, %d)\t blockDim=(%d, %d)\n",
    //            gridDim.x, gridDim.y, blockDim.x, blockDim.y);

    ft theta;       // the current angle, range in [0 45]
    ft t;           // position of current detector
    ft d;           // the distance from rotating center to X-ray source in pixel

    int N =conf->n;// N is of size NxN centering at (N/2,N/2)
    int pC = conf->prjWidth/2;
    ft beamWidth = conf->dSize*conf->effectiveRate;
    ft hbeamW = beamWidth/2;

    theta = thetaIdx*2*PI/conf->prjFull;
    t = (tIdx-pC)*conf->dSize;
    ft tl = t-hbeamW,
       tr = t+hbeamW;
    d = conf->d;

    ft cosT, sinT;  // cosine and sine of theta
    cosT=cos(theta); sinT=sin(theta);

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft bq=sqrt(d*d+t*t);
    ft cosB=d/bq, sinB=t/bq; //, tanB=t/d;
    ft cosR=cosB*cosT-sinB*sinT; // cosR and sinR are based on t
    ft sinR=sinB*cosT+cosB*sinT;
    ft beamwidthCosB=beamWidth*cosB;

    //if(tIdx==pC && thetaIdx==30){
    //    printf("theta=%f,cosT=%f, sinT=%f\n",theta,cosT,sinT);
    //    printf("cosB=%f, sinB=%f\n",cosB,sinB);
    //}

    ft dtl=d*tl, dtr=d*tr;
    ft QxBxl = tl *cosT+d*sinT, QxBxr = tr *cosT+d*sinT;
    ft QyByl =-tl *sinT+d*cosT, QyByr =-tr *sinT+d*cosT;

    ft   xl, xr;
    int dxl,dxr;

    int x,y=(N-1)/2;
    //xc = xc + y*slopeXYc;

    // beamw is based on the position of this pixel
    ft bw, dist, weight;
    int temp;
    ft sinot[8];
    for(int i=0; i<8; i++) sinot[i]=0;
    for(y=(N-1)/2; y>=-(N-1)/2; y--){
        if(QxBxl>0) xl = (dtl - QxBxl *(y+0.5f))/QyByl;
        else xl = (dtl -QxBxl *(y-0.5f))/QyByl;
        if(QxBxr>0) xr = (dtr - QxBxr *(y-0.5f))/QyByr;
        else xr = (dtr - QxBxr *(y+0.5f))/QyByr;

        dxl =MAX((int)round(xl ),-(N-1)/2);
        dxr =MIN((int)round(xr ), (N-1)/2);

        //if(thetaIdx==45 && blockIdx.y==1 && threadIdx.x==0){
        //    printf("\nthetaIdx=%d, t=%f, y=%d, %d->%d\n \t",
        //            thetaIdx,t,y,dxll,dxrr);
        //    for(temp=0; temp<=dxrr-dxll; temp++){
        //        if(shared[1][temp]>0)
        //            printf("%d: %d: %f",temp,temp+dxll,shared[1][temp]);
        //    }
        //}
        for(x=dxl; x<=dxr; x++){

            dist=x*cosR+y*sinR-d*t/bq;
            bw=beamwidthCosB + (x*sinT-y*cosT)*beamwidthCosB/d;
            weight=getWeight(dist,bw,cosR,sinR);

            sinot[0]+=weight*img[( y+N/2)*N+x+N/2];

            //if(thetaIdx==42 && blockIdx.y==4 && threadIdx.x==0){
            //    printf("%d %d %e %e %e\n",y,x,weight, weight*shared[0][x-dxll],sinot[0]);
            //}

            //temp=thetaIdx+conf->prjFull/4;
            sinot[1]+=weight*img[( x+N/2)*N-y+N/2]; //img[imgIdx];
            //temp+=conf->prjFull/4;
            sinot[2]+=weight*img[(-y+N/2)*N-x+N/2]; //img[imgIdx];
            //temp+=conf->prjFull/4;
            sinot[3]+=weight*img[(-x+N/2)*N+y+N/2]; //img[imgIdx];
            //temp=conf->prjFull/2-thetaIdx;
            sinot[4]+=weight*img[(-y+N/2)*N+x+N/2]; //img[imgIdx];
            //temp=3*conf->prjFull/4-thetaIdx;
            sinot[5]+=weight*img[( x+N/2)*N+y+N/2]; //img[imgIdx];
            //temp=conf->prjFull-thetaIdx;
            sinot[6]+=weight*img[( y+N/2)*N-x+N/2]; //img[imgIdx];
            //temp=conf->prjFull/4-thetaIdx;
            sinot[7]+=weight*img[(-x+N/2)*N-y+N/2]; //img[imgIdx];
        }
    }

    if(tIdx>=conf->prjWidth) return;

    if(thetaIdx<conf->np){
        sinoIdx=thetaIdx*conf->prjWidth+tIdx;
        sino[sinoIdx]=sinot[0]/conf->effectiveRate;
    }

    temp = thetaIdx+conf->prjFull/4;
    if(temp<conf->np){
        sino[temp*conf->prjWidth+tIdx]=sinot[1]/conf->effectiveRate;
    }

    temp = thetaIdx+conf->prjFull/2;
    if(temp<conf->np){
        sino[temp*conf->prjWidth+tIdx]=sinot[2]/conf->effectiveRate;
    }

    temp = thetaIdx+3*conf->prjFull/4;
    if(temp<conf->np){
        sino[temp*conf->prjWidth+tIdx]=sinot[3]/conf->effectiveRate;
    }

    if(thetaIdx>0 && thetaIdx<conf->prjFull*0.125f){
        tIdx = 2*pC-tIdx;

        temp = conf->prjFull/2-thetaIdx;
        if(temp<conf->np)
            sino[temp*conf->prjWidth+tIdx]=sinot[4]/conf->effectiveRate;

        temp = 3*conf->prjFull/4-thetaIdx;
        if(temp<conf->np)
            sino[temp*conf->prjWidth+tIdx]=sinot[5]/conf->effectiveRate;

        temp = conf->prjFull-thetaIdx;
        if(temp<conf->np)
            sino[temp*conf->prjWidth+tIdx]=sinot[6]/conf->effectiveRate;

        temp = conf->prjFull/4-thetaIdx;
        if(temp<conf->np)
            sino[temp*conf->prjWidth+tIdx]=sinot[7]/conf->effectiveRate;
    }
}

void setup(int n, int prjWidth, int np, int prjFull, ft dSize, ft 
        effectiveRate, ft d){
    config.n=n; config.prjWidth=prjWidth;
    config.np=np; config.prjFull=prjFull;
    config.dSize=dSize; config.effectiveRate=effectiveRate;
    config.d=d;

    config.imgSize=config.n*config.n;
    config.sinoSize=config.prjWidth*config.np;

    if(config.d>0){
        rayDrive = rayDriveFan;
        pixelDrive = pixelDriveFan;

        fSize = MIN(pConf->np, pConf->prjFull/8+1);
        if(pConf->prjWidth%2==0)
            fSize = fSize*(pConf->prjWidth-1);
        else
            fSize = fSize*pConf->prjWidth;
        bSize = (pConf->n+1)/2;
        bSize = (1+bSize)*bSize/2;
    }else{
        rayDrive = rayDrivePar;
        pixelDrive = pixelDrivePar;

        fSize = MIN(pConf->np, pConf->prjFull/8+1);
        fSize = fSize*((pConf->prjWidth+1)/2);
        bSize = (pConf->n+1)/2;
        bSize = (1+bSize)*bSize/2;
        if(pConf->n%2==0) bSize++;
    }

#if DEBUG
    printf("fsize=%d, bsize=%d, nthread=%d\n", fSize,bSize,nthread);
    printf("setup done\n");
#endif
}

void showSetup(){
    printf("config.n=%d\n",config.n);
    printf("config.prjWidth=%d\n",config.prjWidth);
    printf("config.np=%d\n",config.np);
    printf("config.prjFull=%d\n",config.prjFull);
    printf("config.dSize=%g\n",config.dSize);
    printf("config.effectiveRate=%g\n",config.effectiveRate);
    printf("config.d=%g\n",config.d);
}

void *parPrjFor(void *arg){
    for(int i=(size_t)arg; i<fSize; i+=nthread){
        rayDrive(pImg,pSino,i);
    }
    return 0;
}

void *parPrjBack(void *arg){
    for(int i=(size_t)arg; i<bSize; i+=nthread){
        pixelDrive(pImg,pSino,i);
    }
    return 0;
}

int cpuPrj(ft* img, ft* sino, char cmd){
#if EXE_PROF
    // add some instruction for cpu execution profile
#endif
    pImg = img; pSino = sino;

#if EXE_TIME
    // start the timers
    clock_t start, end;
    start = clock();
#endif

    int res;
    pthread_t *a_thread;
    void *thread_result;

    a_thread=(pthread_t*)calloc(nthread,sizeof(pthread_t));
    if(cmd & FWD_BIT){
#if DEBUG
        printf("Forward projecting ...\n");
#endif
        for(size_t i=0; i<nthread; i++){
            res = pthread_create(&a_thread[i], NULL, parPrjFor, (void *)i);
            if (res != 0) {
                perror("Thread creation failed");
                exit(EXIT_FAILURE);
            }
        }
        /*printf("Waiting for thread to finish...\n");*/
    }else if(cmd & BWD_BIT){
#if DEBUG
        printf("Backward projecting ...\n");
#endif
        for(size_t i=0; i<nthread; i++){
            res = pthread_create(&a_thread[i], NULL, parPrjBack, (void *)i);
            if (res != 0) {
                perror("Thread creation failed");
                exit(EXIT_FAILURE);
            }
        }
    }else if(cmd & FBP_BIT){
#if DEBUG
        printf("Filtered Backprojecting ...\n");
#endif
        pConf->cmd = FBP_BIT; //when put this don't foget to remove it later
        pSino = (ft*) calloc(pConf->sinoSize,sizeof(ft));
        ft bq;
        int pC = pConf->prjWidth/2;

#if DEBUG
        FILE* f = fopen("sinogram_0.data","wb");
        fwrite(sino, sizeof(ft), config.sinoSize, f);
        fclose(f);
#endif

        if(pConf->d>0){
#if DEBUG
            printf("reconstructing by FBP (fan beam) ... \n");
#endif
            for(int j=0; j<(pConf->prjWidth+1)/2; j++){
                bq = sqrt(pConf->d*pConf->d + j*j*pConf->dSize*pConf->dSize);
                for(int i=0, idx1=pC-j, idx2=pC+j; i<pConf->np;
                        i++, idx1+=pConf->prjWidth, idx2+=pConf->prjWidth){
                    pSino[idx1]=sino[idx1]*pConf->d / bq;
                    pSino[idx2]=sino[idx2]*pConf->d / bq;
                }
            }
            if(pConf->prjWidth%2==0){
                bq = sqrt(pConf->d*pConf->d + pC*pC*pConf->dSize*pConf->dSize);
                for(int i=0, idx1=0; i<pConf->np;
                        i++, idx1+=pConf->prjWidth){
                    pSino[idx1]=sino[idx1]*pConf->d / bq;
                }
            }
        }else{
#if DEBUG
            printf("reconstructing by FBP (parallel beam) ... \n");
#endif
            memcpy(pSino,sino,pConf->sinoSize*sizeof(ft));
        }

#if DEBUG
        f = fopen("sinogram_1.data","wb");
        fwrite(pSino, sizeof(ft), config.sinoSize, f);
        fclose(f);
#endif

        for(int i=0; i<pConf->np; i++){
            rampFilter(pSino+i*pConf->prjWidth, pConf->prjWidth, pConf->dSize);
        }

#if DEBUG
        f = fopen("sinogram_2.data","wb");
        fwrite(pSino, sizeof(ft), config.sinoSize, f);
        fclose(f);
#endif
        
        for(size_t i=0; i<nthread; i++){
            res = pthread_create(&a_thread[i], NULL, parPrjBack, (void *)i);
            if (res != 0) {
                perror("Thread creation failed");
                exit(EXIT_FAILURE);
            }
        }

    }
    /*printf("Waiting for thread to finish...\n");*/
    for(size_t i=0; i<nthread; i++){
        res = pthread_join(a_thread[i], &thread_result);
        if (res != 0) {
            perror("Thread join failed");
            exit(EXIT_FAILURE);
        }
        /*printf("Thread joined, it returned %s\n", (char *)thread_result);*/
    }
    pConf->cmd = 0;
    if(pSino!=sino) free(pSino);
    free(a_thread);

#if EXE_PROF
    // add cpu execution profile instructions
#endif

#if EXE_TIME
     double elapsedTime;
     end = clock();
     elapsedTime = 1000*((double) (end - start)) / CLOCKS_PER_SEC;
     printf( "Time taken:  %3.1f ms\n", elapsedTime );
#endif
    //FILE* f = fopen("sinogram.data","wb");
    //fwrite(sino, sizeof(ft), pConf->sinoSize, f);
    //fclose(f);
    return 0;
}

int forwardTest( void ) {
#if SHOWIMG
    CPUBitmap image( config.n, config.n );
    unsigned char *ptr = image.get_ptr();
    CPUBitmap sinogram( config.prjWidth, config.np );
    unsigned char *sinoPtr = sinogram.get_ptr();
    unsigned char value;
#endif

    ft* img = (ft*) malloc(config.imgSize*sizeof(ft));
    ft *sino = (ft *) malloc(config.sinoSize*sizeof(ft));
    int offset;
    int YC = config.n/2, XC = config.n/2;
    for(int i=0; i < config.n; i++){
        for(int j=0; j < config.n; j++){
            offset = i*config.n+j;
            if(((i-YC-0.02*config.n)*(i-YC-0.02*config.n)+(j-XC-0)*(j-XC-0)<=(0.32*config.n)*(0.32*config.n))
                    && ((i-YC-0.12*config.n)*(i-YC-0.12*config.n)+
                        (j-XC-0.12*config.n)*(j-XC-0.12*config.n)>=(0.088*config.n)*(0.088*config.n))
              ){
                img[offset]=1;
            }else
                img[offset]=0;
            //            if(i<5 && j < 5) img[i][j]=1;
#if SHOWIMG
            value = (unsigned char)(0xff * img[offset]);
            offset = (config.n-1-i)*config.n+j;
            ptr[(offset<<2)+0] = value;
            ptr[(offset<<2)+1] = value;
            ptr[(offset<<2)+2] = value;
            ptr[(offset<<2)+3] = 0xff;
#endif
        }
    }
#if DEBUG
    FILE* f = fopen("img.data","w");
    fwrite(img,sizeof(ft), pConf->imgSize,f);
    fclose(f);
#endif
    //image.display_and_exit();

    cpuPrj(img, sino, RENEW_MEM | FWD_BIT);

#if DEBUG
    f = fopen("sinogram.data","wb");
    fwrite(sino, sizeof(ft), config.sinoSize, f);
    fclose(f);
#endif
#if SHOWIMG
    ft temp=0;
    for(int i=0; i < config.np; i++){
        for(int j=0; j < config.prjWidth; j++){
            offset = i*config.prjWidth+j;
            if( sino[offset]>temp)
                temp=sino[offset];
        }
    }
    for(int i=0; i < config.np; i++)
        for(int j=0; j < config.prjWidth; j++){
            offset = i*config.prjWidth+j;
            value = (unsigned char)(255 * sino[offset]/temp);
            offset = (config.np-1-i)*config.prjWidth+j;
            sinoPtr[(offset<<2)+0] = value;
            sinoPtr[(offset<<2)+1] = value;
            sinoPtr[(offset<<2)+2] = value;
            sinoPtr[(offset<<2)+3] = 0xff;
        }
#endif
    free(img); free(sino);

#if SHOWIMG
    //sinogram.display_and_exit();
#endif
    return 0;
}

int backwardTest( void ) {
#if SHOWIMG
    CPUBitmap image( config.n, config.n );
    unsigned char *ptr = image.get_ptr();
    CPUBitmap sinogram( config.prjWidth, config.np );
    unsigned char *sinoPtr = sinogram.get_ptr();
    unsigned char value;
    int offset;
#endif

    ft* img = (ft*) malloc(config.imgSize*sizeof(ft));
    ft *sino = (ft *) malloc(config.sinoSize*sizeof(ft));

#if DEBUG
    FILE* f = fopen("sinogram.data","rb");
    if(!fread(sino,sizeof(ft),config.sinoSize,f))
        perror("cannot read from sinogram.data\n");
    fclose(f);
#endif

    /*  
    int tempI;
    tempI=rand()%config.np;
    for(int i=0; i < config.np; i++){
        for(int j=0; j < config.prjWidth; j++){
            offset = i*config.prjWidth+j;
            if(i==tempI*0 && abs(j-config.prjWidth/2)<=150){
                if(offset%15<6) sino[offset]=1;
                else sino[offset]=0;
            }else
                sino[offset]=0;
            //fscanf(f,"%f ", &sino[offset]);
            //            if(i<5 && j < 5) img[i][j]=1;
#if SHOWIMG
            value = (unsigned char)(255 * sino[offset]);
            offset = (config.np-1-i)*config.prjWidth+j;
            sinoPtr[(offset<<2)+0] = value;
            sinoPtr[(offset<<2)+1] = value;
            sinoPtr[(offset<<2)+2] = value;
            sinoPtr[(offset<<2)+3] = 0xff;
#endif
        }
        //fscanf(f,"\n");
    }*/
    //sinogram.display_and_exit();

    cpuPrj(img,sino, RENEW_MEM | BWD_BIT);

#if DEBUG
    f = fopen("reImg.data","wb");
    fwrite(img,sizeof(ft),config.imgSize,f);
    fclose(f);
#endif

#if SHOWIMG
    ft temp=0;
    for(int i=0; i < config.n; i++)
        for(int j=0; j < config.n; j++){
            offset = i*config.n+j;
            if( img[offset]>temp) temp=img[offset];
        }
    //printf("temp=%g\n",temp);
    for(int i=0; i < config.n; i++){
        for(int j=0; j < config.n; j++){
            offset = i*config.n+j;
            if(img[offset]>0) value = (unsigned char)(255 * img[offset]/temp);
            else value = 0;
            offset = (config.n-1-i)*config.n+j;
            ptr[(offset<<2)+0] = value;
            ptr[(offset<<2)+1] = value;
            ptr[(offset<<2)+2] = value;
            ptr[(offset<<2)+3] = 0xff;
        }
    }
#endif

    free(sino); free(img);
#if SHOWIMG
    image.display_and_exit();
#endif
    return 0;
}

int main(int argc, char *argv[]){
    int N = 512;
    if(argc>1){
        N = atoi(argv[1]);
        printf("N=%d\n",N);
    }

    setup(N,N,720,720,1,1,3*N);
    showSetup();
    forwardTest();
    backwardTest();
    return 0;

    setup(N,N,180,360,1,1,0);
    showSetup();
    forwardTest();
    backwardTest();
    return 0;

    ft signal[1024];
    for(int i=0; i< 1024; i++){
        if(i>=300 && i<600)
            signal[i]=1.0f;
        else
            signal[i]=0;
    }
    rampFilter(signal, 1024, 1);
    return 0;

}

