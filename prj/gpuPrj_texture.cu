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
 *         Author: Renliang Gu (gurenliang@gmail.com)
 *   Organization:  Iowa State University
 *
 * ====================================================================
 */

/*#include <unistd.h>*/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cuda_profiler_api.h>
#include <pthread.h>
#include "gpuPrj.h"
#include "common/cpu_bitmap.h"
#include "common/book.h"
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
ft *img, *prj, *maskIdx;
prjConf *opt;
int blockSz, angleSz;

texture<ft,2> tex_img;
texture<ft,2> tex_sino;

/*
   */
__device__ ft getWeight(ft dist, ft beamWidth, ft cosR, ft sinR){
    // range of gamma is in [-alpha,pi/4+alpha], where alpha is small
    //ft temp=abs(fmod(theta,90)-45)*PI/180; /* test value of temp */
    ft height=1/MAX(cosR,sinR);
    ft bias1, bias2;             //bias2 should be larger.
    ft temp1=abs(cosR+sinR)/2;
    ft temp2=abs(cosR-sinR)/2;
    if(temp1>temp2){
        bias2=temp1; bias1=temp2;
    }else{
        bias1=temp1; bias2=temp2;
    }

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

__global__ void pixelDrive(ft* img, ft* sino, prjConf* conf){
}

__global__ void rayDrive(ft* img, ft* sino, prjConf* conf){
    // detector array is of odd size with the center at the middle
    int tIdx, thetaIdx; // index counts from 0
    int sinoIdx;
    //tIdx = blockIdx.y;
    //thetaIdx = blockIdx.x;
    int idx = (blockIdx.x*gridDim.y+blockIdx.y)*blockDim.x+threadIdx.x;
    tIdx = idx % conf->prjWidth;
    thetaIdx = idx / conf->prjWidth;

    //if(blockIdx.x==0 && threadIdx.x==0)
    //    printf("gridDim=(%d, %d)\t blockDim=(%d, %d)\n",
    //            gridDim.x, gridDim.y, blockDim.x, blockDim.y);

    if(tIdx==0 && conf->prjWidth%2==0){
        //printf("thetaIdx=%d\n",thetaIdx);
        //printf("(x,y)=(%d,%d)\n",gridDim.x,gridDim.y);
        sinoIdx=thetaIdx*conf->prjWidth+tIdx;
        sino[sinoIdx]=0;
        return; //make sure the size is odd
    }

    ft theta;       // the current angle, range in [0 45]
    ft t;           // position of current detector
    ft d;           // the distance from rotating center to X-ray source in pixel
    ft dSize;       // size of each detector in pixel

    int N;   // N is even// size of image is NxN, and centers at (N/2,N/2)
    int pC = conf->prjWidth/2;

    if(conf->nc==conf->nr) N=conf->nc;
    else
        N=0;
    theta = thetaIdx*2*PI/conf->prjFull;
    dSize = conf->dSize;
    t = (tIdx-pC)*dSize;
    d = conf->d;

    ft cosT, sinT;  // cosine and sine of theta
    cosT=cos(theta); sinT=sin(theta);

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft bq=sqrt(d*d+t*t);
    ft cosB=d/bq, sinB=t/bq; //, tanB=t/d;
    ft cosR=cosB*cosT-sinB*sinT;
    ft sinR=sinB*cosT+cosB*sinT;
    ft beamWidth = dSize*conf->effectiveRate;
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
    ft imgVal;

    while(y>-N/2){
        for(x=dxl,
                bwPre = (x*sinT-y*cosT),
                bw = beamwidthCosB+bwPre*beamwidthCosBOverD;
                x<=dxr;
                x++,
                bw+=sinTBeamwidthCosBOverD){

            dist=x*cosR+y*sinR-d*t/bq;
            weight=getWeight(dist,bw,cosR,sinR);
            if(tIdx==pC && thetaIdx==30){
                //printf("(%d,%d)\n",x,y);
                //printf("%f,%f;\n",dist,weight);
            }

            sinoIdx=thetaIdx*conf->prjWidth+tIdx;
            imgVal = tex2D(tex_img,y+N/2,x+N/2);
            imgIdx = (y+N/2)*N+x+N/2;
            //sino[sinoIdx]+=weight*imgVal;
            sino[sinoIdx]+=weight*img[imgIdx];
            
            //temp=thetaIdx+conf->prjFull/4;
            //if(temp<conf->np){
            //    sinoIdx=temp*conf->prjWidth+tIdx;
            //    imgVal = tex2D(tex_img,x+N/2,-y+N/2);
            //    sino[sinoIdx]+=weight*imgVal;
            //}
            //temp+=conf->prjFull/4;
            //if(temp<conf->np){
            //    sinoIdx=temp*conf->prjWidth+tIdx;
            //    imgVal = tex2D(tex_img,-y+N/2,-x+N/2);
            //    sino[sinoIdx]+=weight*imgVal;
            //}
            //temp+=conf->prjFull/4;
            //if(temp<conf->np){
            //    sinoIdx=temp*conf->prjWidth+tIdx;
            //    imgVal = tex2D(tex_img,-x+N/2,y+N/2);
            //    sino[sinoIdx]+=weight*imgVal;
            //}

            //if(thetaIdx!=0 && thetaIdx!=conf->prjFull/8){
            //    temp=conf->prjFull/2-thetaIdx;
            //    if(temp<conf->np){
            //        sinoIdx=temp*conf->prjWidth-tIdx+2*pC;
            //        imgVal = tex2D(tex_img,-y+N/2,x+N/2);
            //        sino[sinoIdx]+=weight*imgVal;
            //    }
            //    temp=3*conf->prjFull/4-thetaIdx;
            //    if(temp<conf->np){
            //        sinoIdx=temp*conf->prjWidth-tIdx+2*pC;
            //        imgVal = tex2D(tex_img,x+N/2,y+N/2);
            //        sino[sinoIdx]+=weight*imgVal;
            //    }
            //    temp=conf->prjFull-thetaIdx;
            //    if(temp<conf->np){
            //        sinoIdx=temp*conf->prjWidth-tIdx+2*pC;
            //        imgVal = tex2D(tex_img,y+N/2,-x+N/2);
            //        sino[sinoIdx]+=weight*imgVal;
            //    }
            //    temp=conf->prjFull/4-thetaIdx;
            //    if(temp<conf->np){
            //        sinoIdx=temp*conf->prjWidth-tIdx+2*pC;
            //        imgVal = tex2D(tex_img,-x+N/2,-y+N/2);
            //        sino[sinoIdx]+=weight*imgVal;
            //    }
            //}
        }
        xl -= slopeXYl; xc -= slopeXYc; xr -= slopeXYr;
        dxl=round(MIN(xl,xr)), dxr=round(MAX(xl,xr));
        dxl=MAX(dxl,-(N-1)/2); dxr=MIN(dxr,(N-1)/2);
        y--;
    }
}

int main( void ) {
    int nthread = 256/2;
    int N = 1024/4;
    prjConf config;
    config.nc=N;
    config.nr=N;
    config.prjWidth=N;
    config.np=360;
    config.prjFull=360;
    config.dSize = 1;
    config.effectiveRate=0.9;
    config.d=80000;
    config.fwd=1;

    if(config.prjWidth%2==1){
    }

    config.imgSize=config.nc*config.nr;
    config.sinoSize=config.prjWidth*config.np;

    CPUBitmap image( config.nr, config.nc );
    unsigned char *ptr = image.get_ptr();

    ft* img;
    img = (ft*) malloc(config.imgSize*sizeof(ft));
    int offset;
    int YC = config.nr/2, XC = config.nc/2;
    unsigned char value;
    for(int i=0; i < config.nr; i++)
        for(int j=0; j < config.nc; j++){
            offset = i*config.nc+j;
            if( (i-YC-30)*(i-YC-30)+(j-XC)*(j-XC) < config.nr*config.nc/16){
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

    cudaEvent_t     start, stop;
    ft           elapsedTime;

    // start the timers
    HANDLE_ERROR( cudaEventCreate( &start ) );
    HANDLE_ERROR( cudaEventCreate( &stop ) );

    ft *sino = (ft *) malloc(config.sinoSize*sizeof(ft));
    ft *dev_img;
    ft *dev_sino;
    prjConf *dev_conf;

    cudaChannelFormatDesc desc = cudaCreateChannelDesc<ft>();

    HANDLE_ERROR(cudaMalloc((void**)&dev_img,config.imgSize*sizeof(ft)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_sino,config.sinoSize*sizeof(ft)));
    HANDLE_ERROR(cudaMalloc((void**)&dev_conf, sizeof(prjConf) ) );

    HANDLE_ERROR(cudaMemset(dev_sino,0,config.sinoSize*sizeof(ft)));

    HANDLE_ERROR(cudaBindTexture2D( NULL, tex_img,
                                   dev_img,
                                   desc, config.nr, config.nc,
                                   sizeof(ft) * config.imgSize ) );
    //HANDLE_ERROR( cudaBindTexture2D( NULL, tex_sino,
    //                               dev_sino,
    //                               desc, config.np, config.prjWidth,
    //                               sizeof(ft) * config.sinoSize ) );

    HANDLE_ERROR(cudaMemcpy(dev_img, img, config.imgSize*sizeof(ft), 
                cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemcpy( dev_conf, &config, sizeof(prjConf),
                cudaMemcpyHostToDevice) );

    int sec = MIN(config.np, config.prjFull/8+1);
    dim3 grid(sec,config.prjWidth/nthread);
    cudaProfilerStart();
    HANDLE_ERROR( cudaEventRecord( start, 0 ) );

    //rayDrive<<<sec,config.prjWidth>>>(dev_img, dev_sino, dev_conf);

    printf("grid=(%d,%d), nthread=%d \n",grid.x,grid.y, nthread);
    rayDrive<<<grid,nthread>>>(dev_img, dev_sino, dev_conf);
    HANDLE_ERROR( cudaEventRecord( stop, 0 ) );

    HANDLE_ERROR( cudaEventSynchronize( stop ) );
    cudaProfilerStop();

    HANDLE_ERROR( cudaMemcpy( sino, dev_sino,
                              config.sinoSize*sizeof(ft),
                              cudaMemcpyDeviceToHost ) );

    HANDLE_ERROR( cudaEventElapsedTime( &elapsedTime,
                                        start, stop ) );
    printf( "Time taken:  %3.1f ms\n", elapsedTime );

    CPUBitmap sinogram( config.prjWidth, config.np );
    unsigned char *sinoPtr = sinogram.get_ptr();

    ft temp=0;
    for(int i=0; i < config.np; i++)
        for(int j=0; j < config.prjWidth; j++){
            offset = i*config.prjWidth+j;
            if( sino[offset]>temp){
                //printf("%f",temp);
                temp=sino[offset];
                //printf(" -> %f\n",temp);
            }
        }
    for(int i=0; i < config.np; i++)
        for(int j=0; j < config.prjWidth; j++){
            offset = i*config.prjWidth+j;
            value = (unsigned char)(255 * sino[offset]/temp);
            sinoPtr[(offset<<2)+0] = value;
            sinoPtr[(offset<<2)+1] = value;
            sinoPtr[(offset<<2)+2] = value;
            sinoPtr[(offset<<2)+3] = 0xff;
        }
    cudaUnbindTexture( tex_img );
    cudaUnbindTexture( tex_sino );
    HANDLE_ERROR( cudaFree( dev_img ) );
    HANDLE_ERROR( cudaFree( dev_sino ) );
    HANDLE_ERROR( cudaFree( dev_conf ) );

    HANDLE_ERROR( cudaEventDestroy( start ) );
    HANDLE_ERROR( cudaEventDestroy( stop ) );

    cudaDeviceReset();
    sinogram.display_and_exit();
    return 0;
}

