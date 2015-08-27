
#ifndef __UTILS_H__
#define __UTILS_H__

#include <stddef.h>
#include <float.h>
#include "prj.h"

#if _WIN32
ft fabs(ft x){
    return fabsf(x);
}
#endif

#if __NVCC__
__device__
#endif
ft getWeight(ft dist, ft beamWidth, ft cosR, ft sinR){
    // range of gamma is in [-alpha,pi/4+alpha], where alpha is small
    //ft temp=abs(fmod(theta,90)-45)*PI/180; /* test value of temp */
    cosR = cosR<0 ? -cosR : cosR;
    sinR = sinR<0 ? -sinR : sinR;
    ft temp;
    if(cosR<sinR){
        temp=cosR; cosR=sinR; sinR=temp;
    }

    ft height=1/cosR;
    ft bias2=(cosR+sinR)/2;    //bias2 should be larger.
    ft bias1=fabs(cosR-sinR)/2;

    ft areaT=0.5*height*(bias2-bias1);
    ft areaR=height*bias1*2;
    ft area=height*(bias1+bias2);
    ft weight=0;

    // if pixel is to the right of the ray, dist > 0, otherwise dist<=0
    // due to the symmetry, take the abs of dist
    dist = dist < 0 ? -dist : dist;
    ft lfoot=dist-bias2, rfoot=dist+bias2;
    ft lshoulder = dist-bias1, rshoulder=dist+bias1;
    //printf("lf ls rs rf: (%f, %f, %f, %f)\n",lfoot, lshoulder, rshoulder, rfoot);
    ft hbeamW = beamWidth/2;
    if(-hbeamW<=lfoot){
        //printf("test1\n");
        if (hbeamW<=lfoot) weight=0;
        else if(hbeamW<=lshoulder){
            temp = (hbeamW-(lfoot));
            weight = 0.5*temp*temp*height/(bias2-bias1);
        }else if(hbeamW<=rshoulder){
            weight=areaT+height*(hbeamW-lshoulder);
            //printf("areaT=%f, height=%f, weight=%f\n",areaT, height, weight);
        }else if(hbeamW<=rfoot){
            temp = rfoot-hbeamW;
            weight=areaT+areaR+areaT
                   -0.5*temp*temp*height/(bias2-bias1);
        }else weight=area;
    }else if(-hbeamW<=lshoulder){
        if(hbeamW<=lshoulder)
            weight=height*(-lfoot)*beamWidth/(bias2-bias1);
        else if(hbeamW<=rshoulder){
            temp = -hbeamW-lfoot;
            weight=areaT +height*(hbeamW-(lshoulder))
                - 0.5*temp*height*temp/(bias2-bias1);
        }else if(hbeamW<=rfoot){
            temp = -hbeamW-lfoot;
            weight=areaT + areaR
                - 0.5*temp*height*temp/(bias2-bias1);
            temp = rfoot - hbeamW;
            weight = weight + areaT
                - 0.5*temp*height*temp/(bias2-bias1);
        }else{
            temp = -hbeamW-lfoot;
            weight=areaT - 0.5*temp*height*temp/(bias2-bias1)
                +areaR+areaT;
        }
    }else if(-hbeamW<=rshoulder){
        if(hbeamW<=rshoulder)
            weight=height*beamWidth;
        else if(hbeamW<=rfoot){
            temp=rfoot-hbeamW;
            weight=height*(rshoulder+hbeamW)+ areaT
                - 0.5*temp*height*temp/(bias2-bias1);
        }else
            weight=height*(rshoulder+hbeamW)+areaT;
    }else{
        if(-hbeamW<=rfoot){
            if(hbeamW<=rfoot)
                weight=(height)*(rfoot)*beamWidth/(bias2-bias1);
            else{
                temp=rfoot+hbeamW;
                weight=0.5*height*temp*temp/(bias2-bias1);
            }
        }
    }
    return weight;
}

void rampFilter(ft *signal, int size, ft Ts){
    int N = 2*size;
    kiss_fft_cfg cfgFFT = kiss_fft_alloc(N,0,0,0);
    kiss_fft_cfg cfgIFFT = kiss_fft_alloc(N,1,0,0);
    kiss_fft_cpx ramp[N];
    kiss_fft_cpx hann[N];
    kiss_fft_cpx proj[N];
    int i;
    for(i=0; i<N; i++){
        //hamming=0.54+0.46*cos(2*pi*n'/N);
        //hann=0.5+0.5*cos(2*pi*n'/N);
        //sinc=sin(pi*n'/N)./(pi*n'/N);
        if(i<N/2){
            ramp[i].r = (ft)i/N;
            proj[i].r = signal[i];
        }else{
            ramp[i].r = (ft)(N-i)/N;
            proj[i].r = 0;
        }
        hann[i].r = 0.5+0.5*cos(2*PI*i/N);
        ramp[i].i = 0;
        proj[i].i = 0;
        hann[i].i = 0;
    }
    // use double length of ramp filter to get the best performance.
    for(i=1; i<N; i++){
        if(i%2==0)
            ramp[i].r=0;
        else{
            if(i<N/2)
                ramp[i].r = -1.0/PI/PI/i/i;
            else
                ramp[i].r = -1.0/PI/PI/(N-i)/(N-i);
        }
        //if(i>=N/4 && i<N*3/4) ramp[i].r=0;
    }
    ramp[0].r=0.25;
    kiss_fft( cfgFFT , ramp, ramp );

    kiss_fft( cfgFFT , proj, proj );
    for(i=0; i<N; i++){
        proj[i].r=proj[i].r*ramp[i].r*hann[i].r;
        proj[i].i=proj[i].i*ramp[i].r*hann[i].r;
    }
    // the kiss_fft inverse fft doesn't divide N
    kiss_fft( cfgIFFT, proj, proj );
    for(i=0; i<N/2; i++){
        signal[i] = proj[i].r/N;
    }

    //kiss_fft( cfgIFFT , ramp, ramp );
    //for(int i=0; i<N; i++){
    //    printf("%d, %g, %g, %g, %g, %g\n",i,ramp[i].r/N,ramp[i].i, hann[i].r, hann[i].i, proj[i].r);
    //}

    free(cfgFFT); free(cfgIFFT);
}

#endif

