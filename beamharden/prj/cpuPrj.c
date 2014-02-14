/*
 *   Description: use CPU to calculate parallel beam, fan beam projection 
 *   and corresponding back projection.
 *
 *   Reference:
 *   Author: Renliang Gu (renliang@iastate.edu)
 *   $Revision: 0.1 $ $Date: Sun 10 Nov 2013 01:23:32 AM CST
 *
 *   v_0.1:     first draft
 */

static int nthread=64;
static int fSize, bSize;

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
        if(-hbeamW<=rfoot)
            if(hbeamW<=rfoot)
                weight=(height)*(rfoot)*beamWidth/(bias2-bias1);
            else{
                temp=rfoot+hbeamW;
                weight=0.5*height*temp*temp/(bias2-bias1);
            }
    }
    return weight;
}

#if GPU
__global__
#endif
void pixelDriveFan(ft* img, ft* sino, int, int){
    //printf("entering pixelDrive\n");
    // detector array is of odd size with the center at the middle
    ft theta, thetal, thetar;       // the current angle, range in [0 45]
    int thetaIdx, thetalIdx, thetarIdx; // index counts from 0
    prjConf* conf = pConf;

    int pC = conf->prjWidth/2;
    ft d=conf->d; // the distance from rotating center to source in pixel
    int N=conf->n;

    int x=0, y=threadIdx;
    while(y>x){ x++; y-=x; }
    if(x==N/2 && N%2==0){
        for(int i=0; i<N; i++){
            img[i]=0; img[i*N]=0;
        }
        return;
    }

    ft cosT, sinT;  // cosine and sine of theta
    ft cosTl, sinTl;  // cosine and sine of theta
    ft cosTr, sinTr;  // cosine and sine of theta
    ft tileSzSinTl, tileSzCosTl, tileSzSinTr, tileSzCosTr;

    // for each point (x,y) on the ray is
    // x= t*cosT + (t*cosT+d*sinT)/(t*sinT-d*cosT)*(y-t*sinT);
    // or x = -t*d/(t*sinT-d*cosT) + y*(t*cosT+d*sinT)/(t*sinT-d*cosT);
    ft qe, oc;
    ft bq;
    ft cosB, sinB; //, tanB=t/d;
    ft cosR, sinR;
    ft beamWidth = conf->dSize*conf->effectiveRate;
    ft bw;

    ft dist, weight;
    int temp, imgIdx;
    ft imgt[8];
    ft t, tl, tr, tll, trr;
    int dt, dtl, dtr, dtll, dtrr;
    for(int i=0; i<8; i++) imgt[i]=0;
    __shared__ volatile ft shared[8][ANG_BLK][4*TILE_SZ];
    for(thetaIdx=threadIdx.y; thetaIdx<conf->prjFull;thetaIdx+=blockDim.y){

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
        tl = min(tl, t); tr=max(tr,t);

        // bottom right
        //qe = (x+0.5)*sinT - (y+0.5)*cosT + d;
        //oc = (x+0.5)*cosT + (y+0.5)*sinT;
        qe = qe-cosT; oc=oc+sinT; t = d*oc/qe;
        tl = min(tl, t); tr=max(tr,t);

        // bottom left
        //qe = (x-0.5)*sinT - (y+0.5)*cosT + d;
        //oc = (x-0.5)*cosT + (y+0.5)*sinT;
        qe = qe-sinT; oc = oc-cosT; t = d*oc/qe;
        tl = min(tl, t); tr=max(tr,t);

        //qe = d+x*sinT-y*cosT; // positive direction is qo
        qe = qe + sinT/2 +cosT/2;
        dtl = max((int)round((tl+EPS)/conf->dSize),-(conf->prjWidth-1)/2);
        dtr = min((int)round((tr-EPS)/conf->dSize), (conf->prjWidth-1)/2);

        for(dt=dtl; dt<=dtr; dt++){
            t = dt*conf->dSize;
            bq=sqrt(d*d+t*t);
            cosB=d/bq; sinB=t/bq; //, tanB=t/d;
            cosR=cosB*cosT-sinB*sinT;
            sinR=sinB*cosT+cosB*sinT;
            dist=x*cosR+y*sinR-d*t/bq;
            bw = qe*beamWidth*cosB/d;
            weight=getWeight(dist,bw,cosR,sinR);

            imgt[0] += weight*sino[thetaIdx*conf->prjWidth+dt+pC];

            temp=(thetaIdx+conf->prjFull/4)%conf->prjFull;
            imgt[1] += weight*sino[temp*conf->prjWidth+dt+pC];

            temp=(thetaIdx+conf->prjFull/2)%conf->prjFull;
            imgt[2] += weight*sino[temp*conf->prjWidth+dt+pC];

            temp=(thetaIdx+3*conf->prjFull/4)%conf->prjFull;
            imgt[3] += weight*sino[temp*conf->prjWidth+dt+pC];

            temp=(3*conf->prjFull/2-thetaIdx)%conf->prjFull;
            imgt[4] += weight*sino[temp*conf->prjWidth-dt+pC];

            temp=(7*conf->prjFull/4-thetaIdx)%conf->prjFull;
            imgt[5] += weight*sino[temp*conf->prjWidth-dt+pC];

            temp=(conf->prjFull-thetaIdx)%conf->prjFull;
            imgt[6] += weight*sino[temp*conf->prjWidth-dt+pC];

            temp=(5*conf->prjFull/4-thetaIdx)%conf->prjFull;
            imgt[7] += weight*sino[temp*conf->prjWidth-dt+pC];

        }
    }
    if(x<=(N-1)/2 && y<=(N-1)/2){
        imgIdx = ( y+N/2)*N+x+N/2; img[imgIdx] = imgt[0];
        imgIdx = ( x+N/2)*N-y+N/2; img[imgIdx] = imgt[1];
        imgIdx = (-y+N/2)*N-x+N/2; img[imgIdx] = imgt[2];
        imgIdx = (-x+N/2)*N+y+N/2; img[imgIdx] = imgt[3];
        imgIdx = (-y+N/2)*N+x+N/2; img[imgIdx] = imgt[4];
        imgIdx = ( x+N/2)*N+y+N/2; img[imgIdx] = imgt[5];
        imgIdx = ( y+N/2)*N-x+N/2; img[imgIdx] = imgt[6];
        imgIdx = (-x+N/2)*N-y+N/2; img[imgIdx] = imgt[7];
    }
}

__global__ void rayDriveFan(ft* img, ft* sino, int tIdx, int thetaIdx){
    //printf("entering rayDrive\n");
    // detector array is of odd size with the center at the middle
    const prjConf * conf = pConf;
    int sinoIdx;
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
    ft tl = t-hbeamW+EPS,
       tr = t+hbeamW-EPS;
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
        if(QxBxl>0) xl = (dtl -QxBxl *(y+0.5f))/QyByl;
        else xl = (dtl -QxBxl *(y-0.5f))/QyByl;

        xr = (dtr - QxBxr *(y-0.5f))/QyByr;

        dxl =max((int)round(xl ),-(N-1)/2);
        dxr =min((int)round(xr ), (N-1)/2);

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

    sinoIdx=thetaIdx*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[0];

    sinoIdx=(thetaIdx+conf->prjFull/4)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[1];

    sinoIdx=(thetaIdx+conf->prjFull/2)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[2];

    sinoIdx=(thetaIdx+3*conf->prjFull/4)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[3];

    temp = (2*pC-tIdx);
    tIdx = temp%conf->prjWidth;
    temp = 1-temp/conf->prjWidth;

    sinoIdx=(conf->prjFull/2-thetaIdx)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[4]*temp;

    sinoIdx=(3*conf->prjFull/4-thetaIdx)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[5]*temp;

    sinoIdx=(conf->prjFull-thetaIdx)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[6]*temp;

    sinoIdx=(conf->prjFull/4-thetaIdx)*conf->prjWidth+tIdx;
    sino[sinoIdx]=sinot[7]*temp;
}

void *parPrjFor(void *arg){
    int tIdx, thetaIdx;
    for(int i=(size_t)arg; i<pConf->fSize; i+=nthread){
        thetaIdx = i/pConf->prjWidth;
        tIdx = i%pConf->prjWidth;
        rayDrive(pImg,pSino,tIdx,thetaIdx);
    }
    return 0;
}

void *parPrjBack(void *arg){
    int 
    for(int i=(size_t)arg; i<pConf->fSize; i+=nthread){
        pixelDrive(pImg,pSino,x,y);
    }
    return 0;
}

extern "C"
int cpuPrj(ft* img, ft* sino, char cmd){
#if EXE_PROF
    // add some instruction for cpu execution profile
#endif

#if EXE_TIME
    // start the timers
#endif

    if(cmd & FWD_BIT)
    else{
        if(pConf->np<pConf->prjFull)
            HANDLE_ERROR(cudaMemset(dev_sino,0,
                        (pConf->prjFull-pConf->np)*pConf->prjWidth*sizeof(ft)));
    }

#if DEBUG
    if(cmd & FWD_BIT) printf("forward project ...\n");
    else printf("backward project ...\n");
    printf("grid=(%d,%d), thread=(%d,%d)\n",
            grid.x,grid.y,thread.x,thread.y);
#endif

    int res;
    pthread_t *a_thread;
    void *thread_result;
    size_t i=0; 
    threadConfig 

    a_thread=calloc(nthread,sizeof(pthread_t));
    angleSz=(opt->np+nthread-1)/(nthread);
    if(cmd & FWD_BIT){
        /*printf("angleSz=%d given nthread=%d\n",angleSz,nthread);*/
        for(i=0; i<nthread; i++){
            res = pthread_create(&a_thread[i], NULL, parPrjFor, (void *)i);
            if (res != 0) {
                perror("Thread creation failed");
                exit(EXIT_FAILURE);
            }
        }
        /*printf("Waiting for thread to finish...\n");*/
    }else{
        blockSz=(totPixels+nthread-1)/(nthread);
        /*printf("blockSz=%d given nthread=%d\n",blockSz,nthread);*/
        for(i=0; i<nthread; i++){
            res = pthread_create(&a_thread[i], NULL, parPrjBack, (void *)i);
            if (res != 0) {
                perror("Thread creation failed");
                exit(EXIT_FAILURE);
            }
        }
    }
    /*printf("Waiting for thread to finish...\n");*/
    for(i=0; i<nthread; i++){
        res = pthread_join(a_thread[i], &thread_result);
        if (res != 0) {
            perror("Thread join failed");
            exit(EXIT_FAILURE);
        }
        /*printf("Thread joined, it returned %s\n", (char *)thread_result);*/
    }
    free(a_thread);
    if(opt->fwd)
        for(i=0; i<opt->np*opt->prjWidth; i++)
            prj[i]/=beamWidth;
    else
        for(i=0; i<totPixels; i++)
            img[i]/=beamWidth;
    return 0;

#if EXE_PROF
    // add cpu execution profile instructions
#endif

#if EXE_TIME
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
    int offset;
    int YC = config.n/2, XC = config.n/2;
    for(int i=0; i < config.n; i++){
        for(int j=0; j < config.n; j++){
            offset = i*config.n+j;
            if(((i-YC-10)*(i-YC-10)+(j-XC-0)*(j-XC-0)<=160*160)
                    && ((i-YC-60)*(i-YC-60)+(j-XC-60)*(j-XC-60)>=44*44)
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
    FILE* f = fopen("img.data","w");
    fwrite(img,sizeof(ft), pConf->imgSize,f);
    fclose(f);
    //image.display_and_exit();

    ft *sino = (ft *) malloc(config.sinoSize*sizeof(ft));

    gpuPrj(img, sino, RENEW_MEM | FWD_BIT);

    f = fopen("sinogram.data","wb");
    fwrite(sino, sizeof(ft), config.sinoSize, f);
    fclose(f);
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
#if SHOWIMG
            offset = (config.np-1-i)*config.prjWidth+j;
            sinoPtr[(offset<<2)+0] = value;
            sinoPtr[(offset<<2)+1] = value;
            sinoPtr[(offset<<2)+2] = value;
            sinoPtr[(offset<<2)+3] = 0xff;
#endif
        }
#endif
    free(sino); free(img);
#if SHOWIMG
    sinogram.display_and_exit();
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
#endif

    ft* img = (ft*) malloc(config.imgSize*sizeof(ft));
    int offset;
    ft *sino = (ft *) malloc(config.sinoSize*sizeof(ft));

    FILE* f = fopen("sinogram.data","rb");
    int tempI;
    tempI=rand()%config.np;

    fread(sino,sizeof(ft),config.sinoSize,f);

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
    }
    fclose(f);
    //sinogram.display_and_exit();

    gpuPrj(img,sino, RENEW_MEM);

    f = fopen("reImg.data","wb");
    fwrite(img,sizeof(ft),config.imgSize,f);
    fclose(f);

#if SHOWIMG
    ft temp=0;
    for(int i=0; i < config.n; i++)
        for(int j=0; j < config.n; j++){
            offset = i*config.n+j;
            if( img[offset]>temp) temp=img[offset];
        }
    for(int i=0; i < config.n; i++){
        for(int j=0; j < config.n; j++){
            offset = i*config.n+j;
            value = (unsigned char)(255 * img[offset]/temp);
#if SHOWIMG
            offset = (config.n-1-i)*config.n+j;
            ptr[(offset<<2)+0] = value;
            ptr[(offset<<2)+1] = value;
            ptr[(offset<<2)+2] = value;
            ptr[(offset<<2)+3] = 0xff;
#endif
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
    int N = 1024;
    if(argc>1){
        N = atoi(argv[1]);
        printf("N=%d\n",N);
    }
    setup(N,N,360,360,1,1,3000);

    forwardTest();
    backwardTest();
    cleanUp();
    return 0;
}

