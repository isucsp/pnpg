#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <netcdf.h>

/* this is an example program writes data in NetCDF */

#define FILE_NAME "volume.nc"
#define ERRCODE 2
#define MAXLINE (1400*20)
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define DATATYPE float
#define FLOAT 1
#define SHORT 0
#define INT 0

/*  call this program by "program filenamebase%03d.sin out.nc" */

int main(int argc, char* argv[])
{
    char* filename=argv[1];
    FILE* f;

    if(argc<3){
        printf("error:   ");
        printf("Call this program by \"%s volfile.vol out.nc\"\n",argv[0]);
        return -1;
    }

    char line[MAXLINE];
    long int nx,ny,nz;
    //char Shortline[100];
    int i,j,max=0,temp;
    int ncid, x_dimid, y_dimid, z_dimid, varid;
    int dimids[3],retval;
    size_t start[3],count[3];
    long int filesize;
    int x, y,z;
    DATATYPE* data,p;

    printf("input filename: %s \n",filename);
    if((f=fopen(filename,"rb"))){
        fseek(f,0,SEEK_END);
        filesize=ftell(f);
        printf("filesize=%ld\n",filesize);
        rewind(f);
        temp=0;
        while(temp!=3){
            i=0;
            while((line[i]=fgetc(f))!=0 && line[i]!=13 && line[i]!=10) i++;
            //printf("%s\n",line);
            temp=sscanf(line,"Volume Size: %ldx%ldx%ld",&ny,&nx,&nz);
        }
        printf("Volume Size: %ldx%ldx%ld\n",nx,ny,nz);
        fseek(f,filesize%(nx*ny*nz),SEEK_SET);

        data = (DATATYPE*)calloc(nx*ny, sizeof(DATATYPE));
        if(data==NULL){ perror("failed calloc call\n"); return -1; }

        if ((retval = nc_create(argv[2], NC_NETCDF4|NC_CLOBBER, &ncid))) ERR(retval);
        if ((retval = nc_def_dim(ncid, "x", nx, &x_dimid))) ERR(retval);
        if ((retval = nc_def_dim(ncid, "y", ny, &y_dimid))) ERR(retval);
        if ((retval = nc_def_dim(ncid, "z", nz, &z_dimid))) ERR(retval);

        dimids[0] = z_dimid;
        dimids[1] = x_dimid;
        dimids[2] = y_dimid;

        printf("Start to convert sin to NetCDF format ...\n");
#if FLOAT
        if ((retval = nc_def_var(ncid, "density", NC_FLOAT, 3, dimids, &varid))) ERR(retval);
#endif
        if ((retval = nc_enddef(ncid))) ERR(retval);

        count[0] = 1;
        count[1] = nx;
        count[2] = ny;
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;

        memset(line,0,MAXLINE);
        for(z=0; z<nz; z++){
            temp=fread(data,sizeof(DATATYPE),nx*ny,f);
            if(temp!=(nx*ny)){
                printf("readed %d bytes intead of %ld\n",temp,nx*ny*sizeof(DATATYPE));
                return -1;
            }
            //printf("data: %f %f %f %f\n", data[0], data[1], data[2], data[3]);
            start[0]=z;
#if FLOAT
            if ((retval = nc_put_vara_float(ncid, varid, start, count, data))) ERR(retval);
#endif
            for(i=0; i<strlen(line); i++) printf("\b");
            sprintf(line,"z=%d/%ld",z,nz);
            printf("%s",line);
            fflush(stdout);
        }
        fclose(f);
        if ((retval = nc_close(ncid))) ERR(retval);
        printf("\n");
    }

    return 0;
}


