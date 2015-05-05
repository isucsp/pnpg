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

/*  call this program by "program filenamebase%03d.sin out.nc" */

int main(int argc, char* argv[])
{
    char filename[100];
    FILE* f;

    if(argc<3){
        printf("error:   ");
        printf("Call this program by \"%s filenamebase%%03d.sin out.nc\"\n",argv[0]);
        return -1;
    }

    char line[MAXLINE];
    //char Shortline[100];
    int i,j,max=0,temp;
    int ncid, x_dimid, y_dimid, z_dimid, varid;
    int dimids[3],retval;
    size_t start[3],count[3];
    int x, y;
    int nx,ny; 
    int datatype=0;
    void* data;
    int* intdata;
    short* shortdata;

    sprintf(filename,argv[1],0);
    printf("%s \n",filename);
    if((f=fopen(filename,"r"))){
        fgets(line,MAXLINE,f);
        fscanf(f,"%d",&nx);
        fscanf(f,"%d",&ny);
        fgets(line,MAXLINE,f);
        printf("Volume Size: %dx%d\n",nx,ny);
        for(j=0; j<3; j++) fgets(line,MAXLINE,f);
        //printf("last line : %s\n",line);
        for(j=0; j<10; j++) for(y=0; y<ny; y++) for(x=0; x<nx; x++){
            fscanf(f,"%d",&temp);
            if(temp>max) max=temp;
        }
    }
    if(max>=SHRT_MAX) datatype=1;

    if ((retval = nc_create(argv[2], NC_NETCDF4|NC_CLOBBER, &ncid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "y", nx, &x_dimid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "x", ny, &y_dimid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "z", NC_UNLIMITED, &z_dimid))) ERR(retval);

    dimids[0] = z_dimid;
    dimids[1] = y_dimid;
    dimids[2] = x_dimid;

    printf("Start to convert sin to NetCDF format ...\n");
    switch(datatype){
        case 0: // short
            printf("Use data type: short ...\n");
            data = calloc(nx*ny, sizeof(short));
            retval = nc_def_var(ncid, "density", NC_SHORT, 3, dimids, &varid);
            break;
        case 1: // int
            printf("Use data type: int ...\n");
            data = calloc(nx*ny, sizeof(int));
            retval = nc_def_var(ncid, "density", NC_INT, 3, dimids, &varid);
            break;
    }
    if (retval) ERR(retval);
    if ((retval = nc_enddef(ncid))) ERR(retval);

    count[0] = 1;
    count[1] = ny;
    count[2] = nx;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    for(i=0; i<100; i++){
        sprintf(filename,argv[1],i);
        printf("%s \n",filename);
        if((f=fopen(filename,"r"))){
            for(j=0; j<5; j++) fgets(line,MAXLINE,f);
            for(j=0; j<10; j++){
                shortdata = (short*)data;
                intdata = (int*)data;

                for(y=0; y<ny; y++)
                    for(x=0; x<nx; x++){
                        fscanf(f,"%d",&temp);
                        switch(datatype){
                            case 0:
                                *shortdata=temp;
                                shortdata++;
                                break;
                            case 1:
                                *intdata=temp;
                                intdata++;
                                break;
                        }
                        if(temp>SHRT_MAX && datatype==0)
                            printf("error: data %d is larger than max %d\n",temp,max);
                        if(temp<0) 
                            printf("error: data %d is negative\n",temp);
                    }
                switch(datatype){
                    case 0:
                        retval = nc_put_vara_short(ncid, varid, start, count, (short*)data);
                        break;
                    case 1:
                        retval = nc_put_vara_int(ncid, varid, start, count, (int*)data);
                        break;
                }
                if(retval) ERR(retval);
                start[0]++;
            }
            fclose(f);
        }else
            break;
    }

    if ((retval = nc_close(ncid))) ERR(retval);
    printf("The maximum data is %d\n",max);

    return 0;

}


