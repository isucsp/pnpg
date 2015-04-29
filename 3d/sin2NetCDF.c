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
#define nx 1400
#define ny 360
#define nz (10*45)
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
    short data[nx][ny];

    if ((retval = nc_create(argv[2], NC_NETCDF4|NC_CLOBBER, &ncid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "x", nx, &x_dimid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "y", ny, &y_dimid))) ERR(retval);
    if ((retval = nc_def_dim(ncid, "z", NC_UNLIMITED, &z_dimid))) ERR(retval);

    dimids[0] = z_dimid;
    dimids[1] = x_dimid;
    dimids[2] = y_dimid;

    printf("Start to convert sin to NetCDF format ...\n");
    if ((retval = nc_def_var(ncid, "density", NC_SHORT, 3, dimids, &varid))) ERR(retval);
    if ((retval = nc_enddef(ncid))) ERR(retval);

    count[0] = 1;
    count[1] = nx;
    count[2] = ny;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    for(i=0; i<100; i++){
        sprintf(filename,argv[1],i);
        printf("%s \n",filename);
        if((f=fopen(filename,"r"))){
            for(j=0; j<5; j++) fgets(line,MAXLINE,f);
            for(j=0; j<10; j++){
                for(y=0; y<ny; y++)
                    for(x=0; x<nx; x++){
                        fscanf(f,"%d",&temp);
                        data[x][y]=temp;
                        if(data[x][y]>max) max=data[x][y];
                        if(temp>=SHRT_MAX)
                            printf("error: data %d is larger than SHRT_MAX %d\n",temp,SHRT_MAX);
                        if(temp<0) 
                            printf("error: data %d is negative\n",temp);
                    }
                if ((retval = nc_put_vara_short(ncid, varid, start, count, &data[0][0]))) ERR(retval);
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


