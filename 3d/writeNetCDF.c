#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <netcdf.h>

/* this program writes data in NetCDF */

#define FILE_NAME "volume.nc"
#define nx 30
#define ny 30
#define nz 30
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int main()
{
    int ncid, x_dimid, y_dimid, z_dimid, varid;
    int dimids[3], retval, i, j, k;
    float data[nx][ny][nz], x, y, z, r, r0 = 0.4;

    for (i = 0; i < nx; i++) {
        x = ((float)(i)+0.5)/(float)(nx);
        for (j = 0; j < ny; j++) {
            y = ((float)(j)+0.5)/(float)(ny);
            for (k = 0; k < nz; k++) {
                z = ((float)(k)+0.5)/(float)(nz);
                r = sqrt(pow(x-0.5,2)+pow(y-0.5,2));
                data[i][j][k] = exp(-sqrt(pow(z-0.5,2)+pow(r-r0,2)));
            }
        }
    }

    if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &ncid)))
        ERR(retval);

    if ((retval = nc_def_dim(ncid, "z", nz, &z_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "y", ny, &y_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "x", nx, &x_dimid)))
        ERR(retval);

    dimids[0] = x_dimid;
    dimids[1] = y_dimid;
    dimids[2] = z_dimid;

    if ((retval = nc_def_var(ncid, "density", NC_FLOAT, 3, dimids, &varid)))
        ERR(retval);

    if ((retval = nc_enddef(ncid)))
        ERR(retval);

    if ((retval = nc_put_var_float(ncid, varid, &data[0][0][0])))
        ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);

    return 0;
}
