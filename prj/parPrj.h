
#ifndef parPrj_h
#define parPrj_h
#define DEBUG   0
#define MULTHREAD 1

#include <stddef.h>
#define EPS (1e-8)      // need to justify

typedef struct {
    double bw;  /* X-ray beam width */
    int nc, nr; /* number of rows and cols of input image */
    int np;     /* number of projections at theta */
    int prjWidth;       /* # of elements per projection */
    double *theta;      /* specify the angles from which to project */
    int fwd;    /* To indicate forward projection by 1 and backward by 0 */
} parConf;

int parPrjConf(double *t_img, double *t_prj, double *t_maskIdx, parConf *t_opt,\
        size_t t_totPixels);
int parPrjThread();
int parPrjRun();
void *parPrj(int istart, int iend, int itheta);

#endif
