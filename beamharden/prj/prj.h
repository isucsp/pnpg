
#ifndef __PRJ_H__
#define __PRJ_H__

#include <stddef.h>
#include <float.h>

#define EPS (1e-8)      // need to justify
#define PI (3.14159265359f)
#define SQRT2 (1.414213562373095f)

#define EXE_TIME 1
#define EXE_PROF 0
#define DEBUG   0
#define SHOWIMG  0

#define DIM   1024

// for forward projection

#define THRD_SZ 64
#define LYR_BLK 1

// for backprojection
#define TILE_SZ 8
#define ANG_BLK 1

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define RENEW_MEM (1<<1)
#define FWD_BIT (1)

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

#if GPU
typedef float ft;
#endif

#if CPU
typedef double ft;
#endif

struct prjConf {
    int n; /* number of rows and cols of input image */
    int prjWidth;       /* # of elements per projection */
    int np;     /* number of projections in need */
    int prjFull;       // # of projections in 360 degree
    float dSize;       // size of a single detector
    float effectiveRate;       // effective area rate of the detector

    float d;   // distance from rotation center to X-ray source in pixel
    // set d be FLT_MAX for parallel projection

    char fwd;    /* To indicate forward projection by 1 and backward by 0 */

    unsigned int imgSize;
    unsigned int sinoSize;
};

void setup(int n, int prjWidth, int np, int prjFull, ft dSize, ft effectiveRate, ft d);
void showSetup(void);
#if GPU
int gpuPrj(ft* img, ft* sino, char cmd);
#endif
#if CPU
#endif
int cpuPrj(ft* img, ft* sino, char cmd);

#endif

