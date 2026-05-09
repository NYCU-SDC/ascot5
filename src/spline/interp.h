/**
 * @file interp.h
 * @brief Spline interpolation library
 *
 * Spline interpolation fits cubic splines on data given on a uniform grid. Each
 * axis may have either natural or periodic boundary condition.
 *
 * There exists two representations for the splines: compact and explicit. Both
 * give identical results but the difference is that compact requires fewer
 * coefficients to be stored (and fetched from the memory at each evaluation)
 * than the explicit. On the other hand, explicit requires less computations
 * per evaluations, but compact is usually faster, and conserves memory,
 * especially in 3D. Therefore compact splines are preferred.
 *
 * To initialize splines, first call corresponding init_coeff() function which
 * evaluates coefficients to a pre-allocated array (i.e. to the offload array).
 * Then (after offloading is done) call init_spline() which assigns the
 * coefficients to a spline struct.
 *
 * In order to allocate the array for storing the coefficients, one needs to
 * know how many coefficients are stored per data grid point:
 *
 * - 1D compact  2, explicit 4
 * - 2D compact  4, explicit 16
 * - 3D compact  8, explicit 64
 */
#ifndef INTERP_H
#define INTERP_H
#include <math.h>
#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"

/**
 * @brief Boundary conditions for the spline interpolation.
 */
enum boundaryCondition {
    NATURALBC  = 0, /**< Second derivative is zero at both ends               */
    PERIODICBC = 1  /**< Function has same value and derivatives at both ends */
};

/**
 * @brief Number of coefficients stored for each data point.
 */
enum splinesize {
    NSIZE_COMP1D =  2,
    NSIZE_COMP2D =  4,
    NSIZE_COMP3D =  8,
    NSIZE_EXPL1D =  4,
    NSIZE_EXPL2D = 16,
    NSIZE_EXPL3D = 64
};

/**
 * @brief Cubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp1D_data;

/**
 * @brief Bicubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int n_y;     /**< number of y grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real y_min;  /**< minimum y coordinate in the grid               */
    real y_max;  /**< maximum y coordinate in the grid               */
    real y_grid; /**< interval between two adjacent points in y grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp2D_data;

/**
 * @brief Tricubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int n_y;     /**< number of y grid points                        */
    int n_z;     /**< number of z grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    int bc_z;    /**< boundary condition for z coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real y_min;  /**< minimum y coordinate in the grid               */
    real y_max;  /**< maximum y coordinate in the grid               */
    real y_grid; /**< interval between two adjacent points in y grid */
    real z_min;  /**< minimum z coordinate in the grid               */
    real z_max;  /**< maximum z coordinate in the grid               */
    real z_grid; /**< interval between two adjacent points in z grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp3D_data;

int interp1Dcomp_init_coeff(real* c, real* f,
                            int n_x, int bc_x,
                            real x_min, real x_max);

int interp2Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

int interp3Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max);

int interp1Dexpl_init_coeff(real* c, real* f,
                            int n_x, int bc_x,
                            real x_min, real x_max);

int interp2Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

int interp3Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max);

void interp1Dcomp_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x,
                              real x_min, real x_max);

void interp2Dcomp_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);

void interp3Dcomp_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max);

void interp1Dexpl_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x,
                              real x_min, real x_max);

void interp2Dexpl_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);

void interp3Dexpl_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max);

int interp1Dcomp_setup(interp1D_data* str, real* f, int n_x, int bc_x,
                       real x_min, real x_max);

int interp2Dcomp_setup(interp2D_data* str, real* f, int n_x, int n_y,
                       int bc_x, int bc_y, real x_min, real x_max,
                       real y_min, real y_max);

int interp3Dcomp_setup(interp3D_data* str, real* f,
                       int n_x, int n_y, int n_z, int bc_x, int bc_y, int bc_z,
                       real x_min, real x_max, real y_min, real y_max,
                       real z_min, real z_max);

static inline a5err interp1Dcomp_eval_f(real* f, interp1D_data* str, real x) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x   = (x-str->x_min) / str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx   = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    /* Helper varibles */
    real dx3  =  dx * (dx*dx - 1.0);
    real dxi  = 1.0 - dx;
    real dxi3 = dxi * (dxi*dxi - 1.0);
    real xg2  = str->x_grid*str->x_grid;

    int n  = i_x*2; /* Index jump to cell       */
    int x1 = 2;     /* Index jump one x forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }

    if(!err) {
        *f =
                      dxi *str->c[n+0]+dx *str->c[n+x1+0]
            +(xg2/6)*(dxi3*str->c[n+1]+dx3*str->c[n+x1+1]);
    }

    return err;
}

static inline a5err interp2Dcomp_eval_f(real* f, interp2D_data* str, real x, real y) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x   = (x - str->x_min) / str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx   = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    /* Helper variables */
    real dx3  =  dx * (dx*dx - 1.0);
    real dxi  = 1.0 - dx;
    real dxi3 = dxi * (dxi*dxi - 1.0);
    real xg2  = str->x_grid*str->x_grid;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y   = (y - str->y_min) / str->y_grid - 1*(y==str->y_max);
    /* Normalized y coordinate in current cell */
    real dy   = ( y - (str->y_min + i_y*str->y_grid) ) / str->y_grid;
    /* Helper variables */
    real dy3  =  dy * (dy*dy - 1.0);
    real dyi  = 1.0 - dy;
    real dyi3 = dyi * (dyi*dyi - 1.0);
    real yg2  = str->y_grid*str->y_grid;

    int n  = i_y*str->n_x*4+i_x*4; /* Index jump to cell       */
    int x1 = 4;                    /* Index jump one x forward */
    int y1 = str->n_x*4;           /* Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the domain. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max) ) {
        err = 1;
    }

    if(!err) {
        *f = (
            dxi*(dyi*str->c[n]+dy*str->c[n+y1])
            +dx*(dyi*str->c[n+x1]+dy*str->c[n+y1+x1]))
            +(xg2/6)*(
                dxi3*(dyi*str->c[n+1] + dy*str->c[n+y1+1])
                +dx3*(dyi*str->c[n+x1+1] + dy*str->c[n+y1+x1+1]))
            +(yg2/6)*(
                dxi*(dyi3*str->c[n+2]+dy3*str->c[n+y1+2])
                +dx*(dyi3*str->c[n+x1+2]+dy3*str->c[n+y1+x1+2]))
            +(xg2*yg2/36)*(
                dxi3*(dyi3*str->c[n+3]+dy3*str->c[n+y1+3])
                +dx3*(dyi3*str->c[n+x1+3]+dy3*str->c[n+y1+x1+3]));
    }

    return err;
}

static inline real interp3Dcomp_eval_f_precomputed(
    const real* c,
    int n, int x1, int y1, int z1,
    real dxi, real dx, real dxi3, real dx3,
    real dyi, real dy, real dyi3, real dy3,
    real dzi, real dz, real dzi3, real dz3,
    real xg2, real yg2, real zg2) {
    /* Preload all 64 coefficients to batch gathers and avoid aliasing reloads. */
    real c0000 = c[n+0];
    real c0001 = c[n+1];
    real c0002 = c[n+2];
    real c0003 = c[n+3];
    real c0004 = c[n+4];
    real c0005 = c[n+5];
    real c0006 = c[n+6];
    real c0007 = c[n+7];

    real c0010 = c[n+x1+0];
    real c0011 = c[n+x1+1];
    real c0012 = c[n+x1+2];
    real c0013 = c[n+x1+3];
    real c0014 = c[n+x1+4];
    real c0015 = c[n+x1+5];
    real c0016 = c[n+x1+6];
    real c0017 = c[n+x1+7];

    real c0100 = c[n+y1+0];
    real c0101 = c[n+y1+1];
    real c0102 = c[n+y1+2];
    real c0103 = c[n+y1+3];
    real c0104 = c[n+y1+4];
    real c0105 = c[n+y1+5];
    real c0106 = c[n+y1+6];
    real c0107 = c[n+y1+7];

    real c1000 = c[n+z1+0];
    real c1001 = c[n+z1+1];
    real c1002 = c[n+z1+2];
    real c1003 = c[n+z1+3];
    real c1004 = c[n+z1+4];
    real c1005 = c[n+z1+5];
    real c1006 = c[n+z1+6];
    real c1007 = c[n+z1+7];

    real c0110 = c[n+y1+x1+0];
    real c0111 = c[n+y1+x1+1];
    real c0112 = c[n+y1+x1+2];
    real c0113 = c[n+y1+x1+3];
    real c0114 = c[n+y1+x1+4];
    real c0115 = c[n+y1+x1+5];
    real c0116 = c[n+y1+x1+6];
    real c0117 = c[n+y1+x1+7];

    real c1010 = c[n+z1+x1+0];
    real c1011 = c[n+z1+x1+1];
    real c1012 = c[n+z1+x1+2];
    real c1013 = c[n+z1+x1+3];
    real c1014 = c[n+z1+x1+4];
    real c1015 = c[n+z1+x1+5];
    real c1016 = c[n+z1+x1+6];
    real c1017 = c[n+z1+x1+7];

    real c1100 = c[n+z1+y1+0];
    real c1101 = c[n+z1+y1+1];
    real c1102 = c[n+z1+y1+2];
    real c1103 = c[n+z1+y1+3];
    real c1104 = c[n+z1+y1+4];
    real c1105 = c[n+z1+y1+5];
    real c1106 = c[n+z1+y1+6];
    real c1107 = c[n+z1+y1+7];

    real c1110 = c[n+z1+y1+x1+0];
    real c1111 = c[n+z1+y1+x1+1];
    real c1112 = c[n+z1+y1+x1+2];
    real c1113 = c[n+z1+y1+x1+3];
    real c1114 = c[n+z1+y1+x1+4];
    real c1115 = c[n+z1+y1+x1+5];
    real c1116 = c[n+z1+y1+x1+6];
    real c1117 = c[n+z1+y1+x1+7];

    return (
        dzi*(
            dxi*(dyi*c0000+dy*c0100)
            +dx*(dyi*c0010+dy*c0110))
        +dz*(
            dxi*(dyi*c1000+dy*c1100)
            +dx*(dyi*c1010+dy*c1110)))
        +xg2/6*(
            dzi*(
                dxi3*(dyi*c0001+dy*c0101)
                +dx3*(dyi*c0011+dy*c0111))
            +dz*(
                dxi3*(dyi*c1001+dy*c1101)
                +dx3*(dyi*c1011+dy*c1111)))
        +yg2/6*(
            dzi*(
                dxi*(dyi3*c0002+dy3*c0102)
                +dx*(dyi3*c0012+dy3*c0112))
            +dz*(
                dxi*(dyi3*c1002+dy3*c1102)
                +dx*(dyi3*c1012+dy3*c1112)))
        +zg2/6*(
            dzi3*(
                dxi*(dyi*c0003+dy*c0103)
                +dx*(dyi*c0013+dy*c0113))
            +dz3*(
                dxi*(dyi*c1003+dy*c1103)
                +dx*(dyi*c1013+dy*c1113)))
        +xg2*yg2/36*(
            dzi*(
                dxi3*(dyi3*c0004+dy3*c0104)
                +dx3*(dyi3*c0014+dy3*c0114))
            +dz*(
                dxi3*(dyi3*c1004+dy3*c1104)
                +dx3*(dyi3*c1014+dy3*c1114)))
        +xg2*zg2/36*(
            dzi3*(
                dxi3*(dyi*c0005+dy*c0105)
                +dx3*(dyi*c0015+dy*c0115))
            +dz3*(
                dxi3*(dyi*c1005+dy*c1105)
                +dx3*(dyi*c1015+dy*c1115)))
        +yg2*zg2/36*(
            dzi3*(
                dxi*(dyi3*c0006+dy3*c0106)
                +dx*(dyi3*c0016+dy3*c0116))
            +dz3*(
                dxi*(dyi3*c1006+dy3*c1106)
                +dx*(dyi3*c1016+dy3*c1116)))
        +xg2*yg2*zg2/216*(
            dzi3*(
                dxi3*(dyi3*c0007+dy3*c0107)
                +dx3*(dyi3*c0017+dy3*c0117))
            +dz3*(
                dxi3*(dyi3*c1007+dy3*c1107)
                +dx3*(dyi3*c1017+dy3*c1117)));
}

static inline a5err interp3Dcomp_eval_f(real* f, interp3D_data* str,
                                        real x, real y, real z) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }
    if(str->bc_z == PERIODICBC) {
        z = fmod(z - str->z_min, str->z_max - str->z_min) + str->z_min;
        z = z + (z < str->z_min) * (str->z_max - str->z_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x   = (x - str->x_min) / str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx   = (x - (str->x_min + i_x*str->x_grid)) / str->x_grid;
    /* Helper varibles */
    real dxi  = 1.0 - dx;
    real dx3  = dx*dx*dx - dx;
    real dxi3 = (1.0 - dx) * (1.0 - dx) * (1.0 - dx) - (1.0 - dx);
    real xg2  = str->x_grid*str->x_grid;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y   = (y - str->y_min) / str->y_grid - 1*(y==str->y_max);
    /* Normalized y coordinate in current cell */
    real dy   = (y - (str->y_min + i_y*str->y_grid)) / str->y_grid;
    /* Helper varibles */
    real dyi  = 1.0 - dy;
    real dy3  = dy*dy*dy - dy;
    real dyi3 = (1.0 - dy) * (1.0 - dy) * (1.0 - dy) - (1.0 - dy);
    real yg2  = str->y_grid*str->y_grid;

    /* Index for z variable. The -1 needed at exactly grid end. */
    int i_z   = (z - str->z_min) / str->z_grid - 1*(z==str->z_max);
    /* Normalized z coordinate in current cell */
    real dz   = (z - (str->z_min + i_z*str->z_grid)) / str->z_grid;
    /* Helper varibles */
    real dzi  = 1.0 - dz;
    real dz3  = dz*dz*dz - dz;
    real dzi3 = (1.0 - dz) * (1.0 - dz) * (1.0 - dz) - (1.0-dz);
    real zg2  = str->z_grid*str->z_grid;

    /**< Index jump to cell */
    int n  = i_z*str->n_y*str->n_x*8 + i_y*str->n_x*8 + i_x*8;
    int x1 = 8;                   /* Index jump one x forward */
    int y1 = str->n_x*8;          /* Index jump one y forward */
    int z1 = str->n_y*str->n_x*8; /* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the domain. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max) ) {
        err = 1;
    }
    if( str->bc_z == PERIODICBC && i_z == str->n_z-1 ) {
        z1 = -(str->n_z-1)*z1;
    }
    else if( str->bc_z == NATURALBC && !(z >= str->z_min && z <= str->z_max) ) {
        err = 1;
    }

    if(!err) {
        /* Evaluate spline value */
        *f = (
            dzi*(
                dxi*(dyi*str->c[n+0]+dy*str->c[n+y1+0])
                +dx*(dyi*str->c[n+x1+0]+dy*str->c[n+y1+x1+0]))
            +dz*(
                dxi*(dyi*str->c[n+z1+0]+dy*str->c[n+y1+z1+0])
                +dx*(dyi*str->c[n+x1+z1+0]+dy*str->c[n+y1+z1+x1+0]))
            +xg2/6*(
                dzi*(
                    dxi3*(dyi*str->c[n+1]+dy*str->c[n+y1+1])
                    +dx3*(dyi*str->c[n+x1+1]+dy*str->c[n+y1+x1+1]))
                +dz*(
                    dxi3*(dyi*str->c[n+z1+1]+dy*str->c[n+y1+z1+1])
                    +dx3*(dyi*str->c[n+x1+z1+1]+dy*str->c[n+y1+z1+x1+1])))
            +yg2/6*(
                dzi*(
                    dxi*(dyi3*str->c[n+2]+dy3*str->c[n+y1+2])
                    +dx*(dyi3*str->c[n+x1+2]+dy3*str->c[n+y1+x1+2]))
                +dz*(
                    dxi*(dyi3*str->c[n+z1+2]+dy3*str->c[n+y1+z1+2])
                    +dx*(dyi3*str->c[n+x1+z1+2]+dy3*str->c[n+y1+z1+x1+2])))
            +zg2/6*(
                dzi3*(
                    dxi*(dyi*str->c[n+3]+dy*str->c[n+y1+3])
                    +dx*(dyi*str->c[n+x1+3]+dy*str->c[n+y1+x1+3]))
                +dz3*(
                    dxi*(dyi*str->c[n+z1+3]+dy*str->c[n+y1+z1+3])
                    +dx*(dyi*str->c[n+x1+z1+3]+dy*str->c[n+y1+z1+x1+3])))
            +xg2*yg2/36*(
                dzi*(
                    dxi3*(dyi3*str->c[n+4]+dy3*str->c[n+y1+4])
                    +dx3*(dyi3*str->c[n+x1+4]+dy3*str->c[n+y1+x1+4]))
                +dz*(
                    dxi3*(dyi3*str->c[n+z1+4]+dy3*str->c[n+y1+z1+4])
                    +dx3*(dyi3*str->c[n+x1+z1+4]+dy3*str->c[n+y1+z1+x1+4])))
            +xg2*zg2/36*(
                dzi3*(
                    dxi3*(dyi*str->c[n+5]+dy*str->c[n+y1+5])
                    +dx3*(dyi*str->c[n+x1+5]+dy*str->c[n+y1+x1+5]))
                +dz3*(
                    dxi3*(dyi*str->c[n+z1+5]+dy*str->c[n+y1+z1+5])
                    +dx3*(dyi*str->c[n+x1+z1+5]+dy*str->c[n+y1+z1+x1+5])))
            +yg2*zg2/36*(
                dzi3*(
                    dxi*(dyi3*str->c[n+6]+dy3*str->c[n+y1+6])
                    +dx*(dyi3*str->c[n+x1+6]+dy3*str->c[n+y1+x1+6]))
                +dz3*(
                    dxi*(dyi3*str->c[n+z1+6]+dy3*str->c[n+y1+z1+6])
                    +dx*(dyi3*str->c[n+x1+z1+6]+dy3*str->c[n+y1+z1+x1+6])))
            +xg2*yg2*zg2/216*(
                dzi3*(
                    dxi3*(dyi3*str->c[n+7]+dy3*str->c[n+y1+7])
                    +dx3*(dyi3*str->c[n+x1+7]+dy3*str->c[n+y1+x1+7]))
                +dz3*(
                    dxi3*(dyi3*str->c[n+z1+7]+dy3*str->c[n+y1+z1+7])
                    +dx3*(dyi3*str->c[n+x1+z1+7]+dy3*str->c[n+y1+z1+x1+7]))));

    }

    return err;
}

static inline a5err interp3Dcomp_eval_f3(real f[3],
                           interp3D_data* str0,
                           interp3D_data* str1,
                           interp3D_data* str2,
                           real x, real y, real z) {
    interp3D_data* str = str0;

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }
    if(str->bc_z == PERIODICBC) {
        z = fmod(z - str->z_min, str->z_max - str->z_min) + str->z_min;
        z = z + (z < str->z_min) * (str->z_max - str->z_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x   = (x - str->x_min) / str->x_grid - 1*(x==str->x_max);
    real dx   = (x - (str->x_min + i_x*str->x_grid)) / str->x_grid;
    real dxi  = 1.0 - dx;
    real dx3  = dx*dx*dx - dx;
    real dxi3 = (1.0 - dx) * (1.0 - dx) * (1.0 - dx) - (1.0 - dx);
    real xg2  = str->x_grid*str->x_grid;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y   = (y - str->y_min) / str->y_grid - 1*(y==str->y_max);
    real dy   = (y - (str->y_min + i_y*str->y_grid)) / str->y_grid;
    real dyi  = 1.0 - dy;
    real dy3  = dy*dy*dy - dy;
    real dyi3 = (1.0 - dy) * (1.0 - dy) * (1.0 - dy) - (1.0 - dy);
    real yg2  = str->y_grid*str->y_grid;

    /* Index for z variable. The -1 needed at exactly grid end. */
    int i_z   = (z - str->z_min) / str->z_grid - 1*(z==str->z_max);
    real dz   = (z - (str->z_min + i_z*str->z_grid)) / str->z_grid;
    real dzi  = 1.0 - dz;
    real dz3  = dz*dz*dz - dz;
    real dzi3 = (1.0 - dz) * (1.0 - dz) * (1.0 - dz) - (1.0-dz);
    real zg2  = str->z_grid*str->z_grid;

    int n  = i_z*str->n_y*str->n_x*8 + i_y*str->n_x*8 + i_x*8;
    int x1 = 8;
    int y1 = str->n_x*8;
    int z1 = str->n_y*str->n_x*8;

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the domain. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max) ) {
        err = 1;
    }
    if( str->bc_z == PERIODICBC && i_z == str->n_z-1 ) {
        z1 = -(str->n_z-1)*z1;
    }
    else if( str->bc_z == NATURALBC && !(z >= str->z_min && z <= str->z_max) ) {
        err = 1;
    }

    if(!err) {
        f[0] = interp3Dcomp_eval_f_precomputed(
            str0->c, n, x1, y1, z1,
            dxi, dx, dxi3, dx3,
            dyi, dy, dyi3, dy3,
            dzi, dz, dzi3, dz3,
            xg2, yg2, zg2);
        f[1] = interp3Dcomp_eval_f_precomputed(
            str1->c, n, x1, y1, z1,
            dxi, dx, dxi3, dx3,
            dyi, dy, dyi3, dy3,
            dzi, dz, dzi3, dz3,
            xg2, yg2, zg2);
        f[2] = interp3Dcomp_eval_f_precomputed(
            str2->c, n, x1, y1, z1,
            dxi, dx, dxi3, dx3,
            dyi, dy, dyi3, dy3,
            dzi, dz, dzi3, dz3,
            xg2, yg2, zg2);
    }

    return err;
}
GPU_DECLARE_TARGET_SIMD_UNIFORM(str0, str1, str2)
a5err interp3Dcomp_eval_df4_3(real f_df0[4],
                              real f_df1[4],
                              real f_df2[4],
                              interp3D_data* str0,
                              interp3D_data* str1,
                              interp3D_data* str2,
                              real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dexpl_eval_f(real* f, interp1D_data* str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dexpl_eval_f(real* f, interp2D_data* str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dexpl_eval_f(real* f, interp3D_data* str,
                          real x, real y, real z);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dcomp_eval_df(real* f_df, interp1D_data* str, real x);
DECLARE_TARGET_END
static inline a5err interp2Dcomp_eval_df(real* f_df, interp2D_data* str,
                                         real x, real y) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable. The -1 needed at exactly grid end. */
    int i_x   = (x - str->x_min) / str->x_grid - 1*(x==str->x_max);
    /* Normalized x coordinate in current cell */
    real dx     = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    /* Helper variables */
    real dx3    =  dx * (dx*dx - 1.0);
    real dx3dx  = 3*dx*dx - 1;
    real dxi    = 1.0 - dx;
    real dxi3   = dxi * (dxi*dxi - 1.0);
    real dxi3dx = -3*dxi*dxi + 1;
    real xg     = str->x_grid;
    real xg2    = xg*xg;
    real xgi    = 1.0/xg;

    /* Index for y variable. The -1 needed at exactly grid end. */
    int i_y   = (y - str->y_min) / str->y_grid - 1*(y==str->y_max);
    /* Normalized y coordinate in current cell */
    real dy     = ( y - (str->y_min + i_y*str->y_grid) ) / str->y_grid;
    /* Helper variables */
    real dy3    =  dy * (dy*dy - 1.0);
    real dy3dy  = 3*dy*dy - 1;
    real dyi    = 1.0 - dy;
    real dyi3   = dyi * (dyi*dyi - 1.0);
    real dyi3dy = -3*dyi*dyi + 1;
    real yg     = str->y_grid;
    real yg2    = yg*yg;
    real ygi    = 1.0/yg;

    int n  = i_y*str->n_x*4+i_x*4; /**< Index jump to cell       */
    int x1 = 4;                    /**< Index jump one x forward */
    int y1 = str->n_x*4;           /**< Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the domain. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max) ) {
        err = 1;
    }

    if(!err) {
        /* Preload all 16 coefficients (4 corners x 4 types = 2 cache lines).
           Avoids repeated reload after each f_df[] write due to pointer aliasing:
           the compiler cannot prove f_df and str->c don't overlap. */
        real c0   = str->c[n],       c1   = str->c[n+1];
        real c2   = str->c[n+2],     c3   = str->c[n+3];
        real cX0  = str->c[n+x1],    cX1  = str->c[n+x1+1];
        real cX2  = str->c[n+x1+2],  cX3  = str->c[n+x1+3];
        real cY0  = str->c[n+y1],    cY1  = str->c[n+y1+1];
        real cY2  = str->c[n+y1+2],  cY3  = str->c[n+y1+3];
        real cYX0 = str->c[n+y1+x1], cYX1 = str->c[n+y1+x1+1];
        real cYX2 = str->c[n+y1+x1+2], cYX3 = str->c[n+y1+x1+3];

        /* f */
        f_df[0] = (
            dxi*(dyi*c0+dy*cY0)
            +dx*(dyi*cX0+dy*cYX0))
            +(xg2/6)*(
                dxi3*(dyi*c1 + dy*cY1)
                +dx3*(dyi*cX1 + dy*cYX1))
            +(yg2/6)*(
                dxi*(dyi3*c2+dy3*cY2)
                +dx*(dyi3*cX2+dy3*cYX2))
            +(xg2*yg2/36)*(
                dxi3*(dyi3*c3+dy3*cY3)
                +dx3*(dyi3*cX3+dy3*cYX3));

        /* df/dx */
        f_df[1] = xgi*(
            -(dyi*c0  +dy*cY0)
            +(dyi*cX0+dy*cYX0))
            +(xg/6)*(
                dxi3dx*(dyi*c1  +dy*cY1)
                +dx3dx*(dyi*cX1+dy*cYX1))
            +(xgi*yg2/6)*(
                -(dyi3*c2  +dy3*cY2)
                +(dyi3*cX2+dy3*cYX2))
            +(xg*yg2/36)*(
                dxi3dx*(dyi3*c3  +dy3*cY3)
                +dx3dx*(dyi3*cX3+dy3*cYX3));

        /* df/dy */
        f_df[2] = ygi*(
            dxi*(-c0  +cY0)
            +dx*(-cX0+cYX0))
            +(xg2*ygi/6)*(
                dxi3*(-c1  +cY1)
                +dx3*(-cX1+cYX1))
            +(yg/6)*(
                dxi*(dyi3dy*c2  +dy3dy*cY2)
                +dx*(dyi3dy*cX2+dy3dy*cYX2))
            +(xg2*yg/36)*(
                dxi3*(dyi3dy*c3  +dy3dy*cY3)
                +dx3*(dyi3dy*cX3+dy3dy*cYX3));

        /* d2f/dx2 */
        f_df[3] = (
            dxi*(dyi*c1  +dy*cY1)
            +dx*(dyi*cX1+dy*cYX1))
            +(yg2/6)*(
                dxi*(dyi3*c3  +dy3*cY3)
                +dx*(dyi3*cX3+dy3*cYX3));

        /* d2f/dy2 */
        f_df[4] = (
              dxi*(dyi*c2  +dy*cY2)
              +dx*(dyi*cX2+dy*cYX2))
        +xg2/6*(
            dxi3*(dyi*c3  +dy*cY3)
            +dx3*(dyi*cX3+dy*cYX3));

        /* d2f/dydx */
        f_df[5] = xgi*ygi*(
            c0  -cY0
            -cX0+cYX0)
            +(xg/6*ygi)*(
                dxi3dx*(-c1  +cY1)
                +dx3dx*(-cX1+cYX1))
            +(xgi/6*yg)*(
                -(dyi3dy*c2  +dy3dy*cY2)
                +(dyi3dy*cX2+dy3dy*cYX2))
            +(xg*yg/36)*(
                dxi3dx*(dyi3dy*c3  +dy3dy*cY3)
                +dx3dx*(dyi3dy*cX3+dy3dy*cYX3));
    }

    return err;
}
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dcomp_eval_df(real* f_df, interp3D_data* str,
                           real x, real y, real z);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dcomp_eval_df4(real* f_df, interp3D_data* str,
                            real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dexpl_eval_df(real* f_df, interp1D_data* str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dexpl_eval_df(real* f_df, interp2D_data* str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dexpl_eval_df(real* f_df, interp3D_data* str,
                           real x, real y, real z);
#endif
