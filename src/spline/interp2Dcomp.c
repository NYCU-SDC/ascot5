/**
 * @file interp2Dcomp.c
 * @brief Bicubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../math.h"
#include "interp.h"
#include "spline.h"

/**
 * @brief Calculate bicubic spline interpolation coefficients for scalar 2D data
 *
 * This function calculates the bicubic spline interpolation coefficients and
 * stores them in a pre-allocated array. Compact cofficients are calculated.
 *
 * For each data point four coefficients are stored for spline-interpolation.
 *
 * @param c allocated array of length n_y*n_x*4 to store the coefficients
 * @param f 2D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis (0) natural (1) periodic
 * @param bc_y boundary condition for y axis (0) natural (1) periodic
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 *
 * @return zero if initialization succeeded
 */
int interp2Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max) {

    /* Check boundary conditions and calculate grid intervals. Grid intervals
       needed because we use normalized grid intervals. For periodic boundary
       condition, grid maximum value and the last data point are not the same.
       Take this into account in grid intervals. */
    real x_grid, y_grid;
    if(bc_x == PERIODICBC || bc_x == NATURALBC) {
        x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    }
    else {
        return 1;
    }

    if(bc_y == PERIODICBC || bc_y == NATURALBC) {
        y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );
    }
    else {
        return 1;
    }

    /* Allocate helper quantities */
    real* f_x = malloc(n_x*sizeof(real));
    real* f_y = malloc(n_y*sizeof(real));
    real* c_x = malloc(n_x*NSIZE_COMP1D*sizeof(real));
    real* c_y = malloc(n_y*NSIZE_COMP1D*sizeof(real));

    if(f_x == NULL || f_y == NULL || c_x == NULL || c_y == NULL) {
        return 1;
    }

    /* Calculate bicubic spline surface coefficients, i.e second derivatives.
       For each grid cell (i_x, i_y), there are four coefficients:
       [f, fxx, fyy, fxxyy]. Note how we account for normalized grid. */

    /* Cubic spline along x for each y, using f values to get fxx */
    for(int i_y=0; i_y<n_y; i_y++) {
        /* fxx */
        for(int i_x=0; i_x<n_x; i_x++) {
            f_x[i_x] = f[i_y*n_x+i_x];
        }
        splinecomp(f_x, n_x, bc_x, c_x);
        for(int i_x=0; i_x<n_x; i_x++) {
            c[i_y*n_x*4 + i_x*4    ] = c_x[i_x*2];
            c[i_y*n_x*4 + i_x*4 + 1] = c_x[i_x*2+1] / (x_grid*x_grid);
        }
    }

    /* Two cubic splines along y for each x, using f and fxx to get fyy and
       fxxyy */
    for(int i_x=0; i_x<n_x; i_x++) {

        /* fyy */
        for(int i_y=0; i_y<n_y; i_y++) {
            f_y[i_y] =  f[i_y*n_x + i_x];
        }
        splinecomp(f_y, n_y, bc_y, c_y);
        for(int i_y=0; i_y<n_y; i_y++) {
            c[i_y*n_x*4+i_x*4+2] = c_y[i_y*2+1]/(y_grid*y_grid);
        }

        /* fxxyy */
        for(int i_y=0; i_y<n_y; i_y++) {
            f_y[i_y] =  c[i_y*n_x*4 + i_x*4 + 1];
        }
        splinecomp(f_y, n_y, bc_y, c_y);
        for(int i_y=0; i_y<n_y; i_y++) {
            c[i_y*n_x*4 + i_x*4 + 3] = c_y[i_y*2 + 1] / (y_grid*y_grid);
        }
    }

    /* Free allocated memory */
    free(f_x);
    free(f_y);
    free(c_x);
    free(c_y);

    return 0;
}

/**
 * @brief Initialize a bicubic spline
 *
 * @param str pointer to spline to be initialized
 * @param c array where coefficients are stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 */
void interp2Dcomp_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max) {

    /* Calculate grid intervals. For periodic boundary condition, grid maximum
       value and the last data point are not the same. Take this into account
       in grid intervals. */
    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );

    /* Initialize the interp2D_data struct */
    str->n_x    = n_x;
    str->n_y    = n_y;
    str->bc_x   = bc_x;
    str->bc_y   = bc_y;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->c      = c;
}

/**
 * @brief Set up splines to interpolate 2D scalar data.
 *
 * This function is equivalent to calling both `interp2Dcomp_init_coeff` and
 * `interp2Dcomp_init_spline`.
 *
 * @param str pointer to spline to be initialized
 * @param f 2D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 *
 * @return zero if initialization succeeded
 */
int interp2Dcomp_setup(interp2D_data* str, real* f,
                       int n_x, int n_y, int bc_x, int bc_y,
                       real x_min, real x_max, real y_min, real y_max) {
    real* c = (real*) malloc(n_y*n_x*NSIZE_COMP2D*sizeof(real));
    int err = interp2Dcomp_init_coeff(c, f, n_x, n_y, bc_x, bc_y, x_min, x_max,
                                      y_min, y_max);
    if(err) {
        return err;
    }
    interp2Dcomp_init_spline(str, c, n_x, n_y, bc_x, bc_y, x_min, x_max,
                             y_min, y_max);
    return 0;
}

