/**
 * @file linint.h
 * @brief Linear interpolation library
 *
 * Linear interpolation interpolates data given on a uniform grid.
 * Each axis may have either natural or periodic boundary condition.
 *
 * The interpolant does not need any initialization on the host. After
 * offloading is done, call linintXD_init() which assigns the values and
 * a pointer to the data to a linint struct.
 */
#ifndef LININT_H
#define LININT_H
#include <math.h>
#include "../ascot5.h"
#include "../offload.h"
#include "../spline/interp.h"

/**
 * @brief 1D interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real* c;     /**< pointer to array with interpolant values       */
} linint1D_data;

/**
 * @brief 2D interpolation struct.
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
    real* c;     /**< pointer to array with interpolant values       */
} linint2D_data;

/**
 * @brief 3D interpolation struct.
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
    real* c;     /**< pointer to array with interpolant values       */
} linint3D_data;

void linint1D_init(linint1D_data* str, real* c,
                   int n_x, int bc_x,
                   real x_min, real x_max);

void linint2D_init(linint2D_data* str, real* c,
                   int n_x, int n_y, int bc_x, int bc_y,
                   real x_min, real x_max,
                   real y_min, real y_max);

void linint3D_init(linint3D_data* str, real* c,
                   int n_x, int n_y, int n_z,
                   int bc_x, int bc_y, int bc_z,
                   real x_min, real x_max,
                   real y_min, real y_max,
                   real z_min, real z_max);

static inline int linint1D_eval_f(real* f, linint1D_data* str, real x) {

    /* Make sure periodic coordinates are within [max, min] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* index for x variable */
    int i_x   = (x-str->x_min) / str->x_grid;
    /**< Normalized x coordinate in current cell */
    real dx = ( x - (str->x_min + i_x*str->x_grid)) / str->x_grid;

    int x1 = 1;   /* Index jump one x forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }

    if(!err) {
        *f = str->c[i_x]*(1 - dx) + str->c[i_x+x1]*dx;
    }
    return err;
}

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
int linint2D_eval_f(real* f, linint2D_data* str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
int linint3D_eval_f(real* f, linint3D_data* str,
                    real x, real y, real z);
DECLARE_TARGET_END
#endif
