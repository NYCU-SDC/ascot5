/**
 * @file linint1D.c
 * @brief Linear interpolation
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../math.h"
#include "../spline/interp.h"
#include "linint.h"

/**
 * @brief Initialize linear interpolation struct for scalar 1D data
 *
 * @param str pointer to struct to be initialized
 * @param c array where data is stored
 * @param n_x number of data points in the x direction
 * @param bc_x boundary condition for x axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 */
void linint1D_init(linint1D_data* str, real* c,
                   int n_x, int bc_x, real x_min, real x_max) {

    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );

    str->n_x    = n_x;
    str->bc_x   = bc_x;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->c      = c;
}

/* linint1D_eval_f is static inline in linint.h */
