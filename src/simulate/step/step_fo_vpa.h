/**
 * @file step_fo_vpa.h
 * @brief Header file for step_fo_vpa.c
 */
#ifndef STEP_FO_VPA_H
#define STEP_FO_VPA_H
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../boozer.h"
#include "../../mhd.h"
#include "../../particle.h"
#include "../../simulate.h"

void step_fo_vpa(particle_simd_fo* p, particle_simd_fo* p0, real* h,
                 B_field_data* Bdata, E_field_data* Edata,
                 sim_data* sim, real cputime, int fuse_endcond,
                 int aldforce, int reverse_time);
void step_fo_vpa_mhd(
    particle_simd_fo* p, particle_simd_fo* p0, real* h,
    B_field_data* Bdata, E_field_data* Edata,
    boozer_data* boozer, mhd_data* mhd,
    sim_data* sim, real cputime, int fuse_endcond,
    int aldforce, int reverse_time);

#endif
