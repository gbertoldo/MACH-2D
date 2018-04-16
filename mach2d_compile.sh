#!/bin/bash


gfortran -O3 -o mach2d.x \
    depp_interface.f90 \
    mach2d_data.f90 \
    msi2d9-v03.f90 \
    msi2d5-v03.f90 \
    tdma2d9.f90 \
    tdma2d5.f90 \
    mach2d_solvers.f90 \
    mach2d_bc5d.f90 \
    mach2d_bc9d.f90 \
    mach2d_coef.f90 \
    thompson2d_hyperbolic.f90 \
    mach2d_grid.f90 \
    mach2d_postp.f90 \
    mach2d_thermo.f90 \
    newton2coef.f90 \
    spline.f90 \
    geometry_splwsp.f90 \
    geometry_splwep.f90 \
    geometry_ncspline.f90\
    geometry_d2cspline.f90\
    geometry_spline_s01.f90\
    mach2d_user.f90 \
    mach2d_main.f90

rm *.mod
