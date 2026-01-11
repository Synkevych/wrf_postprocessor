#!/usr/bin/env bash
#gfortran mod_config.f90 wrf_grid_export.f90 \
#  $(nf-config --fflags) $(nf-config --flibs) -o wrf_grid_export

gfortran mod_config.f90 rh_utils.f90 extract_wrf_fields.f90 \
  $(nf-config --fflags) $(nf-config --flibs) -o extract_wrf_fields
