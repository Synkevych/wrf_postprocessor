#!/usr/bin/env bash
set -euo pipefail

FFLAGS="$(nf-config --fflags) -O2"
FLIBS="$(nf-config --flibs)"

rm -f ./*.o ./*.mod extract_wrf_fields validate_fields

gfortran ${FFLAGS} -c mod_config.f90
gfortran ${FFLAGS} -c mod_netcdf_io.f90
gfortran ${FFLAGS} -c mod_output.f90
gfortran ${FFLAGS} -c extract_wrf_fields.f90
gfortran ${FFLAGS} -c validate_fields.f90

gfortran ${FFLAGS} -o extract_wrf_fields mod_config.o mod_netcdf_io.o mod_output.o extract_wrf_fields.o ${FLIBS}
gfortran ${FFLAGS} -o validate_fields mod_config.o mod_netcdf_io.o validate_fields.o ${FLIBS}
