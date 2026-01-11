program wrf_driver
  use mod_config
  !use extract_wrf_fields
  !use wrf_grid_export

  implicit none

  ! 1) Load configuration from file
  call load_config("config.nml")

  ! 2) Run WRF grid building and export
  !call build_grid()   ! create grid.dat file
  !call export_wrf()   ! create pmsl_*.dat files

  print *, "wrf_driver: all done."
end program wrf_driver
