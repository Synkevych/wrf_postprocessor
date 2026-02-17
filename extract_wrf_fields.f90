program extract_wrf_fields
  use mod_config
  use mod_netcdf_io
  use mod_output
  implicit none

  type(config_type) :: cfg
  integer(kind=4) :: t
  character(len=1024) :: outname
  real(kind=8), allocatable :: xlon(:, :), xlat(:, :)
  real(kind=8), allocatable :: pmsl(:, :), psfc(:, :), u10(:, :), v10(:, :), t2(:, :), rh2(:, :), clc(:, :)

  call read_config(cfg)
  call open_input_and_read_sizes(cfg)

  allocate(xlon(cfg%nx, cfg%ny), xlat(cfg%nx, cfg%ny))
  allocate(pmsl(cfg%nx, cfg%ny), psfc(cfg%nx, cfg%ny), u10(cfg%nx, cfg%ny), v10(cfg%nx, cfg%ny))
  allocate(t2(cfg%nx, cfg%ny), rh2(cfg%nx, cfg%ny), clc(cfg%nx, cfg%ny))

  call read_grid(cfg%ncid, cfg%nx, cfg%ny, xlon, xlat)
  call write_grid_ascii(trim(cfg%grid_outfile), xlon, xlat, cfg%nx, cfg%ny)

  do t = 1, cfg%nt_proc
    call read_surface_fields_t(cfg%ncid, cfg%nx, cfg%ny, t, psfc, u10, v10, t2, rh2, clc, pmsl)
    write(outname, '(A,I3.3,A)') trim(cfg%pmsl_outfile), t, '.dat'
    call write_pmsl_ascii(trim(outname), pmsl, psfc, u10, v10, t2, rh2, clc, cfg%nx, cfg%ny)
  end do

  call close_input(cfg)

end program extract_wrf_fields
