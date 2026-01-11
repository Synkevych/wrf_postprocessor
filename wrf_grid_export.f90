program wrf_grid_export
  use netcdf
  use mod_config
  implicit none

  integer :: ncid, varid_lon, varid_lat
  integer :: dimid_we, dimid_sn
  integer :: nx, ny
  integer :: retval

  integer, dimension(3) :: start, count

  real(kind=8), allocatable :: xlon(:,:), xlat(:,:)

  integer :: i, j

  !character(len=*), parameter :: infile  = "test.nc"
  !character(len=8), parameter :: outfile = "grid.dat"
  
  ! Load configuration from namelist (if exists)
  call load_config("config.nml")

  print *, "Opening NetCDF file:", trim(wrf_infile)

  retval = nf90_open(wrf_infile, NF90_NOWRITE, ncid)
  print *, "wrf_infile", wrf_infile
  print *, "ncid", ncid
  if (retval /= nf90_noerr) then
     print *, nf90_strerror(retval)
     stop
  endif

  print *, "Reading dimensions..."

  retval = nf90_inq_dimid(ncid, "west_east", dimid_we)
  retval = nf90_inq_dimid(ncid, "south_north", dimid_sn)

  retval = nf90_inquire_dimension(ncid, dimid_we, len=nx)
  retval = nf90_inquire_dimension(ncid, dimid_sn, len=ny)

  print *, "Grid size:", nx, "x", ny

  allocate(xlon(nx, ny))
  allocate(xlat(nx, ny))

  retval = nf90_inq_varid(ncid, "XLONG", varid_lon)
  retval = nf90_inq_varid(ncid, "XLAT",  varid_lat)

  ! Read ONLY Time = 1
  start = (/1, 1, 1/)
  count = (/nx, ny, 1/)

  print *, "Reading XLONG..."
  retval = nf90_get_var(ncid, varid_lon, xlon, start=start, count=count)

  print *, "Reading XLAT..."
  retval = nf90_get_var(ncid, varid_lat, xlat, start=start, count=count)

  retval = nf90_close(ncid)

  print *, "Lon min/max:", minval(xlon), maxval(xlon)
  print *, "Lat min/max:", minval(xlat), maxval(xlat)

  open(unit=10, file=grid_outfile, status="replace")
  write(10,'(A)') "# i j lon lat"

  do j = 1, ny
     do i = 1, nx
        write(10,'(I6,1X,I6,1X,F12.6,1X,F12.6)') i, j, xlon(i,j), xlat(i,j)
     end do
  end do
  close(10)

  print *, "DONE"

end program wrf_grid_export
