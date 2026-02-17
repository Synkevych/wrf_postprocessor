module mod_config
  use netcdf
  implicit none

  type :: config_type
    character(len=512) :: infile = 'test.nc'
    character(len=512) :: grid_outfile = 'grid.dat'
    character(len=512) :: pmsl_outfile = 'pmsl_'
    integer(kind=4) :: ntimes1 = 28
    integer(kind=4) :: nx = 0
    integer(kind=4) :: ny = 0
    integer(kind=4) :: nt_file = 0
    integer(kind=4) :: nt_proc = 0
    integer(kind=4) :: ncid = -1
  end type config_type

contains

  subroutine read_config(cfg, filename)
    type(config_type), intent(inout) :: cfg
    character(len=*), intent(in), optional :: filename
    character(len=512) :: nml_file
    character(len=512) :: infile
    character(len=512) :: grid_outfile
    character(len=512) :: pmsl_outfile
    integer(kind=4) :: ntimes1
    integer(kind=4) :: unit_nml
    integer(kind=4) :: ios
    logical :: exists
    namelist /io_nml/ infile, grid_outfile, pmsl_outfile, ntimes1

    nml_file = 'config.nml'
    if (present(filename)) nml_file = trim(filename)

    infile = cfg%infile
    grid_outfile = cfg%grid_outfile
    pmsl_outfile = cfg%pmsl_outfile
    ntimes1 = cfg%ntimes1

    inquire(file=trim(nml_file), exist=exists)
    if (.not. exists) then
      call stop_with_message('Missing namelist file: '//trim(nml_file))
    end if

    unit_nml = 10
    open(unit=unit_nml, file=trim(nml_file), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      call stop_with_message('Cannot open namelist file: '//trim(nml_file))
    end if

    read(unit_nml, nml=io_nml, iostat=ios)
    close(unit_nml)
    if (ios /= 0) then
      call stop_with_message('Cannot read namelist group io_nml in: '//trim(nml_file))
    end if

    cfg%infile = trim(infile)
    cfg%grid_outfile = trim(grid_outfile)
    cfg%pmsl_outfile = trim(pmsl_outfile)
    cfg%ntimes1 = ntimes1
  end subroutine read_config

  subroutine open_input_and_read_sizes(cfg)
    type(config_type), intent(inout) :: cfg
    integer(kind=4) :: rc
    integer(kind=4) :: dimid

    rc = nf90_open(trim(cfg%infile), nf90_nowrite, cfg%ncid)
    if (rc /= nf90_noerr) then
      call stop_with_message('nf90_open failed: '//trim(nf90_strerror(rc)))
    end if

    rc = nf90_inq_dimid(cfg%ncid, 'west_east', dimid)
    if (rc /= nf90_noerr) call stop_with_message('Missing dimension west_east')
    rc = nf90_inquire_dimension(cfg%ncid, dimid, len=cfg%nx)
    if (rc /= nf90_noerr) call stop_with_message('Failed to read west_east length')

    rc = nf90_inq_dimid(cfg%ncid, 'south_north', dimid)
    if (rc /= nf90_noerr) call stop_with_message('Missing dimension south_north')
    rc = nf90_inquire_dimension(cfg%ncid, dimid, len=cfg%ny)
    if (rc /= nf90_noerr) call stop_with_message('Failed to read south_north length')

    rc = nf90_inq_dimid(cfg%ncid, 'Time', dimid)
    if (rc /= nf90_noerr) call stop_with_message('Missing dimension Time')
    rc = nf90_inquire_dimension(cfg%ncid, dimid, len=cfg%nt_file)
    if (rc /= nf90_noerr) call stop_with_message('Failed to read Time length')

    cfg%nt_proc = min(cfg%ntimes1, cfg%nt_file)
  end subroutine open_input_and_read_sizes

  subroutine close_input(cfg)
    type(config_type), intent(inout) :: cfg
    integer(kind=4) :: rc

    if (cfg%ncid >= 0) then
      rc = nf90_close(cfg%ncid)
      if (rc /= nf90_noerr) call stop_with_message('nf90_close failed: '//trim(nf90_strerror(rc)))
      cfg%ncid = -1
    end if
  end subroutine close_input

  subroutine stop_with_message(msg)
    character(len=*), intent(in) :: msg
    write(*,'(A)') trim(msg)
    stop 1
  end subroutine stop_with_message

end module mod_config
