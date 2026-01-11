
module mod_config
  implicit none
  private

  integer, parameter :: STRLEN = 256

  public :: wrf_infile, grid_outfile, ntimes1, load_config

  ! Buffer variables for namelist(fixed line length)
  character(len=STRLEN) :: infile_buf  = "wrfout_d01_2025-10-22_12:00:00"
  character(len=STRLEN) :: grid_outfile_buf = "grid.data"
  integer               :: ntimes1_buf = 28    ! by default

  ! Variables to be used in the program
  character(len=STRLEN) :: wrf_infile
  character(len=STRLEN) :: grid_outfile
  integer               :: ntimes1

  namelist /io_nml/ infile_buf, grid_outfile_buf, ntimes1_buf

contains

  subroutine load_config(nml_path)
    character(len=*), intent(in), optional :: nml_path
    integer :: u, ios
    character(len=STRLEN) :: path
    logical :: ex

    path = merge(trim(nml_path), "config.nml", present(nml_path))

    inquire(file=trim(path), exist=ex)
    if (.not. ex) then
       call finalize_values()
       return
    end if

    open(newunit=u, file=trim(path), status="old", action="read", iostat=ios)
    if (ios == 0) then
       read(u, nml=io_nml, iostat=ios)
       close(u)
    else
       ! if not able to open, use defaults
    end if

    call finalize_values()
  end subroutine load_config

  subroutine finalize_values()
    wrf_infile  = trim(infile_buf)
    grid_outfile = trim(grid_outfile_buf)
    ntimes1 = ntimes1_buf
  end subroutine finalize_values

end module mod_config
