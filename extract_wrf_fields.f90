! gfortran extract_wrf_fields.f90 $(nf-config --fflags) $(nf-config --flibs) -o extract_wrf_fields
program extract_wrf_fields
  use netcdf
  use rh_utils
  use mod_config
  implicit none

  !------------------ constants ------------------
  real, parameter :: Rd = 287.0
  real, parameter :: g  = 9.80665

  !------------------ NetCDF ---------------------
  integer :: ncid
  integer :: nx, ny, nz, nt, varid_cldfra
  integer :: did_we, did_sn, did_bt, did_time
  integer :: vid_psfc, vid_t2, vid_q2
  integer :: vid_u10, vid_v10, vid_hgt
  integer :: vid_qc, vid_qi
  integer :: retval

  !------------------ arrays ---------------------
  real, allocatable :: psfc(:,:), t2(:,:), q2(:,:)
  real, allocatable :: u10(:,:), v10(:,:), hgt(:,:)
  real, allocatable :: qcloud(:,:,:), qice(:,:,:)
  real(kind=4), allocatable :: cldfra4(:,:,:,:)   ! WRF native
  real(kind=8), allocatable :: cloud_col(:)

  !------------------ derived --------------------
  real :: pmsl ! mean sea level pressure (PSFC, T2, HGT, formula)
  real :: rh2 ! 2m relative humidity (PSFC, Q2, T2)
  real :: tv ! virtual temperature (T2, and convert K to C)
  real :: es ! saturation vapor pressure 
  real :: e  ! vapor pressure
  real :: clc ! cloud cover flag (QCLOUD, QICE)
  integer :: i, j, k, t, 
  real :: temperature, psfc_temp
  character(len=256) :: outfile

  ! Load configuration from namelist (if exists)
  call load_config()

  print *, "Opening NetCDF file:", trim(wrf_infile)

  !------------------ open file ------------------
  retval = nf90_open(wrf_infile, NF90_NOWRITE, ncid)
  if (retval /= nf90_noerr) stop nf90_strerror(retval)

  !------------------ dimensions -----------------
  retval = nf90_inq_dimid(ncid,"west_east",did_we)
  if (retval /= nf90_noerr) stop nf90_strerror(retval)

  retval = nf90_inq_dimid(ncid,"south_north",did_sn)
  if (retval /= nf90_noerr) stop nf90_strerror(retval)

  retval = nf90_inq_dimid(ncid,"bottom_top",did_bt)
  if (retval /= nf90_noerr) stop nf90_strerror(retval)

  retval = nf90_inq_dimid(ncid,"Time",did_time)
  if (retval /= nf90_noerr) stop nf90_strerror(retval)


  retval = nf90_inquire_dimension(ncid,did_we,len=nx)
  if (retval /= nf90_noerr) stop nf90_strerror(retval)
  retval = nf90_inquire_dimension(ncid,did_sn,len=ny)
  retval = nf90_inquire_dimension(ncid,did_bt,len=nz)
  retval = nf90_inquire_dimension(ncid,did_time,len=nt)

  !------------------ variables ------------------
  retval = nf90_inq_varid(ncid,"PSFC",vid_psfc)
  retval = nf90_inq_varid(ncid,"T2",vid_t2)
  retval = nf90_inq_varid(ncid,"Q2",vid_q2)
  retval = nf90_inq_varid(ncid,"U10",vid_u10)
  retval = nf90_inq_varid(ncid,"V10",vid_v10)
  retval = nf90_inq_varid(ncid,"HGT",vid_hgt)
  retval = nf90_inq_varid(ncid,"QCLOUD",vid_qc)
  retval = nf90_inq_varid(ncid,"QICE",vid_qi)

  !------------------ allocate -------------------
  allocate(psfc(nx,ny),t2(nx,ny),q2(nx,ny))
  allocate(u10(nx,ny),v10(nx,ny),hgt(nx,ny))
  allocate(qcloud(nx,ny,nz),qice(nx,ny,nz))
  allocate(cldfra4(nx, ny, nz, 1))
  allocate(cloud_col(nz))

  !================== TIME LOOP ==================
  do t = 1, nt

     retval = nf90_get_var(ncid,vid_psfc,psfc,start=(/1,1,t/),count=(/nx,ny,1/))
     retval = nf90_get_var(ncid,vid_t2,t2,start=(/1,1,t/),count=(/nx,ny,1/))
     retval = nf90_get_var(ncid,vid_q2,q2,start=(/1,1,t/),count=(/nx,ny,1/))
     retval = nf90_get_var(ncid,vid_u10,u10,start=(/1,1,t/),count=(/nx,ny,1/))
     retval = nf90_get_var(ncid,vid_v10,v10,start=(/1,1,t/),count=(/nx,ny,1/))
     retval = nf90_get_var(ncid,vid_qc,qcloud,start=(/1,1,1,t/),count=(/nx,ny,nz,1/))
     retval = nf90_get_var(ncid,vid_qi,qice,start=(/1,1,1,t/),count=(/nx,ny,nz,1/))

     write(outfile,"('pmsl_',I3.3,'.dat')") t
     !write(outfile,'(A,I3.3,A)') trim(pmsl_outfile), t, '.dat'
     open(10,file=trim(outfile),status="replace")

     write(10,'(A)') "#pmsl psfc u10 v10 t2 rh2 clc"

     do j=1,ny
       do i=1,nx
         !=====psfc_temp calculation=====
         psfc_temp = psfc(i,j)
         tv = t2(i,j)*(1.0 + 0.61*q2(i,j))
         !===== temperature calculation=====
         temperature = t2(i,j)-273.15_rk
         !=====pmsl calculation=====
         pmsl = psfc_temp*exp(g*hgt(i,j)/(Rd*tv))
         es = 611.2*exp(17.67*(temperature)/(t2(i,j)-29.65))
         e  = q2(i,j)*psfc(i,j)/(0.622+0.378*q2(i,j))
         ! calculate relative humidity in %
         !===== rh2 calculation=====
         ! rh2 = 100.0*e/es
         rh2 = rh_vs_q_p_and_t(real(q2(i,j), rk),&
                   real(t2(i,j), rk)-273.15_rk, &
                   real(psfc(i,j), rk))

         cloud_col(:) = real(cldfra4(i,j,1:nz,1), kind=8)
        !===== clc calculation=====
         clc = real(total_cloud_cover(cloud_col, nz))

         do k=1,nz
           if (qcloud(i,j,k)+qice(i,j,k) > 1e-6) clc = 1.0
         end do

         write(10,'(7F14.6)') pmsl, &
                              psfc_temp, &
                              u10(i,j), &
                              v10(i,j), &
                              temperature, &
                              rh2, &
                              clc
       end do
     end do
     close(10)
  end do

  retval = nf90_close(ncid)
  print *,"DONE"

  call wrf_grid_export()
  ! -----------------------------------------------
contains

  real function total_cloud_cover(cloudfra, n)
    implicit none
    integer, intent(in) :: n
    real(kind=8), intent(in)  :: cloudfra(n)
    real(kind=8), intent(out) :: cl_rand, cl_gh, cl_max

    integer :: k
    real(kind=8) :: prod1, prod2

    ! Maximum overlap
    cl_max = maxval(cloudfra)

    ! Random overlap
    prod1 = 1.0_8
    do k = 1, n
      prod1 = prod1 * (1.0_8 - cloudfra(k))
    end do
    cl_rand = 1.0_8 - prod1

    ! Generalized (maximum–random) overlap
    prod2 = 1.0_8 - cloudfra(1)
    do k = 2, n
      if (cloudfra(k-1) >= 1.0_8) then
          prod2 = 0.0_8
          exit
      end if
      prod2 = prod2 * (1.0_8 - max(cloudfra(k), cloudfra(k-1))) &
                      / (1.0_8 - cloudfra(k-1))
    end do
    cl_gh = 1.0_8 - prod2

  end function total_cloud_cover
  ! -----------------------------------------------

  subroutine wrf_grid_export
    integer :: ncid, varid_lon, varid_lat
    integer :: dimid_we, dimid_sn
    integer :: nx, ny
    integer :: retval

    integer, dimension(3) :: start, count

    real(kind=8), allocatable :: xlon(:,:), xlat(:,:)

    integer :: i, j

    call load_config()

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
    end subroutine wrf_grid_export

end program extract_wrf_fields
