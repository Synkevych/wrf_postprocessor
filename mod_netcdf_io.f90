module mod_netcdf_io
  use netcdf
  use, intrinsic :: ieee_arithmetic
  implicit none

contains

  subroutine check_nc(rc, where)
    integer(kind=4), intent(in) :: rc
    character(len=*), intent(in) :: where
    if (rc /= nf90_noerr) then
      write(*,'(A)') trim(where)//': '//trim(nf90_strerror(rc))
      stop 1
    end if
  end subroutine check_nc

  logical function var_exists(ncid, varname)
    integer(kind=4), intent(in) :: ncid
    character(len=*), intent(in) :: varname
    integer(kind=4) :: rc
    integer(kind=4) :: varid

    rc = nf90_inq_varid(ncid, trim(varname), varid)
    var_exists = (rc == nf90_noerr)
  end function var_exists

  subroutine read_grid(ncid, nx, ny, xlon, xlat)
    integer(kind=4), intent(in) :: ncid
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(out) :: xlon(nx, ny)
    real(kind=8), intent(out) :: xlat(nx, ny)

    call read_2d_at_time(ncid, 'XLONG', nx, ny, 1_4, xlon)
    call read_2d_at_time(ncid, 'XLAT', nx, ny, 1_4, xlat)
  end subroutine read_grid

  subroutine read_surface_fields_t(ncid, nx, ny, tidx, psfc, u10, v10, t2c, rh2, clc, pmsl)
    integer(kind=4), intent(in) :: ncid
    integer(kind=4), intent(in) :: nx, ny
    integer(kind=4), intent(in) :: tidx
    real(kind=8), intent(out) :: psfc(nx, ny)
    real(kind=8), intent(out) :: u10(nx, ny)
    real(kind=8), intent(out) :: v10(nx, ny)
    real(kind=8), intent(out) :: t2c(nx, ny)
    real(kind=8), intent(out) :: rh2(nx, ny)
    real(kind=8), intent(out) :: clc(nx, ny)
    real(kind=8), intent(out) :: pmsl(nx, ny)

    real(kind=8), allocatable :: t2k(:, :)
    real(kind=8), allocatable :: q2(:, :)

    allocate(t2k(nx, ny), q2(nx, ny))

    call read_2d_at_time(ncid, 'PSFC', nx, ny, tidx, psfc)
    call read_2d_at_time(ncid, 'U10', nx, ny, tidx, u10)
    call read_2d_at_time(ncid, 'V10', nx, ny, tidx, v10)
    call read_2d_at_time(ncid, 'T2', nx, ny, tidx, t2k)
    t2c = t2k - 273.16_8

    if (var_exists(ncid, 'RH2')) then
      call read_2d_at_time(ncid, 'RH2', nx, ny, tidx, rh2)
    else
      call read_2d_at_time(ncid, 'Q2', nx, ny, tidx, q2)
      call rh_from_q2(nx, ny, q2, t2k, psfc, rh2)
    end if

    if (var_exists(ncid, 'CLDFRA')) then
      call read_cldfra_max_at_time(ncid, nx, ny, tidx, clc)
      clc = clc * 100.0_8
    else
      clc = 0.0_8
    end if

    if (var_exists(ncid, 'PMSL')) then
      call read_2d_at_time(ncid, 'PMSL', nx, ny, tidx, pmsl)
    else
      pmsl = psfc
    end if

    deallocate(t2k, q2)
  end subroutine read_surface_fields_t

  subroutine read_2d_at_time(ncid, varname, nx, ny, tidx, field)
    integer(kind=4), intent(in) :: ncid
    character(len=*), intent(in) :: varname
    integer(kind=4), intent(in) :: nx, ny
    integer(kind=4), intent(in) :: tidx
    real(kind=8), intent(out) :: field(nx, ny)

    integer(kind=4) :: rc
    integer(kind=4) :: varid
    integer(kind=4) :: ndims
    integer(kind=4) :: dimids(nf90_max_var_dims)
    integer(kind=4) :: start2(2), count2(2)
    integer(kind=4) :: start3(3), count3(3)

    rc = nf90_inq_varid(ncid, trim(varname), varid)
    call check_nc(rc, 'nf90_inq_varid('//trim(varname)//')')

    rc = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    call check_nc(rc, 'nf90_inquire_variable('//trim(varname)//')')

    select case (ndims)
    case (2)
      start2 = (/1_4, 1_4/)
      count2 = (/nx, ny/)
      rc = nf90_get_var(ncid, varid, field, start=start2, count=count2)
      call check_nc(rc, 'nf90_get_var('//trim(varname)//',2D)')
    case (3)
      start3 = (/1_4, 1_4, tidx/)
      count3 = (/nx, ny, 1_4/)
      rc = nf90_get_var(ncid, varid, field, start=start3, count=count3)
      call check_nc(rc, 'nf90_get_var('//trim(varname)//',3D@t)')
    case default
      write(*,'(A,I0,A)') 'Unsupported rank for ', ndims, ': '//trim(varname)
      stop 1
    end select

    call apply_fillvalue_2d(ncid, varid, field, nx, ny)
  end subroutine read_2d_at_time

  subroutine read_cldfra_max_at_time(ncid, nx, ny, tidx, clc_max)
    integer(kind=4), intent(in) :: ncid
    integer(kind=4), intent(in) :: nx, ny
    integer(kind=4), intent(in) :: tidx
    real(kind=8), intent(out) :: clc_max(nx, ny)

    integer(kind=4) :: rc
    integer(kind=4) :: varid
    integer(kind=4) :: ndims
    integer(kind=4) :: dimids(nf90_max_var_dims)
    integer(kind=4) :: dimlen
    integer(kind=4) :: start4(4), count4(4)
    integer(kind=4) :: i, j, k
    integer(kind=4) :: nz
    logical :: has_valid
    real(kind=8), allocatable :: cld(:, :, :, :)

    rc = nf90_inq_varid(ncid, 'CLDFRA', varid)
    call check_nc(rc, 'nf90_inq_varid(CLDFRA)')

    rc = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    call check_nc(rc, 'nf90_inquire_variable(CLDFRA)')

    if (ndims /= 4) then
      write(*,'(A,I0)') 'CLDFRA must have rank 4, got ', ndims
      stop 1
    end if

    rc = nf90_inquire_dimension(ncid, dimids(3), len=dimlen)
    call check_nc(rc, 'nf90_inquire_dimension(bottom_top)')
    nz = dimlen

    allocate(cld(nx, ny, nz, 1))
    start4 = (/1_4, 1_4, 1_4, tidx/)
    count4 = (/nx, ny, nz, 1_4/)
    rc = nf90_get_var(ncid, varid, cld, start=start4, count=count4)
    call check_nc(rc, 'nf90_get_var(CLDFRA,4D@t)')

    call apply_fillvalue_4d_single_t(ncid, varid, cld, nx, ny, nz)

    do j = 1, ny
      do i = 1, nx
        clc_max(i, j) = ieee_value(0.0_8, ieee_quiet_nan)
        has_valid = .false.
        do k = 1, nz
          if (.not. ieee_is_nan(cld(i, j, k, 1))) then
            if (.not. has_valid) then
              clc_max(i, j) = cld(i, j, k, 1)
              has_valid = .true.
            else
              clc_max(i, j) = max(clc_max(i, j), cld(i, j, k, 1))
            end if
          end if
        end do
      end do
    end do

    deallocate(cld)
  end subroutine read_cldfra_max_at_time

  subroutine rh_from_q2(nx, ny, q2, t2k, psfc, rh2)
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in) :: q2(nx, ny)
    real(kind=8), intent(in) :: t2k(nx, ny)
    real(kind=8), intent(in) :: psfc(nx, ny)
    real(kind=8), intent(out) :: rh2(nx, ny)

    real(kind=8), parameter :: eps = 0.622_8
    real(kind=8) :: es
    real(kind=8) :: e
    real(kind=8) :: q
    integer(kind=4) :: i, j

    do j = 1, ny
      do i = 1, nx
        if (ieee_is_nan(q2(i, j)) .or. ieee_is_nan(t2k(i, j)) .or. ieee_is_nan(psfc(i, j))) then
          rh2(i, j) = ieee_value(0.0_8, ieee_quiet_nan)
        else
          q = q2(i, j)
          es = esonT(t2k(i, j))
          if (psfc(i, j) <= 0.0_8 .or. es <= 0.0_8) then
            rh2(i, j) = ieee_value(0.0_8, ieee_quiet_nan)
          else
            e = q * psfc(i, j) / (eps + (1.0_8 - eps) * q)
            rh2(i, j) = 100.0_8 * e / es
          end if
        end if
      end do
    end do
  end subroutine rh_from_q2

  real(kind=8) function esonT(tk)
    real(kind=8), intent(in) :: tk
    real(kind=8) :: log10e

    log10e = 10.79574_8 * (1.0_8 - 273.16_8 / tk) - 5.028_8 * log10(tk / 273.16_8) + &
             1.50475e-4_8 * (1.0_8 - 10.0_8**(-8.2969_8 * (tk / 273.16_8 - 1.0_8))) + &
             0.42873e-3_8 * (10.0_8**(4.76955_8 * (1.0_8 - 273.16_8 / tk)) - 1.0_8) + &
             0.78614_8 + 2.0_8
    esonT = 10.0_8**log10e
  end function esonT

  subroutine apply_fillvalue_2d(ncid, varid, a, nx, ny)
    integer(kind=4), intent(in) :: ncid, varid
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(inout) :: a(nx, ny)

    integer(kind=4) :: rc
    integer(kind=4) :: i, j
    real(kind=8) :: fill
    real(kind=8) :: tol

    rc = nf90_get_att(ncid, varid, '_FillValue', fill)
    if (rc /= nf90_noerr) return

    tol = max(1.0d-12, abs(fill) * 1.0d-10)
    do j = 1, ny
      do i = 1, nx
        if (abs(a(i, j) - fill) <= tol) then
          a(i, j) = ieee_value(0.0_8, ieee_quiet_nan)
        end if
      end do
    end do
  end subroutine apply_fillvalue_2d

  subroutine apply_fillvalue_4d_single_t(ncid, varid, a, nx, ny, nz)
    integer(kind=4), intent(in) :: ncid, varid
    integer(kind=4), intent(in) :: nx, ny, nz
    real(kind=8), intent(inout) :: a(nx, ny, nz, 1)

    integer(kind=4) :: rc
    integer(kind=4) :: i, j, k
    real(kind=8) :: fill
    real(kind=8) :: tol

    rc = nf90_get_att(ncid, varid, '_FillValue', fill)
    if (rc /= nf90_noerr) return

    tol = max(1.0d-12, abs(fill) * 1.0d-10)
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (abs(a(i, j, k, 1) - fill) <= tol) then
            a(i, j, k, 1) = ieee_value(0.0_8, ieee_quiet_nan)
          end if
        end do
      end do
    end do
  end subroutine apply_fillvalue_4d_single_t

end module mod_netcdf_io
