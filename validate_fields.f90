program validate_fields
  use mod_config
  use mod_netcdf_io
  use, intrinsic :: ieee_arithmetic
  implicit none

  type(config_type) :: cfg
  real(kind=8), allocatable :: pmsl(:, :), psfc(:, :), u10(:, :), v10(:, :), t2(:, :), rh2(:, :), clc(:, :)

  call read_config(cfg)
  call open_input_and_read_sizes(cfg)

  allocate(pmsl(cfg%nx, cfg%ny), psfc(cfg%nx, cfg%ny), u10(cfg%nx, cfg%ny), v10(cfg%nx, cfg%ny))
  allocate(t2(cfg%nx, cfg%ny), rh2(cfg%nx, cfg%ny), clc(cfg%nx, cfg%ny))

  call read_surface_fields_t(cfg%ncid, cfg%nx, cfg%ny, 1_4, psfc, u10, v10, t2, rh2, clc, pmsl)

  call print_minmax('PMSL', pmsl)
  call print_minmax('PSFC', psfc)
  call print_minmax('U10 ', u10)
  call print_minmax('V10 ', v10)
  call print_minmax('T2C ', t2)
  call print_minmax('RH2 ', rh2)
  call print_minmax('CLC ', clc)

  call close_input(cfg)

contains

  subroutine print_minmax(name, a)
    character(len=*), intent(in) :: name
    real(kind=8), intent(in) :: a(:, :)

    integer(kind=4) :: i, j
    real(kind=8) :: amin, amax
    logical :: has_val

    has_val = .false.
    amin = 0.0_8
    amax = 0.0_8

    do j = 1, size(a, 2)
      do i = 1, size(a, 1)
        if (.not. ieee_is_nan(a(i, j))) then
          if (.not. has_val) then
            amin = a(i, j)
            amax = a(i, j)
            has_val = .true.
          else
            amin = min(amin, a(i, j))
            amax = max(amax, a(i, j))
          end if
        end if
      end do
    end do

    if (has_val) then
      write(*, '(A,1X,ES14.6,1X,ES14.6)') trim(name), amin, amax
    else
      write(*, '(A,1X,A)') trim(name), 'NaN NaN'
    end if
  end subroutine print_minmax

end program validate_fields
