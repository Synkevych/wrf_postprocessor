module mod_output
  implicit none

contains

  subroutine write_grid_ascii(filename, xlon, xlat, nx, ny)
    character(len=*), intent(in) :: filename
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in) :: xlon(nx, ny), xlat(nx, ny)

    integer(kind=4) :: i, j
    integer(kind=4) :: unit_out

    unit_out = 20
    open(unit=unit_out, file=trim(filename), status='replace', action='write')
    write(unit_out, '(A)') '# i j loni latj'
    do j = 1, ny
      do i = 1, nx
        write(unit_out, '(I0,1X,I0,1X,F0.6,1X,F0.6)') i, j, xlon(i, j), xlat(i, j)
      end do
    end do
    close(unit_out)
  end subroutine write_grid_ascii

  subroutine write_pmsl_ascii(filename, pmsl, psfc, u10, v10, t2, rh2, clc, nx, ny)
    character(len=*), intent(in) :: filename
    integer(kind=4), intent(in) :: nx, ny
    real(kind=8), intent(in) :: pmsl(nx, ny), psfc(nx, ny), u10(nx, ny), v10(nx, ny)
    real(kind=8), intent(in) :: t2(nx, ny), rh2(nx, ny), clc(nx, ny)

    integer(kind=4) :: i, j
    integer(kind=4) :: unit_out

    unit_out = 21
    open(unit=unit_out, file=trim(filename), status='replace', action='write')
    write(unit_out, '(A)') '# pmsl psfc u10 v10 t2 rh2 clc'
    do j = 1, ny
      do i = 1, nx
        write(unit_out, '(F14.6,1X,F14.6,1X,F12.6,1X,F12.6,1X,F12.6,1X,F12.6,1X,F12.6)') &
          pmsl(i, j), psfc(i, j), u10(i, j), v10(i, j), t2(i, j), rh2(i, j), clc(i, j)
      end do
    end do
    close(unit_out)
  end subroutine write_pmsl_ascii

end module mod_output
