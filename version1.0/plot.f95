module plot
  use plplot
  implicit none
  private

  public plot_init, plot_close, plot_points

contains

  subroutine plot_init(boxSize)
    real(8), intent(in) :: boxSize
    call plsdev("xcairo")
    call plinit()

    !call plscol0(0, 255, 255, 255)  ! white
    !call plscol0(1, 255, 0, 0)      ! red
    !call plscol0(2, 0, 255, 0)      ! green
    !call plscol0(3, 0, 0, 255)      ! blue
    !call plscol0(4, 255, 0, 255)    ! magenta
    !call plscol0(5, 0, 255, 255)    ! cyan
    !call plscol0(6, 255, 255, 0)    ! yellow
    !call plscol0(7, 0, 0, 0)        ! black
    !call plscol0(8, 255, 77, 0)     ! orange
    !call plscol0(9, 128, 128, 128)  ! gray

    call pladv(0)
    call plvpor(0d0, 1d0, 0d0, 1d0)
    call plwind(-1d0, 1d0, -2d0/3, 4d0/3)
    call plw3d(1d0, 1d0, 1d0, 0d0, boxSize, 0d0, boxSize, 0d0, boxSize, 45d0, &
    -45d0)
  end subroutine plot_init

  subroutine plot_close()
    call plspause(.false.)
    call plend()
  end subroutine plot_close

  subroutine plot_points(xyz)
    real(8), intent(in) :: xyz(:,:)

    call plclear()
    call plcol0(1)
    call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", 0d0, 0, "bcnmstuv", &
      "z", 0d0, 0)
    call plcol0(2)
    call plpoin3(xyz(1, :), xyz(2, :), xyz(3, :), 4)
    call plflush()
  end subroutine plot_points

end module plot
