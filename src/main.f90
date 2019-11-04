program jump

    use jump_integration
    use eispack
    
    implicit none

    real (kind = 8) interval, x_test
    integer (kind = 4) row, col, test

    real (kind = 8), dimension(:, :), allocatable :: realMatrix, imagMatrix

    real (kind = 8), dimension(:), allocatable :: w
    real (kind = 8), dimension(:, :), allocatable :: zi, zr

    integer ( kind = 4 ) ierr
    real (kind = 8), dimension(:), allocatable ::  points

    complex (kind = 8) tempVar

    integer         :: MODEL_N
    real(kind = 8)  :: MODEL_P_CUTOFF
    real(kind = 8)  :: MODEL_ALPHA
    real(kind = 8)  :: MODEL_X_LOW
    real(kind = 8)  :: MODEL_X_HIGH
    integer         :: MODEL_X_STEPS
    integer         :: MODEL_FUN_TYPE

    integer         :: n_args, iarg
    character(len=16), dimension(7) :: args

    ! Find out how many arguments have been provided. If it's not 6, abort
    n_args = command_argument_count()
    if (n_args /= 7) then
      if (n_args > 0) write(*,'(a,I3)') , "ERROR: Expected 7 arguments, but got", n_args
      write(*,*) , "Usage:"
      write(*,*) , "  ./backflow <n> <p_cutoff> <alpha> <x_low> <x_high> <x_steps> <fun_type>"
      write(*,*) , ""
      write(*,*) , "Example:"
      write(*,*) , "  ./backflow 1000 200.0 -4.0 -1.0 1.0 100 1"
      stop
    end if

    ! If we got the right number of arguments, put them into an array of characters
    do iarg=1, n_args
      call get_command_argument(iarg, args(iarg))
    end do

    ! Read from the array of arguments into the appropriate variables
    read(args(1), *) MODEL_N
    read(args(2), *) MODEL_P_CUTOFF
    read(args(3), *) MODEL_ALPHA
    read(args(4), *) MODEL_X_LOW
    read(args(5), *) MODEL_X_HIGH
    read(args(6), *) MODEL_X_STEPS
    read(args(7), *) MODEL_FUN_TYPE

    ! Write some information about the parameters being used this run
    write(*,'(a)') "Running with parameters:"
    write(*,'(a,I7)')    "  MODEL_N        = ", MODEL_N
    write(*,'(a,F12.4)') "  MODEL_P_CUTOFF = ", MODEL_P_CUTOFF
    write(*,'(a,F12.4)') "  MODEL_ALPHA    = ", MODEL_ALPHA
    write(*,'(a,F12.4)') "  MODEL_X_LOW    = ", MODEL_X_LOW
    write(*,'(a,F12.4)') "  MODEL_X_HIGH   = ", MODEL_X_HIGH
    write(*,'(a,I7)')    "  MODEL_X_STEPS  = ", MODEL_X_STEPS
    write(*,'(a,I7)')    "  FUN_TYPE       = ", MODEL_FUN_TYPE
    write(*,*) , ""

    allocate(points(MODEL_X_STEPS + 1))
    allocate(realMatrix(MODEL_N, MODEL_N))
    allocate(imagMatrix(MODEL_N, MODEL_N))

    allocate(w(MODEL_N))
    allocate(zi(MODEL_N, MODEL_N))
    allocate(zr(MODEL_N, MODEL_N))

    interval = (MODEL_X_HIGH - (MODEL_X_LOW)) / MODEL_X_STEPS

    call JumpModel(MODEL_N, MODEL_P_CUTOFF, MODEL_ALPHA, MODEL_FUN_TYPE)

    write (*,fmt = "(A)") "POINTS = ["

    do test = 0, MODEL_X_STEPS
        x_test = MODEL_X_LOW + test*interval

        write(*, fmt = '(F10.3, A)', advance = "no") x_test, ", "

        call GaussianCurrentKernel(0.1D0, x_test)

        do row = 1, MODEL_N
            do col = 1, row
                tempVar = BasisMatrixElement(row - 1, col - 1)

                realMatrix(row, col) = DREAL(tempVar)
                imagMatrix(row, col) = DIMAG(tempVar)
            end do

!            print*, "Row = ", row
        end do

        call ch(MODEL_N, realMatrix, imagMatrix, w, .FALSE., zr,zi,ierr)

        points(test + 1) = w(1)

        write(*, fmt = '(F18.12, A)') w(1), ";"

    end do

    write (*,*) "];"

    write (*,*)

    write (*,*) "plot(POINTS(:,1), POINTS(:,2));"

end program
