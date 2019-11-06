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

    character(8)  :: date
    character(10) :: time

    character(255) :: datetime
    character(255) :: logfile
    character(255) :: datafile
    character(255) :: gpltfile
    character(255) :: plotfile
    character(255) :: mlabfile
    character(255) :: cwd
    character(255) :: n_char
    character(255) :: p_char
    character(255) :: a_char

    write(*,'(a)') "********************"
    write(*,'(a)') "JUMP MODEL TEST"
    write(*,'(a)') "********************"
    write(*,*) ""

    ! Find out how many arguments have been provided. If it's not 7, warning
    n_args = command_argument_count()
    if (n_args /= 7) then
      write(*,'(a,I3)') "WARNING: Expected 7 arguments, but got", n_args
      write(*,*) ""
      write(*,*) "Usage:"
      write(*,*) "  ./backflow <n> <p_cutoff> <alpha> <x_low> <x_high> <x_steps> <fun_type>"
      write(*,*) ""
      write(*,*) "Example:"
      write(*,*) "  ./backflow 1000 200.0 0.1 -1.0 1.0 100 1"
      write(*,*) ""
      write(*,*) "Using default parameters"
      write(*,*) ""

      MODEL_N          =  150
      MODEL_P_CUTOFF   =   50.0D0
      MODEL_ALPHA      =    0.1D0
      MODEL_X_LOW      = -  0.3D0
      MODEL_X_HIGH     =    0.3D0
      MODEL_X_STEPS    =    5
      MODEL_FUN_TYPE   =    1
    else
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
    end if

    write(n_char, '(i0)') MODEL_N
    write(p_char, '(i0)') int(MODEL_P_CUTOFF)
    write(a_char, '(i0)') int(1000*MODEL_ALPHA)

    call date_and_time(date = date, time = time)
    call getcwd(cwd)

    datetime = date(1:4)//"_"//date(5:6)//"_"//date(7:8)
    datetime = trim(datetime)//"__"//time(1:2)//"_"//time(3:4)//"_"//time(5:6)
    datetime = trim(datetime)//"__N_"//trim(n_char)//"__P_"//trim(p_char)//"__A_"//trim(a_char)

    logfile  = trim(cwd)//"/out/log__"//trim(datetime)//".txt"
    datafile = trim(cwd)//"/out/data__"//trim(datetime)//".csv"
    gpltfile = trim(cwd)//"/out/gplt__"//trim(datetime)//".plt"
    plotfile = trim(cwd)//"/out/plot__"//trim(datetime)//".png"
    mlabfile = trim(cwd)//"/out/mlab__"//trim(datetime)//".m"

    open(1, file = logfile, status = 'new')
    open(2, file = datafile, status = 'new')
    open(4, file = mlabfile, status = 'new')

    ! Write some information about the parameters being used this run
    write(*,'(a)') "********************"
    write(*,'(a)') ""
    write(*,'(a)') "Running with parameters:"
    write(*,'(a,I7)')    "  MODEL_N        = ", MODEL_N
    write(*,'(a,F12.4)') "  MODEL_P_CUTOFF = ", MODEL_P_CUTOFF
    write(*,'(a,F12.4)') "  MODEL_ALPHA    = ", MODEL_ALPHA
    write(*,'(a,F12.4)') "  MODEL_X_LOW    = ", MODEL_X_LOW
    write(*,'(a,F12.4)') "  MODEL_X_HIGH   = ", MODEL_X_HIGH
    write(*,'(a,I7)')    "  MODEL_X_STEPS  = ", MODEL_X_STEPS
    write(*,'(a,I7)')    "  FUN_TYPE       = ", MODEL_FUN_TYPE
    write(*,'(a)') ""
    write(*,'(a)') "********************"
    write(*,*) ""

    ! Write some information about the parameters being used this run
    write(1,'(a)') "********************"
    write(1,'(a)') ""
    write(1,'(a)') "Running with parameters:"
    write(1,'(a,I7)')    "  MODEL_N        = ", MODEL_N
    write(1,'(a,F12.4)') "  MODEL_P_CUTOFF = ", MODEL_P_CUTOFF
    write(1,'(a,F12.4)') "  MODEL_ALPHA    = ", MODEL_ALPHA
    write(1,'(a,F12.4)') "  MODEL_X_LOW    = ", MODEL_X_LOW
    write(1,'(a,F12.4)') "  MODEL_X_HIGH   = ", MODEL_X_HIGH
    write(1,'(a,I7)')    "  MODEL_X_STEPS  = ", MODEL_X_STEPS
    write(1,'(a,I7)')    "  FUN_TYPE       = ", MODEL_FUN_TYPE
    write(1,'(a)') ""
    write(1,'(a)') "********************"
    write(1,*) ""

    allocate(points(MODEL_X_STEPS + 1))
    allocate(realMatrix(MODEL_N, MODEL_N))
    allocate(imagMatrix(MODEL_N, MODEL_N))

    allocate(w(MODEL_N))
    allocate(zi(MODEL_N, MODEL_N))
    allocate(zr(MODEL_N, MODEL_N))

    interval = (MODEL_X_HIGH - (MODEL_X_LOW)) / MODEL_X_STEPS

    call JumpModel(MODEL_N, MODEL_P_CUTOFF, MODEL_ALPHA, MODEL_FUN_TYPE)
    call GeneratePlotFile(gpltfile, datafile, plotfile, MODEL_X_LOW, MODEL_X_HIGH)

    write(*,'(a,F12.4)') "backflow1 (bf1): MODEL_ALPHA = ", - MODEL_ALPHA
    write(*,'(a,F12.4)') "backflow2 (bf2): MODEL_ALPHA = ", MODEL_ALPHA
    write(*,*) ""

    write(1,'(a,F12.4)') "backflow1 (bf1): MODEL_ALPHA = ", - MODEL_ALPHA
    write(1,'(a,F12.4)') "backflow2 (bf2): MODEL_ALPHA = ", MODEL_ALPHA
    write(1,*) ""

    write(4,'(a)') "POINTS = ["


    do test = 0, MODEL_X_STEPS
        x_test = MODEL_X_LOW + test*interval

        write(*, fmt = '(A, F6.2, A)', advance = "no") "x = ", x_test, ", bf1 = "
        write(1, fmt = '(A, F6.2, A)', advance = "no") "x = ", x_test, ", bf1 = "
        write(2, fmt = '(F6.2, A)', advance = "no") x_test, "; "
        write(4, fmt = '(F6.2, A)', advance = "no") x_test, ", "

        call SetAlpha(- MODEL_ALPHA)

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

        write(*, fmt = '(F15.12, A)', advance = "no") w(1), ", bf2 = "
        write(1, fmt = '(F15.12, A)', advance = "no") w(1), ", bf2 = "
        write(2, fmt = '(F15.12, A)', advance = "no") w(1), "; "
        write(4, fmt = '(F15.12, A)', advance = "no") w(1), ", "

        call SetAlpha(MODEL_ALPHA)

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

        write(*, fmt = '(F15.12)') w(1)
        write(1, fmt = '(F15.12)') w(1)
        write(2, fmt = '(F15.12)') w(1)
        write(4, fmt = '(F15.12 a)') w(1), ";"

        !call system("gnuplot -p "//trim(gpltfile))

    end do

    write(4, fmt = '(a)') "];"
    write(4, fmt = '(a)') ""
    write(4, fmt = '(a)') "plot(POINTS(:,1), POINTS(:,2:3));"

    close(1)
    close(2)
    close(4)

end program
