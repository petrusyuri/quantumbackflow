module jump_integration

    use quadpack_double

    implicit none

    real (kind = 8) N, P_CUTOFF, ALPHA
    real (kind = 8) kernel_p, kernel_q

    real (kind = 8), parameter :: PI = 4.D0 * DATAN(1.D0)
    real (kind = 8) x0, sigma

    integer (kind = 4) FUN_TYPE

contains

    subroutine JumpModel(N_in, P_CUTOFF_in, ALPHA_in, FUN_TYPE_in)
        integer (kind = 4) N_in, FUN_TYPE_in
        real (kind = 8) P_CUTOFF_in, ALPHA_in

        N = REAL(N_in, kind = 8)
        P_CUTOFF = P_CUTOFF_in
        ALPHA    = ALPHA_in
        FUN_TYPE = FUN_TYPE_in
    end subroutine

    subroutine SetAlpha(ALPHA_in)
        real (kind = 8) ALPHA_in
        ALPHA = ALPHA_in
    end subroutine

    subroutine GaussianCurrentKernel(sigma_in, x0_in)
        real (kind = 8) sigma_in, x0_in

        sigma = sigma_in
        x0 = x0_in

        return
    end subroutine

    function BasisMatrixElement(i, j)
        integer i, j
        real (kind = 8) step, p, q
        complex (kind = 8) BasisMatrixElement

        step = P_CUTOFF / N

        p = step * (i + 0.5)
        q = step * (j + 0.5)

        BasisMatrixElement = step * kernelValue(p, q)
        return
    end function

    ! -------------------------------
    ! JUMP POTENTIAL MODEL
    ! -------------------------------

    function transmit(k)
        real (kind = 8) k
        complex (kind = 8) transmit

        transmit = DCMPLX(k, ALPHA) / DCMPLX(k, -ALPHA)

        return
    end function

    function reflect(k)
        real (kind = 8) k
        complex (kind = 8) reflect

        reflect = DCMPLX(0*k, 0*k)
        return
    end function

    function phi(k, x)
        real (kind = 8) k, x
        complex (kind = 8) phi, incomingWave, reflectWave, transmitWave

        if (x < 0) then
            incomingWave = EXP(DCMPLX(0, k*x))
            reflectWave  = CONJG(incomingWave)
            phi = reflectWave*reflect(k) + incomingWave
        else
            transmitWave = EXP(DCMPLX(0, k*x))
            phi = transmitWave*transmit(k)
        end if
        return
    end function

    function phiDeriv(k, x)
        real (kind = 8) k, x
        complex (kind = 8) phiDeriv, incomingWave, reflectWave, transmitWave

        if (x < 0) then
            incomingWave = EXP(DCMPLX(0, k*x)) * DCMPLX(0, k)
            reflectWave  = CONJG(incomingWave)
            phiDeriv = reflectWave*reflect(k) + incomingWave
        else
            transmitWave = EXP(DCMPLX(0, k*x)) * DCMPLX(0, k)
            phiDeriv = transmitWave*transmit(k)
        end if
        return
    end function

    ! -------------------------------
    ! SCATTERING MODEL
    ! -------------------------------

    function currentKernel(x, k1, k2)
        real (kind = 8) x, k1, k2
        complex (kind = 8) currentKernel, psi1, psi2, psiD1, psiD2, psiDiff

        psi1  = CONJG(phi(k1, x))

        psi2  = phi(k2, x)
        psiD1 = CONJG(phiDeriv(k1, x))
        psiD2 = phiDeriv(k2, x)
        psiDiff = psiD1*psi2 - psi1*psiD2

        currentKernel = DCMPLX(0, 0.25D0/PI) * psiDiff
        return
    end function

    ! -------------------------------
    ! GAUSSIAN CURRENT KERNEL
    ! -------------------------------

    function stepFunction(x)
        real (kind = 8) stepFunction, x

        if (FUN_TYPE .EQ. 1) then

            stepFunction = EXP(- (((x - x0)/sigma)**2)/2) / (sigma * SQRT(2*PI))

        else if (FUN_TYPE .EQ. 2) then

            stepFunction = (2*(sigma)**3)*(1/PI)/((x - x0)**2 + sigma**2)**2

        else if (FUN_TYPE .EQ. 3) then

            if((x < x0 - 0.5) .OR. (x > x0 + 0.5)) then
                stepFunction = 0
            else
                stepFunction = 1
            end if
        end if

        return
    end function

    ! -------------------------------
    ! SPATIAL CURRENT KERNEL
    ! -------------------------------

    function kernelValue(p, q)
        real (kind = 8) x_a, x_b
        real (kind = 8) p, q
        complex (kind = 8) kernelValue

        if (FUN_TYPE .EQ. 1) then
            x_a = x0 - 8*sigma
            x_b = x0 + 8*sigma
        else if (FUN_TYPE .EQ. 2) then
            x_a = x0 - 8*sigma
            x_b = x0 + 8*sigma
        else if (FUN_TYPE .EQ. 3) then
            x_a = x0 - 0.5
            x_b = x0 + 0.5
        end if

        kernelValue = Integrate(x_a, x_b, p, q)
        return
    end function

    function IntegrandReal(x)
        real (kind = 8) x, k1, k2, IntegrandReal

        k1 = kernel_p
        k2 = kernel_q

        IntegrandReal = stepFunction(x) * DREAL(currentKernel(x, k1, k2))
        return
    end function

    function IntegrandImag(x)
        real (kind = 8) x, k1, k2, IntegrandImag

        k1 = kernel_p
        k2 = kernel_q

        IntegrandImag = stepFunction(x) * DIMAG(currentKernel(x, k1, k2))
        return
    end function

    function Integrate(x_a, x_b, p, q)
        real ( kind = 8 ) x_a, x_b, p, q, realPart, imagPart

        integer ( kind = 4 ), parameter :: limit = 5000
        integer ( kind = 4 ), parameter :: lenw = 4 * limit
        real ( kind = 8 ) abserr
        real ( kind = 8 ), parameter :: epsabs = 1.0D-10
        real ( kind = 8 ), parameter :: epsrel = 1.0D-10
        integer ( kind = 4 ) ier
        integer ( kind = 4 ), parameter :: key = 6
        integer ( kind = 4 ) last
        integer ( kind = 4 ) neval

        integer ( kind = 4 ), parameter :: npts2 = 3
        real (kind = 8), dimension(npts2) :: points, pts
        integer (kind = 4), dimension(npts2) :: ndin

        real (kind = 8), dimension(limit) :: alist, blist, rlist, elist
        integer (kind = 4), dimension(limit) :: level, iord

        complex (kind = 8) Integrate

        kernel_p = p
        kernel_q = q

        if(x_a * x_b <= 0) then
            points(1) = 0

            call dqagpe (IntegrandReal, x_a, x_b, npts2, points, epsabs, epsrel, limit, realPart, abserr, neval, ier, &
            alist,blist,rlist,elist,pts,iord,level,ndin,last)

            call dqagpe (IntegrandImag, x_a, x_b, npts2, points, epsabs, epsrel, limit, imagPart, abserr, neval, ier, &
            alist,blist,rlist,elist,pts,iord,level,ndin,last)
        else
            call dqagpe (IntegrandReal, x_a, x_b, 2, points, epsabs, epsrel, limit, realPart, abserr, neval, ier, &
            alist,blist,rlist,elist,pts,iord,level,ndin,last)

            call dqagpe (IntegrandImag, x_a, x_b, 2, points, epsabs, epsrel, limit, imagPart, abserr, neval, ier, &
            alist,blist,rlist,elist,pts,iord,level,ndin,last)
        end if

        Integrate = DCMPLX(realPart, imagPart)
        return
    end function

    ! -------------------------------
    ! GNUPLOT ROUTINE
    ! -------------------------------

    subroutine GeneratePlotFile(gpltfile, datafile, plotfile, X_LOW, X_HIGH)
        character(255) :: gpltfile, datafile, plotfile
        real(kind = 8)  :: X_LOW, X_HIGH

        open(3, file = gpltfile, status = 'new')

        write(3,'(a)') "set encoding utf8"
        write(3,'(a)') "set terminal pdf size 960, 720 enhanced"
        write(3,'(a)') "set output '"//trim(plotfile)//"'"
        write(3,'(a)') "set datafile separator ';'"
        write(3,'(a)') ""
        write(3,'(a)', advance = "no") "set title 'Backflow in jump defect: "
        write(3,'(a, i4, a, i4, a, F8.2)') "N = ", int(N), ", P_{cutoff} = ", int(P_CUTOFF), ", |α| = ", abs(ALPHA)
        write(3,'(a)') "set style line 1 linecolor rgb 'red' linewidth 2"
        write(3,'(a)') "set style line 2 linecolor rgb 'blue' linewidth 2"
        write(3,'(a)') "set style line 3 linecolor rgb 'forest-green' linewidth 2"
        write(3,'(a)') ""
        write(3,'(a)') "set label ''"
        write(3,'(a)') "set xlabel 'x_{0}'"
        write(3,'(a)') "set ylabel 'β_V(f)'"
        write(3,'(a)') "set key center right"
        write(3,'(a)') ""
        write(3,'(a, f6.2, a, f6.2, a)') "set xrange [", X_LOW, ":", X_HIGH, "]"
        write(3,'(a)') "set yrange [*:0]"
        write(3,'(a)') ""
        write(3,'(a)') "plot '"//trim(datafile)//"' using 1:2 with lines ls 1 title 'α < 0', \"
        write(3,'(a)') "     '"//trim(datafile)//"' using 1:3 with lines ls 2 title 'α > 0'"

        close(3)
    end subroutine
end module
