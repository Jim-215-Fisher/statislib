







module statislib_special
    use stdlib_kinds
    use stdlib_error, only : error_stop

    implicit none
    private
    real(qp), parameter :: DD(0:10) = [2.48574089138753565546e-5_qp,            &
                                     1.05142378581721974210_qp,                 &
                                     -3.45687097222016235469_qp,                &
                                     4.51227709466894823700_qp,                 &
                                     -2.98285225323576655721_qp,                &
                                     1.05639711577126713077_qp,                 &
                                     -1.95428773191645869583e-1_qp,             &
                                     1.70970543404441224307e-2_qp,              &
                                     -5.71926117404305781283e-4_qp,             &
                                     4.63399473359905636708e-6_qp,              &
                                     -2.71994908488607703910e-9_qp]
    ! Coefficients of Lanczos approximation for 16 digits float precision

    real(qp), parameter :: DQ(0:21)= [2.0240434640140357514731512432760e-10_qp, &
                                      1.5333183020199267370932516012553_qp,     &
                                     -1.1640274608858812982567477805332e1_qp,   &
                                      4.0053698000222503376927701573076e1_qp,   &
                                     -8.2667863469173479039227422723581e1_qp,   &
                                      1.1414465885256804336106748692495e2_qp,   &
                                     -1.1135645608449754488425056563075e2_qp,   &
                                      7.9037451549298877731413453151252e1_qp,   &
                                     -4.1415428804507353801947558814560e1_qp,   &
                                      1.6094742170165161102085734210327e1_qp,   &
                                     -4.6223809979028638614212851576524_qp,     &
                                      9.7030884294357827423006360746167e-1_qp,  &
                                     -1.4607332380456449418243363858893e-1_qp,  &
                                      1.5330325530769204955496334450658e-2_qp,  &
                                     -1.0773862404547660506042948153734e-3_qp,  &
                                      4.7911128916072940196391032755132e-5_qp,  &
                                     -1.2437781042887028450811158692678e-6_qp,  &
                                      1.6751019107496606112103160490729e-8_qp,  &
                                     -9.7674656970897286097939311684868e-11_qp, &
                                      1.8326577220560509759575892664132e-13_qp, &
                                     -6.4508377189118502115673823719605e-17_qp, &
                                      1.3382662604773700632782310392171e-21_qp]
    ! Coefficients of Lanczos approximation for 32 digits float precision

    real(qp), parameter :: RD = 10.900511_qp, RQ = 22.618910_qp, HALF = 0.5_qp, &
              sqepi = log(2.0_qp * sqrt(exp(1.0_qp) / acos(-1.0_qp)))
    real(sp), parameter :: tol_sp = epsilon(1.0_sp)
    real(dp), parameter :: tol_dp = epsilon(1.0_dp)
    real(qp), parameter :: tol_qp = epsilon(1.0_qp)
    real(dp), parameter :: dmd = 1.0e-300_dp
    real(qp), parameter :: dmq = 1.0e-4900_qp


    public :: log_gamma, log_factorial
    public :: lower_incomplete_gamma, log_lower_incomplete_gamma
    public :: upper_incomplete_gamma, log_upper_incomplete_gamma
    public :: regularized_gamma_p, regularized_gamma_q
    public :: beta, log_beta, incomplete_beta, regularized_beta


    interface log_gamma
    ! Logrithm of gamma function with real variable
    !
        module procedure l_gamma_rsp
        module procedure l_gamma_rdp
        module procedure l_gamma_rqp
        module procedure l_gamma_csp
        module procedure l_gamma_cdp
        module procedure l_gamma_cqp
    end interface log_gamma



    interface log_factorial
    ! Logrithm of factorial n!, integer variable
    !
        module procedure l_factorial_1_iint8      !1 dummy
        module procedure l_factorial_1_iint16      !1 dummy
        module procedure l_factorial_1_iint32      !1 dummy
        module procedure l_factorial_1_iint64      !1 dummy

        module procedure l_factorial_iint8sp !2 dummy
        module procedure l_factorial_iint8dp !2 dummy
        module procedure l_factorial_iint8qp !2 dummy
        module procedure l_factorial_iint16sp !2 dummy
        module procedure l_factorial_iint16dp !2 dummy
        module procedure l_factorial_iint16qp !2 dummy
        module procedure l_factorial_iint32sp !2 dummy
        module procedure l_factorial_iint32dp !2 dummy
        module procedure l_factorial_iint32qp !2 dummy
        module procedure l_factorial_iint64sp !2 dummy
        module procedure l_factorial_iint64dp !2 dummy
        module procedure l_factorial_iint64qp !2 dummy
    end interface log_factorial



    interface lower_incomplete_gamma
    ! Lower incomplete gamma function
    !
        module procedure ingamma_low_iint8sp
        module procedure ingamma_low_iint8dp
        module procedure ingamma_low_iint8qp
        module procedure ingamma_low_iint16sp
        module procedure ingamma_low_iint16dp
        module procedure ingamma_low_iint16qp
        module procedure ingamma_low_iint32sp
        module procedure ingamma_low_iint32dp
        module procedure ingamma_low_iint32qp
        module procedure ingamma_low_iint64sp
        module procedure ingamma_low_iint64dp
        module procedure ingamma_low_iint64qp
        module procedure ingamma_low_rspsp
        module procedure ingamma_low_rdpdp
        module procedure ingamma_low_rqpqp
    end interface lower_incomplete_gamma



    interface log_lower_incomplete_gamma
    ! Logrithm of lower incomplete gamma function
    !
        module procedure l_ingamma_low_iint8sp
        module procedure l_ingamma_low_iint8dp
        module procedure l_ingamma_low_iint8qp
        module procedure l_ingamma_low_iint16sp
        module procedure l_ingamma_low_iint16dp
        module procedure l_ingamma_low_iint16qp
        module procedure l_ingamma_low_iint32sp
        module procedure l_ingamma_low_iint32dp
        module procedure l_ingamma_low_iint32qp
        module procedure l_ingamma_low_iint64sp
        module procedure l_ingamma_low_iint64dp
        module procedure l_ingamma_low_iint64qp
        module procedure l_ingamma_low_rspsp
        module procedure l_ingamma_low_rdpdp
        module procedure l_ingamma_low_rqpqp
    end interface log_lower_incomplete_gamma



    interface upper_incomplete_gamma
    ! Upper incomplete gamma function
    !
        module procedure ingamma_up_iint8sp
        module procedure ingamma_up_iint8dp
        module procedure ingamma_up_iint8qp
        module procedure ingamma_up_iint16sp
        module procedure ingamma_up_iint16dp
        module procedure ingamma_up_iint16qp
        module procedure ingamma_up_iint32sp
        module procedure ingamma_up_iint32dp
        module procedure ingamma_up_iint32qp
        module procedure ingamma_up_iint64sp
        module procedure ingamma_up_iint64dp
        module procedure ingamma_up_iint64qp
        module procedure ingamma_up_rspsp
        module procedure ingamma_up_rdpdp
        module procedure ingamma_up_rqpqp
    end interface upper_incomplete_gamma



    interface log_upper_incomplete_gamma
    ! Logrithm of upper incomplete gamma function
        module procedure l_ingamma_up_iint8sp
        module procedure l_ingamma_up_iint8dp
        module procedure l_ingamma_up_iint8qp
        module procedure l_ingamma_up_iint16sp
        module procedure l_ingamma_up_iint16dp
        module procedure l_ingamma_up_iint16qp
        module procedure l_ingamma_up_iint32sp
        module procedure l_ingamma_up_iint32dp
        module procedure l_ingamma_up_iint32qp
        module procedure l_ingamma_up_iint64sp
        module procedure l_ingamma_up_iint64dp
        module procedure l_ingamma_up_iint64qp
        module procedure l_ingamma_up_rspsp
        module procedure l_ingamma_up_rdpdp
        module procedure l_ingamma_up_rqpqp
    end interface log_upper_incomplete_gamma



    interface regularized_gamma_p
    ! Regularized (normalized) lower incomplete gamma function, P
    !
        module procedure regamma_p_iint8sp
        module procedure regamma_p_iint8dp
        module procedure regamma_p_iint8qp
        module procedure regamma_p_iint16sp
        module procedure regamma_p_iint16dp
        module procedure regamma_p_iint16qp
        module procedure regamma_p_iint32sp
        module procedure regamma_p_iint32dp
        module procedure regamma_p_iint32qp
        module procedure regamma_p_iint64sp
        module procedure regamma_p_iint64dp
        module procedure regamma_p_iint64qp
        module procedure regamma_p_rspsp
        module procedure regamma_p_rdpdp
        module procedure regamma_p_rqpqp
    end interface regularized_gamma_p



    interface regularized_gamma_q
    ! Regularized (normalized) upper incomplete gamma function, Q
    !
        module procedure regamma_q_iint8sp
        module procedure regamma_q_iint8dp
        module procedure regamma_q_iint8qp
        module procedure regamma_q_iint16sp
        module procedure regamma_q_iint16dp
        module procedure regamma_q_iint16qp
        module procedure regamma_q_iint32sp
        module procedure regamma_q_iint32dp
        module procedure regamma_q_iint32qp
        module procedure regamma_q_iint64sp
        module procedure regamma_q_iint64dp
        module procedure regamma_q_iint64qp
        module procedure regamma_q_rspsp
        module procedure regamma_q_rdpdp
        module procedure regamma_q_rqpqp
    end interface regularized_gamma_q



    interface gpx
    ! Evaluation of incomplete gamma function
    !
        module procedure gpx_rsp
        module procedure gpx_rdp
        module procedure gpx_rqp

        module procedure gpx_iint8sp
        module procedure gpx_iint8dp
        module procedure gpx_iint8qp
        module procedure gpx_iint16sp
        module procedure gpx_iint16dp
        module procedure gpx_iint16qp
        module procedure gpx_iint32sp
        module procedure gpx_iint32dp
        module procedure gpx_iint32qp
        module procedure gpx_iint64sp
        module procedure gpx_iint64dp
        module procedure gpx_iint64qp
    end interface gpx



    interface beta
    ! Beta function
    !
        module procedure beta_rsp
        module procedure beta_rdp
        module procedure beta_rqp
        module procedure beta_csp
        module procedure beta_cdp
        module procedure beta_cqp
    end interface beta



    interface log_beta
    ! Logrithm of beta function
    !
        module procedure l_beta_rsp
        module procedure l_beta_rdp
        module procedure l_beta_rqp
        module procedure l_beta_csp
        module procedure l_beta_cdp
        module procedure l_beta_cqp
    end interface log_beta



    interface incomplete_beta
    ! Incomplete beta function
    !
        module procedure inbeta_rsp
        module procedure inbeta_rdp
        module procedure inbeta_rqp
    end interface incomplete_beta




    interface regularized_beta
    ! Regularized incomplete beta function
    !
        module procedure regbeta_rsp
        module procedure regbeta_rdp
        module procedure regbeta_rqp
    end interface regularized_beta





contains

    impure elemental function l_gamma_rsp(x) result (res)
    !
    ! Log gamma function for any positive real number i,e, {R+}
    ! Algorithm is based on Glendon Pugh, "An Analysis of The Lanczos Gamma
    ! Approximation", The University of British Columbia, 2004
    !
    ! Fortran 90 program by Jim-215-Fisher
    !
        real(sp), intent(in) :: x
        real(sp) :: res
        real(qp) :: q, sum
        integer :: i

        if(x <= 0._sp) call error_stop("Error(l_gamma): Logrithm of Gamma "&
             //"function augument must be greater than zero.")
        if(x == 1.0_sp .or. x == 2.0_sp) then
            res = 0.0_sp
        else
            q = x - HALF
            ! single precision use 16 digits coefficients
            sum = DD(0)
            do i = 1, 10
                sum = sum + DD(i) / (x - 1.0_qp + i)
            end do
            res = real(sqepi + log(sum) - q + q * log(q + RD), kind=sp)
        endif
    end function l_gamma_rsp

    impure elemental function l_gamma_rdp(x) result (res)
    !
    ! Log gamma function for any positive real number i,e, {R+}
    ! Algorithm is based on Glendon Pugh, "An Analysis of The Lanczos Gamma
    ! Approximation", The University of British Columbia, 2004
    !
    ! Fortran 90 program by Jim-215-Fisher
    !
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: q, sum
        integer :: i

        if(x <= 0._dp) call error_stop("Error(l_gamma): Logrithm of Gamma "&
             //"function augument must be greater than zero.")
        if(x == 1.0_dp .or. x == 2.0_dp) then
            res = 0.0_dp
        else
            q = x - HALF
            ! double/quadruple precision use 32 digits coef.
            sum = DQ(0)
            do i=1, 21
                sum = sum + DQ(i) / (x - 1.0_qp + i)
            end do
            res = real(sqepi + log(sum) - q + q * log(q + RQ), kind=dp)
        endif
    end function l_gamma_rdp

    impure elemental function l_gamma_rqp(x) result (res)
    !
    ! Log gamma function for any positive real number i,e, {R+}
    ! Algorithm is based on Glendon Pugh, "An Analysis of The Lanczos Gamma
    ! Approximation", The University of British Columbia, 2004
    !
    ! Fortran 90 program by Jim-215-Fisher
    !
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: q, sum
        integer :: i

        if(x <= 0._qp) call error_stop("Error(l_gamma): Logrithm of Gamma "&
             //"function augument must be greater than zero.")
        if(x == 1.0_qp .or. x == 2.0_qp) then
            res = 0.0_qp
        else
            q = x - HALF
            ! double/quadruple precision use 32 digits coef.
            sum = DQ(0)
            do i=1, 21
                sum = sum + DQ(i) / (x - 1.0_qp + i)
            end do
            res = real(sqepi + log(sum) - q + q * log(q + RQ), kind=qp)
        endif
    end function l_gamma_rqp




    impure elemental function l_gamma_csp(x) result (res)
    !
    ! Log gamma function for any complex number with Re(x) > 0
    !
        complex(sp), intent(in) :: x
        complex(sp) :: res
        real(qp) :: q0, q, p, lgr, sumr, sumi, tr, ti, theta1, theta2, t0, t1
        integer :: i

        if(x % re <= 0._sp) call error_stop("Error(l_gamma): Logrithm of " &
            //"Gamma function complex augument must have a positive real part")
         
        q0 = x % re - HALF; p = x % im
        q = q0 + RD; theta1 = atan2(p, q); lgr = HALF * log(q * q + p * p)
        tr = sqepi + q0 * (lgr - 1.0_qp) - theta1 * p
        ti = p * (lgr - 1.0_qp) + theta1 * q0
        sumr = DD(0); sumi = 0.0_qp
        do i = 1, 10
            t0 = q0 - HALF + i
            t1 = t0 * t0 + p * p; t1 = 1.0_qp / t1
            sumr = sumr + DD(i) * t0 * t1
            sumi = sumi + DD(i) * p * t1
        end do
        sumi = - sumi
        theta2 = atan2(sumi, sumr)
        tr = tr + HALF * log(sumr * sumr + sumi * sumi)
        ti = ti + theta2
        res = cmplx(tr, ti, kind = sp)
    end function l_gamma_csp

    impure elemental function l_gamma_cdp(x) result (res)
    !
    ! Log gamma function for any complex number with Re(x) > 0
    !
        complex(dp), intent(in) :: x
        complex(dp) :: res
        real(qp) :: q0, q, p, lgr, sumr, sumi, tr, ti, theta1, theta2, t0, t1
        integer :: i

        if(x % re <= 0._dp) call error_stop("Error(l_gamma): Logrithm of " &
            //"Gamma function complex augument must have a positive real part")
         
        q0 = x % re - HALF; p = x % im
        q = q0 + RQ; theta1 = atan2(p, q); lgr = HALF * log(q * q + p * p)
        tr = sqepi + q0 * (lgr - 1.0_qp) - theta1 * p
        ti = p * (lgr - 1.0_qp) + theta1 * q0
        sumr = DQ(0); sumi = 0.0_qp
        do i = 1, 21
            t0 = q0 - HALF + i
            t1 = t0 * t0 + p * p; t1 = 1.0_qp / t1
            sumr = sumr + DQ(i) * t0 * t1
            sumi = sumi + DQ(i) * p * t1
        end do
        sumi = - sumi
        theta2 = atan2(sumi, sumr)
        tr = tr + HALF * log(sumr * sumr + sumi * sumi)
        ti = ti + theta2
        res = cmplx(tr, ti, kind = dp)
    end function l_gamma_cdp

    impure elemental function l_gamma_cqp(x) result (res)
    !
    ! Log gamma function for any complex number with Re(x) > 0
    !
        complex(qp), intent(in) :: x
        complex(qp) :: res
        real(qp) :: q0, q, p, lgr, sumr, sumi, tr, ti, theta1, theta2, t0, t1
        integer :: i

        if(x % re <= 0._qp) call error_stop("Error(l_gamma): Logrithm of " &
            //"Gamma function complex augument must have a positive real part")
         
        q0 = x % re - HALF; p = x % im
        q = q0 + RQ; theta1 = atan2(p, q); lgr = HALF * log(q * q + p * p)
        tr = sqepi + q0 * (lgr - 1.0_qp) - theta1 * p
        ti = p * (lgr - 1.0_qp) + theta1 * q0
        sumr = DQ(0); sumi = 0.0_qp
        do i = 1, 21
            t0 = q0 - HALF + i
            t1 = t0 * t0 + p * p; t1 = 1.0_qp / t1
            sumr = sumr + DQ(i) * t0 * t1
            sumi = sumi + DQ(i) * p * t1
        end do
        sumi = - sumi
        theta2 = atan2(sumi, sumr)
        tr = tr + HALF * log(sumr * sumr + sumi * sumi)
        ti = ti + theta2
        res = cmplx(tr, ti, kind = qp)
    end function l_gamma_cqp




    impure elemental function l_factorial_1_iint8(n) result(res)
    !
    ! Log(n!) with single precision result, n is integer
    !
        integer(int8), intent(in) :: n
        real :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
            //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0
        case (1)
            res = 0.0
        case (2:)
            res = log_gamma(real(n + 1, dp))
        end select
    end function l_factorial_1_iint8
     
    impure elemental function l_factorial_1_iint16(n) result(res)
    !
    ! Log(n!) with single precision result, n is integer
    !
        integer(int16), intent(in) :: n
        real :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
            //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0
        case (1)
            res = 0.0
        case (2:)
            res = log_gamma(real(n + 1, dp))
        end select
    end function l_factorial_1_iint16
     
    impure elemental function l_factorial_1_iint32(n) result(res)
    !
    ! Log(n!) with single precision result, n is integer
    !
        integer(int32), intent(in) :: n
        real :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
            //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0
        case (1)
            res = 0.0
        case (2:)
            res = log_gamma(real(n + 1, dp))
        end select
    end function l_factorial_1_iint32
     
    impure elemental function l_factorial_1_iint64(n) result(res)
    !
    ! Log(n!) with single precision result, n is integer
    !
        integer(int64), intent(in) :: n
        real :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
            //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0
        case (1)
            res = 0.0
        case (2:)
            res = log_gamma(real(n + 1, dp))
        end select
    end function l_factorial_1_iint64
     



    impure elemental function l_factorial_iint8sp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int8), intent(in) :: n
        real(sp), intent(in) :: x
        real(sp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_sp
        case (1)
            res = 0.0_sp
        case (2:)
            res = log_gamma(real(n + 1, kind = sp))
        end select
    end function l_factorial_iint8sp
     
    impure elemental function l_factorial_iint8dp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int8), intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_dp
        case (1)
            res = 0.0_dp
        case (2:)
            res = log_gamma(real(n + 1, kind = dp))
        end select
    end function l_factorial_iint8dp
     
    impure elemental function l_factorial_iint8qp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int8), intent(in) :: n
        real(qp), intent(in) :: x
        real(qp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_qp
        case (1)
            res = 0.0_qp
        case (2:)
            res = log_gamma(real(n + 1, kind = qp))
        end select
    end function l_factorial_iint8qp
     
    impure elemental function l_factorial_iint16sp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int16), intent(in) :: n
        real(sp), intent(in) :: x
        real(sp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_sp
        case (1)
            res = 0.0_sp
        case (2:)
            res = log_gamma(real(n + 1, kind = sp))
        end select
    end function l_factorial_iint16sp
     
    impure elemental function l_factorial_iint16dp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int16), intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_dp
        case (1)
            res = 0.0_dp
        case (2:)
            res = log_gamma(real(n + 1, kind = dp))
        end select
    end function l_factorial_iint16dp
     
    impure elemental function l_factorial_iint16qp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int16), intent(in) :: n
        real(qp), intent(in) :: x
        real(qp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_qp
        case (1)
            res = 0.0_qp
        case (2:)
            res = log_gamma(real(n + 1, kind = qp))
        end select
    end function l_factorial_iint16qp
     
    impure elemental function l_factorial_iint32sp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int32), intent(in) :: n
        real(sp), intent(in) :: x
        real(sp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_sp
        case (1)
            res = 0.0_sp
        case (2:)
            res = log_gamma(real(n + 1, kind = sp))
        end select
    end function l_factorial_iint32sp
     
    impure elemental function l_factorial_iint32dp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int32), intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_dp
        case (1)
            res = 0.0_dp
        case (2:)
            res = log_gamma(real(n + 1, kind = dp))
        end select
    end function l_factorial_iint32dp
     
    impure elemental function l_factorial_iint32qp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int32), intent(in) :: n
        real(qp), intent(in) :: x
        real(qp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_qp
        case (1)
            res = 0.0_qp
        case (2:)
            res = log_gamma(real(n + 1, kind = qp))
        end select
    end function l_factorial_iint32qp
     
    impure elemental function l_factorial_iint64sp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int64), intent(in) :: n
        real(sp), intent(in) :: x
        real(sp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_sp
        case (1)
            res = 0.0_sp
        case (2:)
            res = log_gamma(real(n + 1, kind = sp))
        end select
    end function l_factorial_iint64sp
     
    impure elemental function l_factorial_iint64dp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int64), intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_dp
        case (1)
            res = 0.0_dp
        case (2:)
            res = log_gamma(real(n + 1, kind = dp))
        end select
    end function l_factorial_iint64dp
     
    impure elemental function l_factorial_iint64qp(n,x) result(res)
    !
    ! Log(n!) with required prescision for result, n is integer, x is a real
    ! for specified kind
    !
        integer(int64), intent(in) :: n
        real(qp), intent(in) :: x
        real(qp) :: res

        if(n < 0) call error_stop("Error(l_factorial): Factorial function"      &
                  //" augument must be non-negative")
        select case(n)
        case (0)
            res = 0.0_qp
        case (1)
            res = 0.0_qp
        case (2:)
            res = log_gamma(real(n + 1, kind = qp))
        end select
    end function l_factorial_iint64qp
     



	    ! single precision use double variables
    impure elemental function gpx_rsp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with real auguments s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
    ! Fortran 90 program by Jim-215-Fisher
    !
        real(sp), intent(in) :: x, s
        real(dp) :: res
        real(sp) :: p_lim
        real(dp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmd
        if(x < -9.0_sp) then
            p_lim = 5.0_sp * (sqrt(abs(x)) - 1.0_sp)
        elseif(x >= -9.0_sp .and. x <= 0.0_sp) then
            p_lim = 0.0_sp
        else
            p_lim = x
        endif
        if(x < 0._sp) then
            call error_stop("Error(gpx): Incomplete gamma function with"        &
                //" negative x must come with integer of s")
        elseif(s >= p_lim) then
            a = 1.0_dp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        else
            a = 1.0_dp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_dp
                d = d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        endif
        res = g
    end function gpx_rsp

	    ! double/quadruple precision use quadruple variables
    impure elemental function gpx_rdp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with real auguments s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
    ! Fortran 90 program by Jim-215-Fisher
    !
        real(dp), intent(in) :: x, s
        real(qp) :: res
        real(dp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_dp) then
            p_lim = 5.0_dp * (sqrt(abs(x)) - 1.0_dp)
        elseif(x >= -9.0_dp .and. x <= 0.0_dp) then
            p_lim = 0.0_dp
        else
            p_lim = x
        endif
        if(x < 0._dp) then
            call error_stop("Error(gpx): Incomplete gamma function with"        &
                //" negative x must come with integer of s")
        elseif(s >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        endif
        res = g
    end function gpx_rdp

	    ! double/quadruple precision use quadruple variables
    impure elemental function gpx_rqp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with real auguments s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
    ! Fortran 90 program by Jim-215-Fisher
    !
        real(qp), intent(in) :: x, s
        real(qp) :: res
        real(qp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_qp) then
            p_lim = 5.0_qp * (sqrt(abs(x)) - 1.0_qp)
        elseif(x >= -9.0_qp .and. x <= 0.0_qp) then
            p_lim = 0.0_qp
        else
            p_lim = x
        endif
        if(x < 0._qp) then
            call error_stop("Error(gpx): Incomplete gamma function with"        &
                //" negative x must come with integer of s")
        elseif(s >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        endif
        res = g
    end function gpx_rqp




		  ! single precision use double variables
    impure elemental function gpx_iint8sp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(dp) :: res
        real(sp) :: p_lim
        real(dp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmd
        if(x < -9.0_sp) then
            p_lim = 5.0_sp * (sqrt(abs(x)) - 1.0_sp)
        elseif(x >= -9.0_sp .and. x <= 0.0_sp) then
            p_lim = 0.0_sp
        else
            p_lim = x
        endif
        if(real(s, sp) >= p_lim) then
            a = 1.0_dp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        elseif(x >= 0.0_sp) then
            a = 1.0_dp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_dp
                d = d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        else
            a = -x
            c = 1.0_dp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_dp) / (a * a)
                d = d - 2.0_dp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int8) > (s - 2_int8) / 2_int8 .or. y < b *     &
                    tol_dp) exit
            end do
            if(y >= b * tol_dp .and. mod(s, 2_int8) /= 0_int8)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, dp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint8sp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint8dp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(qp) :: res
        real(dp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_dp) then
            p_lim = 5.0_dp * (sqrt(abs(x)) - 1.0_dp)
        elseif(x >= -9.0_dp .and. x <= 0.0_dp) then
            p_lim = 0.0_dp
        else
            p_lim = x
        endif
        if(real(s, dp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_dp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int8) > (s - 2_int8) / 2_int8 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int8) /= 0_int8)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint8dp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint8qp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_qp) then
            p_lim = 5.0_qp * (sqrt(abs(x)) - 1.0_qp)
        elseif(x >= -9.0_qp .and. x <= 0.0_qp) then
            p_lim = 0.0_qp
        else
            p_lim = x
        endif
        if(real(s, qp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_qp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int8) > (s - 2_int8) / 2_int8 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int8) /= 0_int8)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint8qp

		  ! single precision use double variables
    impure elemental function gpx_iint16sp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(dp) :: res
        real(sp) :: p_lim
        real(dp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmd
        if(x < -9.0_sp) then
            p_lim = 5.0_sp * (sqrt(abs(x)) - 1.0_sp)
        elseif(x >= -9.0_sp .and. x <= 0.0_sp) then
            p_lim = 0.0_sp
        else
            p_lim = x
        endif
        if(real(s, sp) >= p_lim) then
            a = 1.0_dp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        elseif(x >= 0.0_sp) then
            a = 1.0_dp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_dp
                d = d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        else
            a = -x
            c = 1.0_dp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_dp) / (a * a)
                d = d - 2.0_dp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int16) > (s - 2_int16) / 2_int16 .or. y < b *     &
                    tol_dp) exit
            end do
            if(y >= b * tol_dp .and. mod(s, 2_int16) /= 0_int16)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, dp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint16sp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint16dp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(qp) :: res
        real(dp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_dp) then
            p_lim = 5.0_dp * (sqrt(abs(x)) - 1.0_dp)
        elseif(x >= -9.0_dp .and. x <= 0.0_dp) then
            p_lim = 0.0_dp
        else
            p_lim = x
        endif
        if(real(s, dp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_dp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int16) > (s - 2_int16) / 2_int16 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int16) /= 0_int16)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint16dp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint16qp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_qp) then
            p_lim = 5.0_qp * (sqrt(abs(x)) - 1.0_qp)
        elseif(x >= -9.0_qp .and. x <= 0.0_qp) then
            p_lim = 0.0_qp
        else
            p_lim = x
        endif
        if(real(s, qp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_qp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int16) > (s - 2_int16) / 2_int16 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int16) /= 0_int16)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint16qp

		  ! single precision use double variables
    impure elemental function gpx_iint32sp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(dp) :: res
        real(sp) :: p_lim
        real(dp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmd
        if(x < -9.0_sp) then
            p_lim = 5.0_sp * (sqrt(abs(x)) - 1.0_sp)
        elseif(x >= -9.0_sp .and. x <= 0.0_sp) then
            p_lim = 0.0_sp
        else
            p_lim = x
        endif
        if(real(s, sp) >= p_lim) then
            a = 1.0_dp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        elseif(x >= 0.0_sp) then
            a = 1.0_dp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_dp
                d = d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        else
            a = -x
            c = 1.0_dp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_dp) / (a * a)
                d = d - 2.0_dp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int32) > (s - 2_int32) / 2_int32 .or. y < b *     &
                    tol_dp) exit
            end do
            if(y >= b * tol_dp .and. mod(s, 2_int32) /= 0_int32)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, dp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint32sp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint32dp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(qp) :: res
        real(dp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_dp) then
            p_lim = 5.0_dp * (sqrt(abs(x)) - 1.0_dp)
        elseif(x >= -9.0_dp .and. x <= 0.0_dp) then
            p_lim = 0.0_dp
        else
            p_lim = x
        endif
        if(real(s, dp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_dp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int32) > (s - 2_int32) / 2_int32 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int32) /= 0_int32)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint32dp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint32qp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_qp) then
            p_lim = 5.0_qp * (sqrt(abs(x)) - 1.0_qp)
        elseif(x >= -9.0_qp .and. x <= 0.0_qp) then
            p_lim = 0.0_qp
        else
            p_lim = x
        endif
        if(real(s, qp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_qp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int32) > (s - 2_int32) / 2_int32 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int32) /= 0_int32)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint32qp

		  ! single precision use double variables
    impure elemental function gpx_iint64sp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(dp) :: res
        real(sp) :: p_lim
        real(dp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmd
        if(x < -9.0_sp) then
            p_lim = 5.0_sp * (sqrt(abs(x)) - 1.0_sp)
        elseif(x >= -9.0_sp .and. x <= 0.0_sp) then
            p_lim = 0.0_sp
        else
            p_lim = x
        endif
        if(real(s, sp) >= p_lim) then
            a = 1.0_dp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        elseif(x >= 0.0_sp) then
            a = 1.0_dp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_dp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_dp
                d = d * a + b
                if(d == 0.0_dp) d = dm
                c = b + a / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
        else
            a = -x
            c = 1.0_dp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_dp) / (a * a)
                d = d - 2.0_dp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int64) > (s - 2_int64) / 2_int64 .or. y < b *     &
                    tol_dp) exit
            end do
            if(y >= b * tol_dp .and. mod(s, 2_int64) /= 0_int64)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, dp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint64sp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint64dp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(qp) :: res
        real(dp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_dp) then
            p_lim = 5.0_dp * (sqrt(abs(x)) - 1.0_dp)
        elseif(x >= -9.0_dp .and. x <= 0.0_dp) then
            p_lim = 0.0_dp
        else
            p_lim = x
        endif
        if(real(s, dp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_dp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int64) > (s - 2_int64) / 2_int64 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int64) /= 0_int64)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint64dp

          ! double/quadruple precision use quadruple variables
    impure elemental function gpx_iint64qp(s, x) result(res)
    !
    ! Approximation of incomplete gamma G function with integer augument s.
    !
    ! Based on Rémy Abergel and Lionel Moisan "Algorithm 1006, Fast and
    ! Accurate Evaluation of a Generalized Incomplete Gamma Function", ACM
    ! Transactions on Mathematical Software, March 2020.
    !
        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: p_lim
        real(qp) :: a, b, g, c, d, y, dm
        integer :: n

        dm = dmq
        if(x < -9.0_qp) then
            p_lim = 5.0_qp * (sqrt(abs(x)) - 1.0_qp)
        elseif(x >= -9.0_qp .and. x <= 0.0_qp) then
            p_lim = 0.0_qp
        else
            p_lim = x
        endif
        if(real(s, qp) >= p_lim) then
            a = 1.0_qp
            b = s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                if(mod(n, 2) == 0) then
                    a = (1 - s - n / 2) * x
                else
                    a = (n / 2) * x
                end if
                b = s - 1 + n
                d =  d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        elseif(x >= 0.0_qp) then
            a = 1.0_qp
            b = x + 1 - s
            g = a / b
            c = a / dm
            d = 1.0_qp / b
            n = 2
            do
                a = -(n - 1) * (n - 1 - s)
                b = b + 2.0_qp
                d = d * a + b
                if(d == 0.0_qp) d = dm
                c = b + a / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
        else
            a = -x
            c = 1.0_qp / a
            d = s - 1
            b = c * (a - d)
            n = 1
            do
                c = d * (d - 1.0_qp) / (a * a)
                d = d - 2.0_qp
                y = c * ( a - d)
                b = b + y
                n = n + 1
                if(int(n, int64) > (s - 2_int64) / 2_int64 .or. y < b *     &
                    tol_qp) exit
            end do
            if(y >= b * tol_qp .and. mod(s, 2_int64) /= 0_int64)         &
                       b = b + d * c / a
            g = ((-1) ** s * exp(-a + log_gamma(real(s, qp)) - (s - 1) *   &
                log(a)) + b ) / a
        endif
        res = g
    end function gpx_iint64qp




		  ! single precision use double variables
    impure elemental function ingamma_low_iint8sp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint8sp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint8dp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint8dp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint8qp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint8qp
		  ! single precision use double variables
    impure elemental function ingamma_low_iint16sp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint16sp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint16dp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint16dp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint16qp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint16qp
		  ! single precision use double variables
    impure elemental function ingamma_low_iint32sp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint32sp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint32dp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint32dp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint32qp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint32qp
		  ! single precision use double variables
    impure elemental function ingamma_low_iint64sp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint64sp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint64dp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint64dp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_iint64qp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(ingamma_low): Lower"           &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_iint64qp
		  ! single precision use double variables
    impure elemental function ingamma_low_rspsp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        real(sp), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0._sp) call error_stop("Error(ingamma_low): Lower"          &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_rspsp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_rdpdp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        real(dp), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0._dp) call error_stop("Error(ingamma_low): Lower"          &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_rdpdp
          ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_low_rqpqp(s, x)          &
        result(res)
    !
    ! Approximation of lower incomplete gamma function.
    !
        real(qp), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0._qp) call error_stop("Error(ingamma_low): Lower"          &
            //" incomplete gamma function input s value must be greater than 0")
        
        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = -xx + ss * log(-xx)
            res = (-1) ** s * gpx(s,x) * exp(s1)
        endif
    end function ingamma_low_rqpqp



		    ! single precision use double variables
    impure elemental function l_ingamma_low_iint8sp(s, x)        &
        result(res)

        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0_sp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint8sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint8dp(s, x)        &
        result(res)

        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0_dp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint8dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint8qp(s, x)        &
        result(res)

        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint8qp

		    ! single precision use double variables
    impure elemental function l_ingamma_low_iint16sp(s, x)        &
        result(res)

        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0_sp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint16sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint16dp(s, x)        &
        result(res)

        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0_dp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint16dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint16qp(s, x)        &
        result(res)

        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint16qp

		    ! single precision use double variables
    impure elemental function l_ingamma_low_iint32sp(s, x)        &
        result(res)

        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0_sp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint32sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint32dp(s, x)        &
        result(res)

        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0_dp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint32dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint32qp(s, x)        &
        result(res)

        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint32qp

		    ! single precision use double variables
    impure elemental function l_ingamma_low_iint64sp(s, x)        &
        result(res)

        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0_sp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint64sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint64dp(s, x)        &
        result(res)

        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0_dp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint64dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_iint64qp(s, x)        &
        result(res)

        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(l_ingamma_low): Logrithm of "  &
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_iint64qp

		    ! single precision use double variables
    impure elemental function l_ingamma_low_rspsp(s, x)        &
        result(res)

        real(sp), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0.0_sp) call error_stop("Error(l_ingamma_low): Logrithm of "&
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0_sp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_rspsp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_rdpdp(s, x)        &
        result(res)

        real(dp), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0.0_dp) call error_stop("Error(l_ingamma_low): Logrithm of "&
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0_dp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_rdpdp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_low_rqpqp(s, x)        &
        result(res)

        real(qp), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0.0_qp) call error_stop("Error(l_ingamma_low): Logrithm of "&
            //"lower incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_low): Logrithm of"&
            //" lower incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        elseif(x > real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        endif
    end function l_ingamma_low_rqpqp




		    ! single precision use double variables
    impure elemental function ingamma_up_iint8sp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_dp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint8sp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint8dp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint8dp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint8qp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint8qp

		    ! single precision use double variables
    impure elemental function ingamma_up_iint16sp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_dp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint16sp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint16dp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint16dp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint16qp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint16qp

		    ! single precision use double variables
    impure elemental function ingamma_up_iint32sp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_dp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint32sp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint32dp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint32dp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint32qp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint32qp

		    ! single precision use double variables
    impure elemental function ingamma_up_iint64sp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_dp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint64sp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint64dp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint64dp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_iint64qp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(ingamma_up): Upper"            &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_iint64qp

		    ! single precision use double variables
    impure elemental function ingamma_up_rspsp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        real(sp), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0._sp) call error_stop("Error(ingamma_up): Upper"           &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_sp) then
            res = 0.0
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_dp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_rspsp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_rdpdp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        real(dp), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0._dp) call error_stop("Error(ingamma_up): Upper"           &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_dp) then
            res = 0.0
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_rdpdp

            ! double/quadruple precision use quadruple variables
    impure elemental function ingamma_up_rqpqp(s, x)           &
        result(res)
    !
    ! Approximation of upper incomplete gamma function
    !
        real(qp), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0._qp) call error_stop("Error(ingamma_up): Upper"           &
          //" incomplete gamma function input s value must be greater than 0")

        xx = x; ss = s
        if(x == 0.0_qp) then
            res = 0.0
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = gpx(s,x) * exp(s1)
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = exp(s1 + log(y))
        else
            s1 = log_gamma(ss)
            y = -xx + ss * log(-xx) - s1
            y = 1.0_qp - (-1) ** s * gpx(s,x) * exp(y)
            res = exp(s1 + log(y))
        endif
    end function ingamma_up_rqpqp




		    ! single precision use double variables
    impure elemental function l_ingamma_up_iint8sp(s, x)          &
        result(res)

        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint8sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint8dp(s, x)          &
        result(res)

        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint8dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint8qp(s, x)          &
        result(res)

        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int8) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint8qp

		    ! single precision use double variables
    impure elemental function l_ingamma_up_iint16sp(s, x)          &
        result(res)

        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint16sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint16dp(s, x)          &
        result(res)

        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint16dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint16qp(s, x)          &
        result(res)

        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int16) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint16qp

		    ! single precision use double variables
    impure elemental function l_ingamma_up_iint32sp(s, x)          &
        result(res)

        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint32sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint32dp(s, x)          &
        result(res)

        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint32dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint32qp(s, x)          &
        result(res)

        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int32) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint32qp

		    ! single precision use double variables
    impure elemental function l_ingamma_up_iint64sp(s, x)          &
        result(res)

        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint64sp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint64dp(s, x)          &
        result(res)

        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint64dp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_iint64qp(s, x)          &
        result(res)

        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0_int64) call error_stop("Error(l_ingamma_up): Logrithm of "   &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_iint64qp

		    ! single precision use double variables
    impure elemental function l_ingamma_up_rspsp(s, x)          &
        result(res)

        real(sp), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, y, xx, ss

        if(s < 0.0_sp) call error_stop("Error(l_ingamma_up): Logrithm of " &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_sp .and. x <= real(s, sp)) then
            s1 = log_gamma(ss)
            y = 1.0_dp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, sp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_rspsp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_rdpdp(s, x)          &
        result(res)

        real(dp), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0.0_dp) call error_stop("Error(l_ingamma_up): Logrithm of " &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_dp .and. x <= real(s, dp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, dp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_rdpdp

            ! double/quadruple precision use quadruple variables
    impure elemental function l_ingamma_up_rqpqp(s, x)          &
        result(res)

        real(qp), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, y, xx, ss

        if(s < 0.0_qp) call error_stop("Error(l_ingamma_up): Logrithm of " &
            //"upper incomplete gamma function s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(l_ingamma_up): Logrithm of"  &
            //" upper incomplete gamma function is not defined at x < 0")

        xx = x; ss = s
        if(x > 0.0_qp .and. x <= real(s, qp)) then
            s1 = log_gamma(ss)
            y = 1.0_qp - exp(-xx + ss * log(xx) - s1) * gpx(s,x)
            res = s1 + log(y)
        elseif(x > real(s, qp)) then
            s1 = -xx + ss * log(xx)
            res = log(gpx(s,x)) + s1
        endif
    end function l_ingamma_up_rqpqp




		    ! single precision use double variables
    impure elemental function regamma_p_iint8sp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int8) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 0.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint8sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint8dp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int8) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 0.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint8dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint8qp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int8) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint8qp

		    ! single precision use double variables
    impure elemental function regamma_p_iint16sp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int16) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 0.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint16sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint16dp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int16) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 0.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint16dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint16qp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int16) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint16qp

		    ! single precision use double variables
    impure elemental function regamma_p_iint32sp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int32) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 0.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint32sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint32dp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int32) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 0.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint32dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint32qp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int32) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint32qp

		    ! single precision use double variables
    impure elemental function regamma_p_iint64sp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int64) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 0.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint64sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint64dp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int64) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 0.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint64dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_iint64qp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int64) call error_stop("Error(regamma_p): Regularized gamma" &
            //"_p function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_iint64qp

		    ! single precision use double variables
    impure elemental function regamma_p_rspsp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        real(sp), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0.0_sp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 0.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_rspsp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_rdpdp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        real(dp), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0.0_dp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 0.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_rdpdp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_p_rqpqp(s, x) result(res)
    !
    ! Approximation of regularized incomplet gamma function P(s,x)
    !
        real(qp), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0.0_qp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_p): Regularized "    &
            //"gamma_p function is not defined at x < 0")
        
        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        elseif(x > real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        endif
    end function regamma_p_rqpqp




		    ! single precision use double variables
    impure elemental function regamma_q_iint8sp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int8), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int8) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 1.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        elseif(x > real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint8sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint8dp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int8), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int8) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 1.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint8dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint8qp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int8), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int8) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 1.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint8qp

		    ! single precision use double variables
    impure elemental function regamma_q_iint16sp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int16), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int16) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 1.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        elseif(x > real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint16sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint16dp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int16), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int16) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 1.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint16dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint16qp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int16), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int16) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 1.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint16qp

		    ! single precision use double variables
    impure elemental function regamma_q_iint32sp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int32), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int32) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 1.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        elseif(x > real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint32sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint32dp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int32), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int32) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 1.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint32dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint32qp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int32), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int32) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 1.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint32qp

		    ! single precision use double variables
    impure elemental function regamma_q_iint64sp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int64), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0_int64) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 1.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        elseif(x > real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint64sp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint64dp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int64), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int64) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 1.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint64dp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_iint64qp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        integer(int64), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0_int64) call error_stop("Error(regamma_q): Regularized gamma" &
            //"_q function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 1.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_iint64qp

		    ! single precision use double variables
    impure elemental function regamma_q_rspsp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        real(sp), intent(in) :: s
        real(sp), intent(in) :: x
        real(sp) :: res
        real(dp) :: s1, xx, ss

        if(s < 0.0_sp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function input s value must be non-negative")
        if(x < 0.0_sp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_sp) then
            res = 1.0_dp
        elseif(x > 0.0_sp .and. x <= real(s, sp)) then
            res = 1.0_dp - exp(s1) * gpx(s,x)
        elseif(x > real(s, sp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_rspsp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_rdpdp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        real(dp), intent(in) :: s
        real(dp), intent(in) :: x
        real(dp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0.0_dp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function input s value must be non-negative")
        if(x < 0.0_dp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_dp) then
            res = 1.0_qp
        elseif(x > 0.0_dp .and. x <= real(s, dp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, dp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_rdpdp

            ! double/quadruple precision use quadruple variables
    impure elemental function regamma_q_rqpqp(s, x)             &
        result(res)
    !
    ! Approximation of regularized incomplet gamma function Q(s,x)
    !
        real(qp), intent(in) :: s
        real(qp), intent(in) :: x
        real(qp) :: res
        real(qp) :: s1, xx, ss

        if(s < 0.0_qp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function input s value must be non-negative")
        if(x < 0.0_qp) call error_stop("Error(regamma_q): Regularized "    &
            //"gamma_q function is not defined at x < 0")

        xx = x; ss = s
        s1 = -xx + ss * log(xx) - log_gamma(ss)
        if(x == 0.0_qp) then
            res = 1.0_qp
        elseif(x > 0.0_qp .and. x <= real(s, qp)) then
            res = 1.0_qp - exp(s1) * gpx(s,x)
        elseif(x > real(s, qp)) then
            res = exp(log(gpx(s,x)) + s1)
        endif
     end function regamma_q_rqpqp




    impure elemental function beta_rsp(a, b) result(res)
    !
    ! Evaluation of beta function through gamma function
    !
        real(sp), intent(in) :: a, b
        real(sp) :: res

        if(a <= 0._sp .or. b <= 0._sp) call error_stop("Error(beta):"  &
            //" Beta function auguments a, b values must be greater than 0")
        res = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function beta_rsp

    impure elemental function beta_rdp(a, b) result(res)
    !
    ! Evaluation of beta function through gamma function
    !
        real(dp), intent(in) :: a, b
        real(dp) :: res

        if(a <= 0._dp .or. b <= 0._dp) call error_stop("Error(beta):"  &
            //" Beta function auguments a, b values must be greater than 0")
        res = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function beta_rdp

    impure elemental function beta_rqp(a, b) result(res)
    !
    ! Evaluation of beta function through gamma function
    !
        real(qp), intent(in) :: a, b
        real(qp) :: res

        if(a <= 0._qp .or. b <= 0._qp) call error_stop("Error(beta):"  &
            //" Beta function auguments a, b values must be greater than 0")
        res = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function beta_rqp

    impure elemental function beta_csp(a, b) result(res)
    !
    ! Evaluation of beta function through gamma function
    !
        complex(sp), intent(in) :: a, b
        complex(sp) :: res

		if(a % re <= 0._sp .or. b % re <= 0._sp)                       &
		    call error_stop("Error(beta): Beta function complex auguments"     &
			   //" a, b must have a positive real part")
        res = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function beta_csp

    impure elemental function beta_cdp(a, b) result(res)
    !
    ! Evaluation of beta function through gamma function
    !
        complex(dp), intent(in) :: a, b
        complex(dp) :: res

		if(a % re <= 0._dp .or. b % re <= 0._dp)                       &
		    call error_stop("Error(beta): Beta function complex auguments"     &
			   //" a, b must have a positive real part")
        res = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function beta_cdp

    impure elemental function beta_cqp(a, b) result(res)
    !
    ! Evaluation of beta function through gamma function
    !
        complex(qp), intent(in) :: a, b
        complex(qp) :: res

		if(a % re <= 0._qp .or. b % re <= 0._qp)                       &
		    call error_stop("Error(beta): Beta function complex auguments"     &
			   //" a, b must have a positive real part")
        res = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function beta_cqp




    impure elemental function l_beta_rsp(a, b) result(res)
    !
    ! Logrithm of beta function through log(gamma)
    !
        real(sp), intent(in) :: a, b
        real(sp) :: res

        if(a <= 0._sp .or. b <= 0._sp) call error_stop("Error(l_beta):"&
            //" Beta function auguments a, b values must be greater than 0")
        res = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function l_beta_rsp

    impure elemental function l_beta_rdp(a, b) result(res)
    !
    ! Logrithm of beta function through log(gamma)
    !
        real(dp), intent(in) :: a, b
        real(dp) :: res

        if(a <= 0._dp .or. b <= 0._dp) call error_stop("Error(l_beta):"&
            //" Beta function auguments a, b values must be greater than 0")
        res = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function l_beta_rdp

    impure elemental function l_beta_rqp(a, b) result(res)
    !
    ! Logrithm of beta function through log(gamma)
    !
        real(qp), intent(in) :: a, b
        real(qp) :: res

        if(a <= 0._qp .or. b <= 0._qp) call error_stop("Error(l_beta):"&
            //" Beta function auguments a, b values must be greater than 0")
        res = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function l_beta_rqp

    impure elemental function l_beta_csp(a, b) result(res)
    !
    ! Logrithm of beta function through log(gamma)
    !
        complex(sp), intent(in) :: a, b
        complex(sp) :: res

		if(a % re <= 0._sp .or. b % re <= 0._sp)                       &
		    call error_stop("Error(l_beta): log_beta function complex"         &
			   //" auguments a, b must have a positive real part")
        res = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function l_beta_csp

    impure elemental function l_beta_cdp(a, b) result(res)
    !
    ! Logrithm of beta function through log(gamma)
    !
        complex(dp), intent(in) :: a, b
        complex(dp) :: res

		if(a % re <= 0._dp .or. b % re <= 0._dp)                       &
		    call error_stop("Error(l_beta): log_beta function complex"         &
			   //" auguments a, b must have a positive real part")
        res = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function l_beta_cdp

    impure elemental function l_beta_cqp(a, b) result(res)
    !
    ! Logrithm of beta function through log(gamma)
    !
        complex(qp), intent(in) :: a, b
        complex(qp) :: res

		if(a % re <= 0._qp .or. b % re <= 0._qp)                       &
		    call error_stop("Error(l_beta): log_beta function complex"         &
			   //" auguments a, b must have a positive real part")
        res = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function l_beta_cqp




	    ! single precision use double variables
    impure elemental function inbeta_rsp(x, a, b) result(res)
    !
    ! Evaluation of incomplete beta function using continued fractions
    ! "Computation of Special Functions" by S. Zhang and J. Jin, 1996
    !
        real(sp), intent(in) :: x, a, b
        real(sp) :: res, s0
        integer :: n, k
        real(dp) :: an, bn, g, c, d, y, ak, ak2, dm

        dm = dmd
        if(a <= 0._sp .or. b <= 0._sp) call error_stop("Error(inbeta):"&
            //" Incomplete beta function auguments a, b must be positive")
        if(x < 0.0_sp .or. x > 1.0_sp) call error_stop("Error(inbeta):"&
            //" Incomplete beta function augument x must be in range [0,1]")
        if(x == 0.0_sp) then
            res = 0.0_sp
            return
        elseif(x == 1.0_sp) then
            res = beta(a, b)
            return
        end if
        s0 = (a + 1) / (a + b + 2)
        an = 1.0_dp
        bn = 1.0_dp
        g = an / bn
        c = an / dm
        d = 1.0_dp / bn
        n = 1
        if(x < s0) then
            do
                if(mod(n, 2) == 0) then
                    k = n / 2; ak = a + 2 * k
                    an = real(k, dp) * x * (b - k) / (ak * ak - ak)
                else
                    k = (n - 1) / 2; ak = a + k; ak2 = ak + k
                    an = - (ak + b) * ak * x / (ak2 * ak2 + ak2)
                endif
                d = d * an + bn
                if(d == 0.0_dp) d = dm
                c = bn + an / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
            g = exp(a * log(x) + b * log(1.0_dp - x) + log(g) -log(a))
        else
            do
                if(mod(n, 2) == 0) then
                    k = n / 2; ak = b + 2 * k
                    an = k * (1.0_dp - x) * (a - k) / (ak * ak - ak)
                else
                    k = (n - 1) / 2; ak = b + k; ak2 = ak + k
                    an = - ak * (1.0_dp - x) * (a + ak) / (ak2 * ak2 + ak2)
                endif
                d = d * an + bn
                if(d == 0.0_dp) d = dm
                c = bn + an / c
                if(c == 0.0_dp) c = dm
                d = 1.0_dp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_dp) < tol_dp) exit
            end do
            g = a * log(x) + b * log(1.0_dp - x) + log(g) - log(b)
            g = beta(a, b) - exp(g)
        endif
        res = g
    end function inbeta_rsp

	    ! double/quadruple precision use quadruple variables
    impure elemental function inbeta_rdp(x, a, b) result(res)
    !
    ! Evaluation of incomplete beta function using continued fractions
    ! "Computation of Special Functions" by S. Zhang and J. Jin, 1996
    !
        real(dp), intent(in) :: x, a, b
        real(dp) :: res, s0
        integer :: n, k
        real(qp) :: an, bn, g, c, d, y, ak, ak2, dm

        dm = dmq
        if(a <= 0._dp .or. b <= 0._dp) call error_stop("Error(inbeta):"&
            //" Incomplete beta function auguments a, b must be positive")
        if(x < 0.0_dp .or. x > 1.0_dp) call error_stop("Error(inbeta):"&
            //" Incomplete beta function augument x must be in range [0,1]")
        if(x == 0.0_dp) then
            res = 0.0_dp
            return
        elseif(x == 1.0_dp) then
            res = beta(a, b)
            return
        end if
        s0 = (a + 1) / (a + b + 2)
        an = 1.0_qp
        bn = 1.0_qp
        g = an / bn
        c = an / dm
        d = 1.0_qp / bn
        n = 1
        if(x < s0) then
            do
                if(mod(n, 2) == 0) then
                    k = n / 2; ak = a + 2 * k
                    an = real(k, qp) * x * (b - k) / (ak * ak - ak)
                else
                    k = (n - 1) / 2; ak = a + k; ak2 = ak + k
                    an = - (ak + b) * ak * x / (ak2 * ak2 + ak2)
                endif
                d = d * an + bn
                if(d == 0.0_qp) d = dm
                c = bn + an / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
            g = exp(a * log(x) + b * log(1.0_qp - x) + log(g) -log(a))
        else
            do
                if(mod(n, 2) == 0) then
                    k = n / 2; ak = b + 2 * k
                    an = k * (1.0_qp - x) * (a - k) / (ak * ak - ak)
                else
                    k = (n - 1) / 2; ak = b + k; ak2 = ak + k
                    an = - ak * (1.0_qp - x) * (a + ak) / (ak2 * ak2 + ak2)
                endif
                d = d * an + bn
                if(d == 0.0_qp) d = dm
                c = bn + an / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
            g = a * log(x) + b * log(1.0_qp - x) + log(g) - log(b)
            g = beta(a, b) - exp(g)
        endif
        res = g
    end function inbeta_rdp

	    ! double/quadruple precision use quadruple variables
    impure elemental function inbeta_rqp(x, a, b) result(res)
    !
    ! Evaluation of incomplete beta function using continued fractions
    ! "Computation of Special Functions" by S. Zhang and J. Jin, 1996
    !
        real(qp), intent(in) :: x, a, b
        real(qp) :: res, s0
        integer :: n, k
        real(qp) :: an, bn, g, c, d, y, ak, ak2, dm

        dm = dmq
        if(a <= 0._qp .or. b <= 0._qp) call error_stop("Error(inbeta):"&
            //" Incomplete beta function auguments a, b must be positive")
        if(x < 0.0_qp .or. x > 1.0_qp) call error_stop("Error(inbeta):"&
            //" Incomplete beta function augument x must be in range [0,1]")
        if(x == 0.0_qp) then
            res = 0.0_qp
            return
        elseif(x == 1.0_qp) then
            res = beta(a, b)
            return
        end if
        s0 = (a + 1) / (a + b + 2)
        an = 1.0_qp
        bn = 1.0_qp
        g = an / bn
        c = an / dm
        d = 1.0_qp / bn
        n = 1
        if(x < s0) then
            do
                if(mod(n, 2) == 0) then
                    k = n / 2; ak = a + 2 * k
                    an = real(k, qp) * x * (b - k) / (ak * ak - ak)
                else
                    k = (n - 1) / 2; ak = a + k; ak2 = ak + k
                    an = - (ak + b) * ak * x / (ak2 * ak2 + ak2)
                endif
                d = d * an + bn
                if(d == 0.0_qp) d = dm
                c = bn + an / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
            g = exp(a * log(x) + b * log(1.0_qp - x) + log(g) -log(a))
        else
            do
                if(mod(n, 2) == 0) then
                    k = n / 2; ak = b + 2 * k
                    an = k * (1.0_qp - x) * (a - k) / (ak * ak - ak)
                else
                    k = (n - 1) / 2; ak = b + k; ak2 = ak + k
                    an = - ak * (1.0_qp - x) * (a + ak) / (ak2 * ak2 + ak2)
                endif
                d = d * an + bn
                if(d == 0.0_qp) d = dm
                c = bn + an / c
                if(c == 0.0_qp) c = dm
                d = 1.0_qp / d
                y = c * d
                g = g * y
                n = n + 1
                if(abs(y - 1.0_qp) < tol_qp) exit
            end do
            g = a * log(x) + b * log(1.0_qp - x) + log(g) - log(b)
            g = beta(a, b) - exp(g)
        endif
        res = g
    end function inbeta_rqp




    impure elemental function regbeta_rsp(x, a, b) result(res)

        real(sp), intent(in) :: x, a, b
        real(sp) :: res

        if(a <= 0._sp .or. b <= 0._sp) call error_stop("Error(regbeta)"&
            //": Regularized incomplete beta function auguments a, b must be " &
            //"positive")
        if(x < 0.0_sp .or. x > 1.0_sp) call error_stop("Error(regbeta)"&
            //": Regularized incomplete beta function augument x must be in "  &
            //"range [0,1]")
        if(x == 0.0_sp) then
            res = 0.0_sp
        elseif(x == 1.0_sp) then
            res = 1.0_sp
        else
            res = incomplete_beta(x, a, b) / beta(a, b)
        end if
    end function regbeta_rsp

    impure elemental function regbeta_rdp(x, a, b) result(res)

        real(dp), intent(in) :: x, a, b
        real(dp) :: res

        if(a <= 0._dp .or. b <= 0._dp) call error_stop("Error(regbeta)"&
            //": Regularized incomplete beta function auguments a, b must be " &
            //"positive")
        if(x < 0.0_dp .or. x > 1.0_dp) call error_stop("Error(regbeta)"&
            //": Regularized incomplete beta function augument x must be in "  &
            //"range [0,1]")
        if(x == 0.0_dp) then
            res = 0.0_dp
        elseif(x == 1.0_dp) then
            res = 1.0_dp
        else
            res = incomplete_beta(x, a, b) / beta(a, b)
        end if
    end function regbeta_rdp

    impure elemental function regbeta_rqp(x, a, b) result(res)

        real(qp), intent(in) :: x, a, b
        real(qp) :: res

        if(a <= 0._qp .or. b <= 0._qp) call error_stop("Error(regbeta)"&
            //": Regularized incomplete beta function auguments a, b must be " &
            //"positive")
        if(x < 0.0_qp .or. x > 1.0_qp) call error_stop("Error(regbeta)"&
            //": Regularized incomplete beta function augument x must be in "  &
            //"range [0,1]")
        if(x == 0.0_qp) then
            res = 0.0_qp
        elseif(x == 1.0_qp) then
            res = 1.0_qp
        else
            res = incomplete_beta(x, a, b) / beta(a, b)
        end if
    end function regbeta_rqp

end module statislib_special
