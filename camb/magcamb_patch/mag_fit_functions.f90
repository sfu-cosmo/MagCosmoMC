!---------------------------------------------------------------------------
!
!> This module contains the fitting functions of the magnetic power spectra
!> the fitting functions are defined as t^(2n+6) * int where t = k/k_D and int is the integral
!
!   @author: Alex Zucca, azucca@sfu.ca
!   @author: Yun Li, yun_li_3@sfu.ca
!
!   MagCAMB:
!   version 1.0 - 2016: introducing fits.
!   version 2.0 - 2019: created functions for the fitting function
!                       adding the helical contribution
!
!---------------------------------------------------------------------------

module MagCAMB_fitting_functions

    !> common parameters
    use Precision
    implicit none


contains

    !---------------------------------------------------------------------------
    !> this functions computes the scalar delta-delta integral for non-helical modes
    function MagInt_Scalar_DeltaDelta(x,y)
        implicit none

        real(dl) :: MagInt_Scalar_DeltaDelta
        real(dl) :: x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        real(dl) :: lx,lx2, lx3, lx4, lx5, lx6, lx7, lx8, lx9
        !T VS P, N=-2.05
        !    f(x) = p1*x9 + p2*x8 + p3*x7 + p4*x6 +
        !                  p5*x5 + p6*x4 + p7*x3 + p8*x2 + p9*x + p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: p1 =   1.609e+08  ! (1.545e+08, 1.672e+08)
        real(dl) :: p2 =  -1.176e+08  !(-1.219e+08, -1.133e+08)
        real(dl) :: p3 =   3.627e+07  !(3.508e+07, 3.746e+07)
        real(dl) :: p4 =  -6.151e+06  !(-6.329e+06, -5.973e+06)
        real(dl) :: p5 =   6.279e+05  !(6.123e+05, 6.436e+05)
        real(dl) :: p6 =   -3.99e+04  !(-4.071e+04, -3.908e+04)
        real(dl) :: p7 =        1601  !(1576, 1625)
        real(dl) :: p8 =         -45  !(-45.39, -44.61)
        real(dl) :: p9 =      -1.983  !(-1.986, -1.98)
        real(dl) :: p10 =       7.329 ! (7.329, 7.329)

        !Goodness of fit:
        ! SSE: 4.804e-06
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 3.821e-05

        !!!!!!!!!!!!! n=-1.55, !!!x=logt vs y= p
        !Linear model Poly9:
        !   f(x) = sp1*lx9 + sp2*lx8 + sp3*lx7 + sp4*lx6 +
        !                 sp5*lx5 + sp6*lx4 + sp7*lx3 + sp8*lx2 + sp9*lx + sp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sp1 =  -5.917e-10  !(-5.963e-10, -5.871e-10)
        real(dl) :: sp2 =   -5.88e-08  !(-5.92e-08, -5.839e-08)
        real(dl) :: sp3 =   -2.53e-06  !(-2.545e-06, -2.515e-06)
        real(dl) :: sp4 =  -6.195e-05  !(-6.225e-05, -6.164e-05)
        real(dl) :: sp5 =  -0.0009543  !(-0.000958, -0.0009505)
        real(dl) :: sp6 =   -0.009714  !(-0.009742, -0.009685)
        real(dl) :: sp7 =    -0.06965  !(-0.06978, -0.06951)
        real(dl) :: sp8 =     -0.4654  !(-0.4658, -0.4651)
        real(dl) :: sp9 =      -4.656  !(-4.656, -4.655)
        real(dl) :: sp10 =      0.8963  !(0.8959, 0.8967)

        !Goodness of fit:
        ! SSE: 2.21e-06
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 2.597e-05

        !!!!!!!N=-0.005
        !!!! t<0.00164
        !Linear model Poly9:
        !  f(x) = sddp1*x9 + sddp2*x8 + sddp3*x7 + sddp4*x6 +
        !                sddp5*x5 + sddp6*x4 + sddp7*x3 + sddp8*x2 + sddp9*x + sddp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sddp1 =  -1.428e+13 ! (-1.023e+14, 7.37e+13)
        real(dl) :: sddp2 =   1.109e+11  !(-5.357e+11, 7.575e+11)
        real(dl) :: sddp3 =  -3.571e+08  !(-2.347e+09, 1.632e+09)
        real(dl) :: sddp4 =   6.151e+05  !(-2.709e+06, 3.94e+06)
        real(dl) :: sddp5 =      -612.1  !(-3883, 2659)
        real(dl) :: sddp6 =     -0.6486  !(-2.568, 1.271)
        real(dl) :: sddp7 =       1.338  !(1.337, 1.338)
        real(dl) :: sddp8 =   1.756e-08  !(-9.613e-08, 1.313e-07)
        real(dl) :: sddp9 =  -1.014e-12  !(-9.119e-12, 7.091e-12)
        real(dl) :: sddp10 =   4.967e-18  !(-8.722e-17, 9.715e-17)

        !Goodness of fit:
        ! SSE: 3.349e-30
        !R-square: 1
        !djusted R-square: 1
        !RMSE: 1.878e-16

        !t>0.00164
        !Linear model Poly9:
        !  f(x) = sdd2p1*x9 + sdd2p2*x8 + sdd2p3*x7 + sdd2p4*x6 +
        !                sdd2p5*x5 + sdd2p6*x4 + sdd2p7*x3 + sdd2p8*x2 + sdd2p9*x + sdd2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sdd2p1 =  -1.591e+04 ! (-3.134e+04, -473.3)
        real(dl) :: sdd2p2 =        4245  !(-38.94, 8529)
        real(dl) :: sdd2p3 =      -474.6  !(-976.8, 27.66)
        real(dl) :: sdd2p4 =       30.42  !(-1.943, 62.79)
        real(dl) :: sdd2p5 =      -2.386  !(-3.636, -1.136)
        real(dl) :: sdd2p6 =     -0.9775  !(-1.007, -0.9479)
        real(dl) :: sdd2p7 =       1.338  !(1.337, 1.338)
        real(dl) :: sdd2p8 =   2.053e-06  !(-1.336e-06, 5.443e-06)
        real(dl) :: sdd2p9 =  -7.035e-09  !(-2.057e-08, 6.503e-09)
        real(dl) :: sdd2p10 =   8.407e-12  !(-1.117e-11, 2.798e-11)
        !Goodness of fit:
        ! SSE: 7.367e-19
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.591e-11
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        lx = log(x)
        lx2 = lx*lx
        lx3 = lx*lx2
        lx4 = lx*lx3
        lx5 = lx*lx4
        lx6 = lx*lx5
        lx7 = lx*lx6
        lx8 = lx*lx7
        lx9 = lx*lx8
        !!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.5) then !n=-2.05
                MagInt_Scalar_DeltaDelta = (x**(1.9_dl))*(p1*x9 + p2*x8 + p3*x7 + p4*x6 &
                + p5*x5 + p6*x4 + p7*x3 + p8*x2 + p9*x + p10)
                !write(*,*) "non-he,t, MagInt_Scalar_DeltaDelta  = ", x, MagInt_Scalar_DeltaDelta  ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.5) then !n=-1.55
                MagInt_Scalar_DeltaDelta = (x**(2.9_dl))*( sp1*lx9 + sp2*lx8 + sp3*lx7 + sp4*lx6 +&
                sp5*lx5 + sp6*lx4 + sp7*lx3 + sp8*lx2 + sp9*lx + sp10)
                !write(*,*) "non-he,t, MagInt_Scalar_DeltaDelta  = ", x, MagInt_Scalar_DeltaDelta  ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaDelta = (x**(2._dl*y+6._dl))*((4*x)/3 + 4/x + x2/4 - 5)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.00164) then
                MagInt_Scalar_DeltaDelta = sddp1*x9 + sddp2*x8 + sddp3*x7 + sddp4*x6 +&
                sddp5*x5 + sddp6*x4 + sddp7*x3 + sddp8*x2 + sddp9*x + sddp10
                !write(*,*) "non-he,t, MagInt_Scalar_DeltaDelta  = ", x, MagInt_Scalar_DeltaDelta  ! good
            else if ( x .GE. 0.00164 .and. x .le. 0.5) then
                MagInt_Scalar_DeltaDelta = sdd2p1*x9 + sdd2p2*x8 + sdd2p3*x7 + sdd2p4*x6 + &
                sdd2p5*x5 + sdd2p6*x4 + sdd2p7*x3 + sdd2p8*x2 + sdd2p9*x + sdd2p10
                !write(*,*) "non-he,2, MagInt_Scalar_DeltaDelta  = ", x, MagInt_Scalar_DeltaDelta  ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaDelta = (x**(2._dl*y+6._dl))*(4/(15*x) + 1/(4*x2) - 1/x4 + 4/(5*x5) - 1/5)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaDelta = (x**(2._dl*y+6._dl))*(8/(15*x5) - 1/(24*x2) - 1/x6 + 4/(7*x7) + 11/2240)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaDelta = (x**(2._dl*y+6._dl))*(4/(315*x3) + 4/(75*x5) - 5/(12*x6) + 20/(21*x7)&
                - 1/x8 + 4/(9*x9) - 1/525)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Scalar_DeltaDelta


    !---------------------------------------------------------------------------
    !> this functions computes the scalar delta-delta integral for helical modes
    function MagInt_Scalar_Helical_DeltaDelta(x,y)

        implicit none

        !real(dl), parameter :: Pi=3.14159265d0
        real(dl) :: MagInt_Scalar_Helical_DeltaDelta, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        real(dl) :: lx,lx2, lx3, lx4, lx5, lx6, lx7, lx8, lx9
        !! N=-1.5     !!!! logt vs p
        !Linear model Poly9:
        !   f(x) = hddp1*lx9 + hddp2*lx8 + hddp3*lx7 + hddp4*lx6 +
        !                 hddp5*lx5 + hddp6*lx4 + hddp7*lx3 + hddp8*lx2 + hddp9*lx + hddp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hddp1 =   1.526e-10  !(1.516e-10, 1.535e-10)
        real(dl) :: hddp2 =   1.594e-08  !(1.585e-08, 1.603e-08)
        real(dl) :: hddp3 =   7.254e-07  !(7.218e-07, 7.289e-07)
        real(dl) :: hddp4 =   1.889e-05  !(1.881e-05, 1.897e-05)
        real(dl) :: hddp5 =   0.0003108  !(0.0003098, 0.0003119)
        real(dl) :: hddp6 =    0.003363  !(0.003354, 0.003372)
        real(dl) :: hddp7 =     0.02406  !(0.02401, 0.02411)
        real(dl) :: hddp8 =      0.1105  !(0.1104, 0.1107)
        real(dl) :: hddp9 =       2.299  !(2.299, 2.299)
        real(dl) :: hddp10 =      0.4525  !(0.4523, 0.4527)

        !Goodness of fit:
        ! SSE: 1.903e-07
        ! R-square: 1
        !Adjusted R-square: 1
        !RMSE: 5.544e-06

        !!!!!!!!!!!!! n=-1.005
        !!!!!! t<0.0014
        !General model:
        !  f(x) = hhddp4*x6 + hhddp5*x5 + hhddp6*x4 + hhddp7*x3 + hhddp8*x^(3.8)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hhddp4 =      0.1576 ! (0.1211, 0.194)
        real(dl) :: hhddp5 =     -0.7739  !(-0.7821, -0.7656)
        real(dl) :: hhddp6 =       2.949  !(2.947, 2.951)
        real(dl) :: hhddp7 =       -2.02  !(-2.02, -2.02)
        real(dl) :: hhddp8 =     0.06198  !(0.06101, 0.06295)

        !Goodness of fit:
        ! SSE: 1.859e-18
        !R-square: 1
        !djusted R-square: 1
        !RMSE: 1.732e-11

        !General model:
        !   f(x) = -hhdda*x3 -hhddb*x**(2.93_dl)
        !Coefficients (with 95% confidence bounds):
        !    real(dl) :: hhdda =      0.7495 ! (-5.639, 7.138)
        !   real(dl) :: hhddb =       0.408  !(-3.612, 4.428)

        !Goodness of fit:
        ! SSE: 1.147e-16
        !R-square: 0.884
        !Adjusted R-square: 0.8837
        !RMSE: 5.685e-10

        !!!!!t> 0.0014
        !Linear model Poly9:
        !    f(x) = hhdd2p1*x9 + hhdd2p2*x8 + hhdd2p3*x7 + hhdd2p4*x6 +
        !                  hhdd2p5*x5 + hhdd2p6*x4 + hhdd2p7*x3 + hhdd2p8*x2 + hhdd2p9*x + hhdd2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hhdd2p1 =       620.7 ! (-1.119e+04, 1.243e+04)
        real(dl) :: hhdd2p2 =       67.15  !(-3198, 3333)
        real(dl) :: hhdd2p3 =      -50.98  !(-432.2, 330.2)
        real(dl) :: hhdd2p4 =       8.308  !(-16.13, 32.75)
        real(dl) :: hhdd2p5 =      -1.542  !(-2.48, -0.605)
        real(dl) :: hhdd2p6 =       3.077  !(3.055, 3.099)
        real(dl) :: hhdd2p7 =       -2.02  !(-2.02, -2.019)
        real(dl) :: hhdd2p8 =  -1.784e-06  !(-4.233e-06, 6.655e-07)
        real(dl) :: hhdd2p9 =   4.724e-09  !(-4.823e-09, 1.427e-08)
        real(dl) :: hhdd2p10 =   -4.74e-12  !(-1.804e-11, 8.56e-12)
        !Goodness of fit:
        ! SSE: 1.857e-18
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.781e-11
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        lx = log(x)
        lx2 = lx*lx
        lx3 = lx*lx2
        lx4 = lx*lx3
        lx5 = lx*lx4
        lx6 = lx*lx5
        lx7 = lx*lx6
        lx8 = lx*lx7
        lx9 = lx*lx8
        !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta = (x**(2._dl*y+6._dl))*(2*x + x2/2 - 2)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta = (x**(2._dl*y+6._dl))*(hddp1*lx9 + hddp2*lx8 + hddp3*lx7 + hddp4*lx6 + &
                hddp5*lx5 + hddp6*lx4 + hddp7*lx3 + hddp8*lx2 + hddp9*lx + hddp10)
                !write(*,*) "helic,t, MagInt_Scalar_Helical_DeltaDelta  = ", x, MagInt_Scalar_Helical_DeltaDelta   !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.0014) then
                MagInt_Scalar_Helical_DeltaDelta = (hhddp4*x6 + hhddp5*x5 + hhddp6*x4 + hhddp7*x3 + hhddp8*x**(3.8_dl))
                ! -hhdda*x3 -hhddb*x**(2.93_dl)
                !write(*,*) "helic,t, MagInt_Scalar_Helical_DeltaDelta  = ", x, MagInt_Scalar_Helical_DeltaDelta   ! good
            else if ( x .GE. 0.0014 .and. x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta =  hhdd2p1*x9 + hhdd2p2*x8 + hhdd2p3*x7 + hhdd2p4*x6 +&
                hhdd2p5*x5 + hhdd2p6*x4 + hhdd2p7*x3 + hhdd2p8*x2 + hhdd2p9*x + hhdd2p10
                !write(*,*) "helic,2, MagInt_Scalar_Helical_DeltaDelta  = ", x, MagInt_Scalar_Helical_DeltaDelta   !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta = (x**(2._dl*y+6._dl))*(2/(3*x) + 1/(2*x2) - 2/(3*x3) - 1/2)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta = (x**(2._dl*y+6._dl))*(1/(2*x4) - 1/(8*x2)&
                - 2/(5*x5) + 1/80)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta = (x**(2._dl*y+6._dl))*(2/(45*x3) - 4/(15*x5)&
                + 1/(2*x6) - 2/(7*x7) - 1/210)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaDelta = (x**(2._dl*y+6._dl))*(5/(24*x6) - 1/(48*x4) &
                - 10/(21*x7) + 1/(2*x8) - 2/(9*x9) + 1/4032)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Scalar_Helical_DeltaDelta


    !---------------------------------------------------------------------------
    !> this functions computes the scalar Pi-Pi integral for non-helical modes
    function MagInt_Scalar_PiPi(x,y)

        implicit none

        real(dl) :: MagInt_Scalar_PiPi, x, y, try !check
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        real(dl) :: lx,lx2, lx3, lx4, lx5, lx6, lx7, lx8, lx9
        ! T VS P , N=-2.05,
        ! t<0.0005,  !!!t vs p !t<0.0005
        !Linear model Poly4:
        !    f(x) = sp1p1*x4 + sp1p2*x3 + sp1p3*x2 + sp1p4*x + sp1p5
        !oefficients (with 95% confidence bounds):
        real(dl) :: sp1p1 =  -9.086e+14 ! (-1.617e+15, -2.005e+14)
        real(dl) :: sp1p2 =   1.142e+12  !(3.913e+11, 1.893e+12)
        real(dl) :: sp1p3 =  -4.718e+08  !(-7.293e+08, -2.144e+08)
        real(dl) :: sp1p4 =   6.421e+04  !(3.334e+04, 9.508e+04)
        real(dl) :: sp1p5 =       10.88  !(9.672, 12.1)
        !Goodness of fit:
        ! SSE: 5.663
        !R-square: 0.7579
        !Adjusted R-square: 0.701
        !RMSE: 0.5772

        !! t vs p,   t>0.0005
        !Linear model Poly9:
        !    f(x) = spp1*x9 + spp2*x8 + spp3*x7 + spp4*x6 +
        !                  spp5*x5 + spp6*x4 + spp7*x3 + spp8*x2 + spp9*x + spp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: spp1 =  -2.314e+12  !(-2.949e+12, -1.679e+12)
        real(dl) :: spp2 =   6.648e+11  !(4.918e+11, 8.379e+11)
        real(dl) :: spp3 =  -8.083e+10  !(-1.007e+11, -6.098e+10)
        real(dl) :: spp4 =    5.41e+09  !(4.165e+09, 6.655e+09)
        real(dl) :: spp5 =  -2.172e+08  !(-2.636e+08, -1.708e+08)
        real(dl) :: spp6 =   5.337e+06  !(4.288e+06, 6.386e+06)
        real(dl) :: spp7 =  -7.831e+04  !(-9.233e+04, -6.429e+04)
        real(dl) :: spp8 =       622.1  !(519, 725.3)
        real(dl) :: spp9 =      -5.692  !(-6.051, -5.333)
        real(dl) :: spp10 =       10.96  !(10.96, 10.96)

        !Goodness of fit:
        ! SSE: 0.00181
        !R-square: 0.9999
        !Adjusted R-square: 0.9999
        !RMSE: 0.0007811

        !!!!!!!!!!!!! n=-1.55,
        !!!!!!! t vs p , t<0.0012  !!t vs p
        !Linear model Poly3:
        !   f(x) = sspp1*x3 + sspp2*x2 + sspp3*x + sspp4
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sspp1 =  -4.336e+10  !(-1.003e+11, 1.355e+10)
        real(dl) :: sspp2 =   1.251e+08  !(1.288e+07, 2.372e+08)
        real(dl) :: sspp3 =  -1.167e+05  !(-1.836e+05, -4.984e+04)
        real(dl) :: sspp4 =       66.03  !(53.14, 78.92)

        !Goodness of fit:
        ! SSE: 2026
        ! R-square: 0.4311
        ! Adjusted R-square: 0.3923
        ! RMSE: 6.785

        !!!!!!logt vs p, t>0.0012
        !Linear model Poly9:
        !   f(x) = p1*x^9 + p2*x^8 + p3*x^7 + p4*x^6 +
        !                 p5*x^5 + p6*x^4 + p7*x^3 + p8*x^2 + p9*x + p10
        !     ssp2p1*lx9 + ssp2p2*lx8 + ssp2p3*lx7 + ssp2p4*lx6 +
        !                ssp2p5*lx5 + ssp2p6*lx4 + ssp2p7*lx3 + ssp2p8*lx2 + ssp2p9*lx + ssp2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: ssp2p1 =   2.596e-05  !(8.279e-06, 4.364e-05)
        real(dl) :: ssp2p2 =   0.0009522  !(0.0002736, 0.001631)
        real(dl) :: ssp2p3 =     0.01525  !(0.003879, 0.02662)
        real(dl) :: ssp2p4 =      0.1399  !(0.0307, 0.249)
        real(dl) :: ssp2p5 =      0.8083  !(0.1474, 1.469)
        real(dl) :: ssp2p6 =       3.046  !(0.4322, 5.66)
        real(dl) :: ssp2p7 =       7.455  !(0.7027, 14.21)
        real(dl) :: ssp2p8 =       11.09  !(0.1203, 22.07)
        real(dl) :: ssp2p9 =        3.97  !(-6.204, 14.14)
        real(dl) :: ssp2p10 =       6.425  !(2.326, 10.52)

        !Goodness of fit:
        ! SSE: 0.01856
        ! R-square: 1
        ! Adjusted R-square: 1
        !RMSE: 0.002439

        !!!!!!!  N=-0.005
        !!!! t< 0.0018
        !General model:
        !   f(x) = ssspp1*x6 + sssp2*x5 + sssp3*x4 + sssp4*x3
        !Coefficients (with 95% confidence bounds):
        real(dl) :: ssspp1 =      0.4116  !(-1.005, 1.828)
        real(dl) :: ssspp2 =     -0.6527  !(-0.8447, -0.4607)
        real(dl) :: ssspp3 =          -1  !(-1.009, -0.9917)
        real(dl) :: ssspp4 =       1.873  !(1.873, 1.873)

        !Goodness of fit:
        ! SSE: 1.612e-14
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 2.331e-09

        !!!! t> 0.0018
        !    f(x) = sssp2p1*x9 + sssp2p2*x8 + sssp2p3*x7 + sssp2p4*x6 +
        !                  sssp2p5*x5 + sssp2p6*x4 + sssp2p7*x3 + sssp2p8*x2 + sssp2p9*x + sssp2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sssp2p1 =  -1.306e+05  !(-2.445e+05, -1.665e+04)
        real(dl) :: sssp2p2 =   3.872e+04  !(6989, 7.045e+04)
        real(dl) :: sssp2p3 =       -4891  !(-8626, -1156)
        real(dl) :: sssp2p4 =       343.9  !(102.1, 585.8)
        real(dl) :: sssp2p5 =      -15.33  !(-24.73, -5.934)
        real(dl) :: sssp2p6 =     -0.6079  !(-0.832, -0.3838)
        real(dl) :: sssp2p7 =       1.866  !(1.863, 1.87)
        real(dl) :: sssp2p8 =   6.341e-05  !(3.707e-05, 8.976e-05)
        real(dl) :: sssp2p9 =  -3.284e-07  !(-4.36e-07, -2.209e-07)
        real(dl) :: sssp2p10 =    6.89e-10  !(5.283e-10, 8.497e-10)

        !Goodness of fit:
        ! SSE: 3.749e-17
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.137e-10
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        lx = log(x)
        lx2 = lx*lx
        lx3 = lx*lx2
        lx4 = lx*lx3
        lx5 = lx*lx4
        lx6 = lx*lx5
        lx7 = lx*lx6
        lx8 = lx*lx7
        lx9 = lx*lx8
        !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.0005) then !n=-2.05
                MagInt_Scalar_PiPi = (x**(1.9_dl))*(sp1p1*x4 + sp1p2*x3 + sp1p3*x2 + sp1p4*x + sp1p5)
                ! write(*,*) "non-h,t, MagInt_Scalar_PiPi  = ", x, MagInt_Scalar_PiPi   ! good
            else if ( x .GE. 0.0005 .and. x .le. 0.5) then
                MagInt_Scalar_PiPi = (x**(1.9_dl))*(spp1*x9 + spp2*x8 + spp3*x7 + spp4*x6 +&
                spp5*x5 + spp6*x4 + spp7*x3 + spp8*x2 + spp9*x + spp10)
                !write(*,*) "non-h,2, MagInt_Scalar_PiPi  = ", x, MagInt_Scalar_PiPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.0012) then !n=-1.55,  2.9_dl, 2._dl*y+6._dl
                MagInt_Scalar_PiPi = (x**(2.9_dl))*(sspp1*x3 + sspp2*x2 + sspp3*x + sspp4)
                !write(*,*) "non-he,t, lx,MagInt_Scalar_PiPi  = ", x, MagInt_Scalar_PiPi   ! good
            else if ( x .GE. 0.0012 .and. x .le. 0.5) then
                MagInt_Scalar_PiPi = (x**(2.9_dl))*(ssp2p1*lx9 + ssp2p2*lx8 + ssp2p3*lx7 + ssp2p4*lx6 +&
                ssp2p5*lx5 + ssp2p6*lx4 + ssp2p7*lx3 + ssp2p8*lx2 + ssp2p9*lx + ssp2p10)
                !write(*,*) "non-he,2, MagInt_Scalar_PiPi  = ", x, MagInt_Scalar_PiPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.5) then
                MagInt_Scalar_PiPi = (x**(2._dl*y+6._dl))*((68*x)/105 + 28/(5*x) + x2/16 - 5)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.0018) then
                MagInt_Scalar_PiPi = ssspp1*x6 + ssspp2*x5 + ssspp3*x4 + ssspp4*x3
                ! write(*,*) "non-he,t, MagInt_Scalar_PiPi  = ", x, MagInt_Scalar_PiPi !good,
            else if ( x .GE. 0.0018 .and. x .le. 0.5) then
                MagInt_Scalar_PiPi = sssp2p1*x9 + sssp2p2*x8 + sssp2p3*x7 + sssp2p4*x6 +&
                sssp2p5*x5 + sssp2p6*x4 + sssp2p7*x3 + sssp2p8*x2 + sssp2p9*x + sssp2p10
                !write(*,*) "non-he,2, MagInt_Scalar_PiPi  = ", x, MagInt_Scalar_PiPi   !  good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Scalar_PiPi = (x**(2._dl*y+6._dl))*(4/(105*x) + 1/(16*x2)&
                + 16/(105*x3) - 1/x4 + 28/(25*x5) - 1/25)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Scalar_PiPi = (x**(2._dl*y+6._dl))*(8/(15*x5) - 1/(8*x4) &
                - 1/(240*x2) - 1/x6 + 4/(5*x7) + 1/640)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Scalar_PiPi = (x**(2._dl*y+6._dl))*(4/(3465*x3) + 52/(525*x5)&
                - 7/(16*x6) + 628/(735*x7) - 1/x8 + 28/(45*x9) - 1/1225)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Scalar_PiPi

    !---------------------------------------------------------------------------
    !> this functions computes the scalar Pi-Pi integral for non-helical modes
    function MagInt_Scalar_Helical_PiPi(x,y)

        implicit none

        real(dl) :: MagInt_Scalar_Helical_PiPi, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        !n=  T VS TP, x=t, y=tp , N=-1.5
        !!!!!t<  0.0012
        !Linear model Poly9:
        !   f(x) = hsp1*x9 + hsp2*x8 + hsp3*x7 + hsp4*x6 +
        !                 hsp5*x5 + hsp6*x4 + hsp7*x3 + hsp8*x2 + hsp9*x + hsp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hsp1 =   1.172e+18  !(1.093e+18, 1.251e+18)
        real(dl) :: hsp2 =  -7.041e+15  !(-7.461e+15, -6.621e+15)
        real(dl) :: hsp3 =   1.826e+13  !(1.733e+13, 1.92e+13)
        real(dl) :: hsp4 =  -2.708e+10  !(-2.821e+10, -2.594e+10)
        real(dl) :: hsp5 =    2.61e+07  !(2.529e+07, 2.691e+07)
        real(dl) :: hsp6 =  -1.947e+04  !(-1.981e+04, -1.912e+04)
        real(dl) :: hsp7 =       37.79  !(37.7, 37.87)
        real(dl) :: hsp8 =   0.0003971  !(0.0003865, 0.0004076)
        real(dl) :: hsp9 =   -6.98e-09  !(-7.523e-09, -6.437e-09)
        real(dl) :: hsp10 =   2.637e-14  !(2.297e-14, 2.977e-14)

        !Goodness of fit:
        ! SSE: 4.88e-26
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.297e-14
        !!!!t>0.0012
        !Linear model Poly9:
        !    f(x) = hs2p1*x9 + hs2p2*x8 + hs2p3*x7 + hs2p4*x6 +
        !                  hs2p5*x5 + hs2p6*x4 + hs2p7*x3 + hs2p8*x2 + hs2p9*x + hs2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hs2p1 =   3.867e+07  !(3.821e+07, 3.913e+07)
        real(dl) :: hs2p2 =   -1.22e+07  !(-1.233e+07, -1.208e+07)
        real(dl) :: hs2p3 =   1.684e+06  !(1.669e+06, 1.699e+06)
        real(dl) :: hs2p4 =  -1.352e+05  !(-1.361e+05, -1.343e+05)
        real(dl) :: hs2p5 =        7245  !(7209, 7281)
        real(dl) :: hs2p6 =      -313.6  !(-314.5, -312.8)
        real(dl) :: hs2p7 =       21.06  !(21.05, 21.07)
        real(dl) :: hs2p8 =     0.02821  !(0.02812, 0.02831)
        real(dl) :: hs2p9 =   -4.77e-05  !(-4.805e-05, -4.735e-05)
        real(dl) :: hs2p10 =   3.702e-08  !(3.656e-08, 3.749e-08)

        !Goodness of fit:
        !! SSE: 2.985e-15
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 7.13e-10

        !!!!!!!!!!!!! n=-1.005
        !!! t<0.00595
        ! t vs tp
        !General model:
        !  f(x) = hhsp3*x7 + hhsp4*x6 + hhsp5*x5 + hhsp6*x4 + hhsp7*x3 + hhsp8*x**(3.8_dl)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hhsp3 =      0.7147  !(0.3985, 1.031)
        real(dl) :: hhsp4 =     -0.4131  !(-0.4933, -0.3328)
        real(dl) :: hhsp5 =       0.387  !(0.3753, 0.3986)
        real(dl) :: hhsp6 =      -2.967  !(-2.97, -2.965)
        real(dl) :: hhsp7 =        4.04  !(4.04, 4.04)
        real(dl) :: hhsp8 =    -0.06261  !(-0.06381, -0.06142)

        !Goodness of fit:
        ! SSE: 2.582e-18
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 2.047e-11

        !!!! t>0.00595
        !Linear model Poly9:
        !   f(x) = hhs2p1*x9 + hhs2p2*x8 + hhs2p3*x7 + hhs2p4*x6 +
        !                 hhs2p5*x5 + hhs2p6*x4 + hhs2p7*x3 + hhs2p8*x2 + hhs2p9*x + hhs2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hhs2p1 =   1.665e+04  !(-1.439e+04, 4.769e+04)
        real(dl) :: hhs2p2 =       -5016  !(-1.423e+04, 4199)
        real(dl) :: hhs2p3 =       653.1  !(-517.3, 1823)
        real(dl) :: hhs2p4 =      -49.03  !(-132.1, 34.07)
        real(dl) :: hhs2p5 =       2.808  !(-0.8062, 6.422)
        real(dl) :: hhs2p6 =      -3.138  !(-3.237, -3.039)
        real(dl) :: hhs2p7 =       4.041  !(4.039, 4.042)
        real(dl) :: hhs2p8 =  -4.421e-06  !(-2.209e-05, 1.325e-05)
        real(dl) :: hhs2p9 =   2.653e-08  !(-7.216e-08, 1.252e-07)
        real(dl) :: hhs2p10 =  -5.815e-11  !(-2.834e-10, 1.671e-10)

        !Goodness of fit:
        ! SSE: 2.579e-18
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 2.185e-11
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = (x**(2._dl*y+6._dl))*(8 - x2/2 - 4*x)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.0012 ) then
                MagInt_Scalar_Helical_PiPi = hsp1*x9 + hsp2*x8 + hsp3*x7 + hsp4*x6 + hsp5*x5&
                + hsp6*x4 + hsp7*x3 + hsp8*x2 + hsp9*x + hsp10
                !write(*,*) "he,t, MagInt_Scalar_Helical_PiPi  = ", x, MagInt_Scalar_Helical_PiPi  !good
            else if ( x .GE. 0.0012  .and. x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = hs2p1*x9 + hs2p2*x8 + hs2p3*x7 + hs2p4*x6 + hs2p5*x5&
                + hs2p6*x4 + hs2p7*x3 + hs2p8*x2 + hs2p9*x + hs2p10
                !write(*,*) "he,2, MagInt_Scalar_Helical_PiPi  = ", x, MagInt_Scalar_Helical_PiPi  !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.00595) then
                MagInt_Scalar_Helical_PiPi = (hhsp3*x7 + hhsp4*x6 + hhsp5*x5 + hhsp6*x4 + hhsp7*x3 + hhsp8*x**(3.8_dl))
                ! hhsp1*x9 + hhsp2*x8 + hhsp3*x7 + hhsp4*x6 + hhsp5*x5&
                !          + hhsp6*x4 + hhsp7*x3 + hhsp8*x2 + hhsp9*x + hhsp10
                !write(*,*) "he,t, MagInt_Scalar_Helical_PiPi  = ", x, MagInt_Scalar_Helical_PiPi  !good
            else if ( x .GE. 0.00595 .and. x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = hhs2p1*x9 + hhs2p2*x8 + hhs2p3*x7 + hhs2p4*x6 + hhs2p5*x5&
                + hhs2p6*x4 + hhs2p7*x3 + hhs2p8*x2 + hhs2p9*x + hhs2p10
                !write(*,*) "he,2, MagInt_Scalar_Helical_PiPi  = ", x, MagInt_Scalar_Helical_PiPi  !g
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = (x**(2._dl*y+6._dl))*(4/(3*x3) - 1/(2*x2) &
                - 4/(15*x) + 3._dl/36028797018963968._dl)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = (x**(2._dl*y+6._dl))*(4/(5*x5) - 1/(2*x4) + 1/160)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = (x**(2._dl*y+6._dl))*(4/(315*x3) + 8/(75*x5)&
                - 1/(2*x6) + 4/(7*x7) - 2/525)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_PiPi = (x**(2._dl*y+6._dl))*(4/(21*x7) - 1/(96*x4) &
                - 1/(2*x8) + 4/(9*x9) + 1/4032 )
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Scalar_Helical_PiPi


    !---------------------------------------------------------------------------
    !> this functions computes the scalar Delta-Pi integral for non-helical modes
    function MagInt_Scalar_DeltaPi(x,y)

        implicit none

        real(dl) :: MagInt_Scalar_DeltaPi, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        !n=  ! T VS P, x=t, y=p , N=-2.05

        !f(x) = sdp1*x5 + sdp2*x4 + sdp3*x3 + sdp4*x2 + sdp5*x + sdp6

        ! f(x) = sdpa*exp(sdpb*x) + sdpc*exp(sdpd*x)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sdpa =   0.0003878 ! (-0.001127, 0.001902)
        real(dl) :: sdpb =       18.99  !(-9.348, 47.34)
        real(dl) :: sdpc =      -3.963  !(-3.965, -3.962)
        real(dl) :: sdpd =  -0.0001521  !(-0.005138, 0.004834)

        !Goodness of fit:
        ! SSE: 6.322e-05
        !R-square: 0.7511
        !Adjusted R-square: 0.7509
        !RMSE: 0.0001434

        !!!!!!!!!!! n=-1.55, T VS P

        !f(x) = ssdpa*exp(ssdpb*x) + ssdpc*exp(ssdpd*x)
        !!!! t vs p ,  t<1.000200E-05
        !General model Exp2:
        !    f(x) = a*exp(b*x) + c*exp(d*x)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: ssdpa =     -0.1279 ! (-1.021, 0.765)
        real(dl) :: ssdpb =  -4.591e+05  !(-6.13e+06, 5.211e+06)
        real(dl) :: ssdpc =      -2.411  !(-3.139, -1.683)
        real(dl) :: ssdpd =       -6194  !(-3.326e+04, 2.087e+04)

        !Goodness of fit:
        ! SSE: 2.678
        !R-square: 0.0467
        !Adjusted R-square: -0.04267
        !RMSE: 0.2893

        !t>0.00001
        !Linear model Poly9:
        !    f(x) = ssd2p1*x9 + ssd2p2*x8 + ssd2p3*x7 + ssd2p4*x6 +
        !                 ssd2p5*x5 + ssd2p6*x4 + ssd2p7*x3 + ssd2p8*x2 + ssd2p9*x + ssd2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: ssd2p1 =  -8.945e+04  !(-1.28e+05, -5.088e+04)
        real(dl) :: ssd2p2 =   1.259e+05  !(7.527e+04, 1.765e+05)
        real(dl) :: ssd2p3 =  -7.451e+04  !(-1.022e+05, -4.687e+04)
        real(dl) :: ssd2p4 =   2.413e+04  !(1.599e+04, 3.228e+04)
        real(dl) :: ssd2p5 =       -4662  !(-6064, -3259)
        real(dl) :: ssd2p6 =       551.5  !(408.5, 694.4)
        real(dl) :: ssd2p7 =      -39.67  !(-48.01, -31.32)
        real(dl) :: ssd2p8 =        2.41  !(2.154, 2.666)
        real(dl) :: ssd2p9 =      0.1409  !(0.1374, 0.1444)
        real(dl) :: ssd2p10 =      -2.305  !(-2.305, -2.305)

        !Goodness of fit:
        ! SSE: 2.675e-05
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 8.767e-05

        !!!!!!!  N=-0.005, t vs tp
        !   f(x) = sdp1*x8 + sdp2*x7 + sdp3*x6 + sdp4*x5 +
        !                 sdp5*x4 + sdp6*x3 + sdp7*x2 + sdp8*x + sdp9
        !Coefficients (with 95% confidence bounds):
        real(dl) :: sdp1 =      0.8822 ! (-0.7014, 2.466)
        real(dl) :: sdp2 =     -0.3097  !(-0.6919, 0.07245)
        real(dl) :: sdp3 =      0.7403  !(0.7024, 0.7782)
        real(dl) :: sdp4 =      -1.077  !(-1.078, -1.075)
        real(dl) :: sdp5 =        0.25  !(0.25, 0.2501)
        real(dl) :: sdp6 =  -2.456e-07  !(-1.25e-06, 7.592e-07)
        real(dl) :: sdp7 =   1.947e-09  !(-7.124e-09, 1.102e-08)
        real(dl) :: sdp8 =  -6.578e-12  !(-4.415e-11, 3.1e-11)
        real(dl) :: sdp9 =   2.093e-15  !(-4.886e-14, 5.305e-14)

        !Goodness of fit:
        ! SSE: 5.541e-23
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.364e-13

        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.5) then !n=-2.05, 1.9_dl
                MagInt_Scalar_DeltaPi = (x**(1.9_dl))*(sdpa*exp(sdpb*x) + sdpc*exp(sdpd*x))
                !MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*(sdpa*exp(sdpb*x) + sdpc*exp(sdpd*x))
                !write(*,*) "non-he,t, MagInt_Scalar_DeltaPi  = ", x, MagInt_Scalar_DeltaPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.00001) then!n=-1.55, 2.9_dl， 2._dl*y+6._dl
                MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*(ssdpa*exp(ssdpb*x) + ssdpc*exp(ssdpd*x))
                !write(*,*) "non-he,t, MagInt_Scalar_DeltaPi  = ", x, MagInt_Scalar_DeltaPi   ! good
            else if ( x .GE. 0.00001 .and. x .le. 0.5) then
                MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*(ssd2p1*x9 + ssd2p2*x8 + ssd2p3*x7 + ssd2p4*x6 +&
                ssd2p5*x5 + ssd2p6*x4 + ssd2p7*x3 + ssd2p8*x2 + ssd2p9*x + ssd2p10)
                !write(*,*) "non-he,2, MagInt_Scalar_DeltaPi  = ", x, MagInt_Scalar_DeltaPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*((16*x)/15 + x2/8 - 7/4)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaPi = sdp1*x8 + sdp2*x7 + sdp3*x6 + sdp4*x5 + sdp5*x4 + sdp6*x3 + sdp7*x2 + sdp8*x + sdp9
                !write(*,*) "non-he,t, MagInt_Scalar_DeltaPi  = ", x, MagInt_Scalar_DeltaPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*(8/(105*x) + 1/(8*x2) &
                - 8/(15*x3) + 1/(4*x4) - 1._dl/144115188075855872._dl)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*(1/(4*x4) - 1/(192*x2) &
                - 8/(15*x5) + 1/(4*x6) - 1/640)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Scalar_DeltaPi = (x**(2._dl*y+6._dl))*(11/(24*x6) - 64/(525*x5) &
                - 64/(105*x7) + 1/(4*x8) + 1/1050)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Scalar_DeltaPi

    !---------------------------------------------------------------------------
    !> this functions computes the scalar Delta-Pi integral for helical modes
    function IntSDPH(x,y)

        implicit none

        real(dl) :: MagInt_Scalar_Helical_DeltaPi, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        !n=  , x=t, y=tp, N=-1.5  ! T VS P

        !   f(x) = hdpa*exp(hdpb*x) + hdpc*exp(hdpd*x)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hdpa =  -0.0001159 ! (-0.06067, 0.06044)
        real(dl) :: hdpb =       17.22  !(-3462, 3497)
        real(dl) :: hdpc =       1.334  !(1.273, 1.394)
        real(dl) :: hdpd =      0.1797  !(-0.3362, 0.6956)

        !Goodness of fit:
        ! SSE: 0.1712
        !R-square: 0.3726
        ! Adjusted R-square: 0.3723
        !RMSE: 0.005343

        !!!!!!!!!!!!! n=-1.005
        !!!T<0.004

        !t vs tp General model:

        !General model:
        !  f(x) =hhdp4*x6 + hhdp5*x5 + hhdp6*x4 + hhdp7*x3 + hhdp8*x**(3.8_dl)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hhdp4 =     0.01238 ! (0.009807, 0.01495)
        real(dl) :: hhdp5 =     -0.5847  !(-0.5853, -0.5842)
        real(dl) :: hhdp6 =       1.466  !(1.466, 1.466)
        real(dl) :: hhdp7 =  -3.145e-05  !(-3.189e-05, -3.101e-05)
        real(dl) :: hhdp8 =     0.03053  !(0.03046, 0.0306)

        !p2*x^(3.9)-p3*x^(5)+p4*x^(6.7)+p1*x^(4.7)-p5*x^4-p6*x^(4.1)+p7*x^(3.7)
        !Goodness of fit:
        ! SSE: 8.871e-21
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.217e-12

        !!!! T>0.004
        !Linear model Poly9:
        !  f(x) = hhd2p1*x9 + hhd2p2*x8 + hhd2p3*x7 + hhd2p4*x6 +&
        !                hhd2p5*x5 + hhd2p6*x4 + hhd2p7*x3 + hhd2p8*x2 + hhd2p9*x + hhd2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hhd2p1 =       -1276 ! (-2444, -106.9)
        real(dl) :: hhd2p2 =       434.1  !(96.83, 771.4)
        real(dl) :: hhd2p3 =      -67.34  !(-108.8, -25.9)
        real(dl) :: hhd2p4 =       6.404  !(3.577, 9.232)
        real(dl) :: hhd2p5 =       -1.04  !(-1.157, -0.9224)
        real(dl) :: hhd2p6 =        1.53  !(1.527, 1.533)
        real(dl) :: hhd2p7 =   0.0001978  !(0.0001492, 0.0002464)
        real(dl) :: hhd2p8 =  -8.883e-07  !(-1.347e-06, -4.299e-07)
        real(dl) :: hhd2p9 =   2.822e-09  !(5.442e-10, 5.1e-09)
        real(dl) :: hhd2p10 =   -3.21e-12 !(-7.693e-12, 1.272e-12)

        !Goodness of fit:
        ! SSE: 7.137e-21
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.13e-12
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = (x**(2._dl*y+6._dl))*(x2/4 + 2)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = (x**(2._dl*y+6._dl))*(hdpa*exp(hdpb*x) + hdpc*exp(hdpd*x))
                !MagInt_Scalar_Helical_DeltaPi = hdp1*x5 + hdp2*x4 + hdp3*x3 + hdp4*x2 + hdp5*x + hdp6
                ! write(*,*) "he,t, MagInt_Scalar_Helical_DeltaPi  = ", x, MagInt_Scalar_Helical_DeltaPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.004) then
                MagInt_Scalar_Helical_DeltaPi = hhdp4*x6 + hhdp5*x5 + hhdp6*x4 + hhdp7*x3 + hhdp8*x**(3.8_dl)
                !hhdp5*x5 + hhdp6*x4 + hhdp7*x3
                ! hhdp1*x9 + hhdp2*x8 + hhdp3*x7 + hhdp4*x6 +&
                !       hhdp5*x5 + hhdp6*x4 + hhdp7*x3 + hhdp8*x2 + hhdp9*x + hhdp10
                ! write(*,*) "he,t, MagInt_Scalar_Helical_DeltaPi  = ", x, MagInt_Scalar_Helical_DeltaPi   ! good,  check back negetive
            else if ( x .GE. 0.004 .and. x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = hhd2p1*x9 + hhd2p2*x8 + hhd2p3*x7 + hhd2p4*x6 +&
                hhd2p5*x5 + hhd2p6*x4 + hhd2p7*x3 + hhd2p8*x2 + hhd2p9*x + hhd2p10
            !write(*,*) "he,2, MagInt_Scalar_Helical_DeltaPi  = ", x, MagInt_Scalar_Helical_DeltaPi   ! good,
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = (x**(2._dl*y+6._dl))*(8/(15*x) + 1/(4*x2) - 1/2)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = (x**(2._dl*y+6._dl))*(1/(4*x4) - 1/(8*x2) + 1/64)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = (x**(2._dl*y+6._dl))*(16/(315*x3) &
                - 16/(75*x5) + 1/(4*x6) - 1/150)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Scalar_Helical_DeltaPi = (x**(2._dl*y+6._dl))*(5/(24*x6) - 5/(192*x4) &
                - 8/(21*x7) + 1/(4*x8) + 1/2688)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Scalar_Helical_DeltaPi

    !---------------------------------------------------------------------------
    !> this functions computes the vector Pi-Pi integral for non-helical modes
    function MagInt_Vector_PiPi(x,y)
        implicit none
        !variables
        real(dl) :: MagInt_Vector_PiPi, x,y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9

        !  N=-2.005   !!! t vs p

        !  f(x) = va*exp(vb*x) + vc*exp(vd*x)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: va =       3.697 ! (3.696, 3.698)
        real(dl) :: vb =     -0.5098  !(-0.5133, -0.5063)
        real(dl) :: vc =    -0.07685  !(-0.08638, -0.06732)
        real(dl) :: vd =       -1469  !(-1666, -1273)

        !Goodness of fit:
        ! SSE: 0.07956
        ! R-square: 0.9778
        !Adjusted R-square: 0.9777
        !RMSE: 0.006323

        !!!!!!!!!!!!!!!!! n=-1.5
        !!!n=-1.5 t vs p, t<0.0008
        !Linear model Poly6:
        !   f(x) = vvp1*x6 + vvp2*x5 + vvp3*x4 + vvp4*x3 + vvp5*x2 +
        !                 vvp6*x + vvp7
        !Coefficients (with 95% confidence bounds):
        real(dl) :: vvp1 =   1.697e+10 ! (1.492e+10, 1.902e+10)
        real(dl) :: vvp2 =  -3.365e+09  !(-3.737e+09, -2.992e+09)
        real(dl) :: vvp3 =   2.607e+08  !(2.347e+08, 2.867e+08)
        real(dl) :: vvp4 =  -9.972e+06  !(-1.084e+07, -9.101e+06)
        real(dl) :: vvp5 =   1.967e+05  !(1.825e+05, 2.109e+05)
        real(dl) :: vvp6 =       -1969  !(-2071, -1868)
        real(dl) :: vvp7 =       16.97  !(16.74, 17.2)

        !Goodness of fit:
        ! SSE: 1.918e+04
        !R-square: 0.6375
        !Adjusted R-square: 0.6373
        !RMSE: 1.391

        !!!!! t>0.0008  t vs p
        !Linear model Poly9:
        !  f(x) = p1*x^9 + p2*x^8 + p3*x^7 + p4*x^6 +
        !                p5*x^5 + p6*x^4 + p7*x^3 + p8*x^2 + p9*x + p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: vv2p1 =  -6.373e+14 ! (-6.855e+14, -5.89e+14)
        real(dl) :: vv2p2 =   1.869e+14  !(1.737e+14, 2.002e+14)
        real(dl) :: vv2p3 =  -2.334e+13  !(-2.486e+13, -2.182e+13)
        real(dl) :: vv2p4 =   1.618e+12  !(1.522e+12, 1.714e+12)
        real(dl) :: vv2p5 =  -6.817e+10  !(-7.179e+10, -6.456e+10)
        real(dl) :: vv2p6 =   1.799e+09  !(1.717e+09, 1.882e+09)
        real(dl) :: vv2p7 =  -2.965e+07  !(-3.078e+07, -2.853e+07)
        real(dl) :: vv2p8 =   2.986e+05  !(2.901e+05, 3.07e+05)
        real(dl) :: vv2p9 =       -1862  !(-1892, -1831)
        real(dl) :: vv2p10 =       15.46  !(15.42, 15.5)

        !Goodness of fit:
        ! SSE: 106.3
        ! R-square: 0.996
        !Adjusted R-square: 0.996
        !RMSE: 0.1038


        !!!!!!!  N=-0.005
        !!! t< 0.0008
        ! t vs tp
        !General model:
        !   f(x) = vvv1p1*x6 + vvv1p2*x5 + vvv1p3*x4 + vvv1p4*x3
        !Coefficients (with 95% confidence bounds):
        real(dl) :: vvv1p1 =     0.03989  !(0.02159, 0.05819)
        real(dl) :: vvv1p2 =     -0.1563  !(-0.1588, -0.1539)
        real(dl) :: vvv1p3 =     -0.4166  !(-0.4167, -0.4165)
        real(dl) :: vvv1p4 =      0.6243  !(0.6243, 0.6243)

        !Goodness of fit:
        ! SSE: 2.969e-17
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 5.493e-11

        !!!!t>0.0008
        !   f(x) = vvvp1*x9 + vvvp2*x8 + vvvp3*x7 + vvvp4*x6 +
        !             vvvp5*x5 + vvvp6*x4 + vvvp7*x3 + vvvp8*x2 + vvvp9*x + vvvp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: vvvp1 =    3.78e+04  !(2.419e+04, 5.141e+04)
        real(dl) :: vvvp2 =  -1.131e+04  !(-1.511e+04, -7501)
        real(dl) :: vvvp3 =        1442  !(991.9, 1891)
        real(dl) :: vvvp4 =      -102.2  !(-131.5, -72.92)
        real(dl) :: vvvp5 =       4.256  !(3.111, 5.402)
        real(dl) :: vvvp6 =     -0.5358  !(-0.5633, -0.5082)
        real(dl) :: vvvp7 =      0.6263  !(0.6259, 0.6267)
        real(dl) :: vvvp8 =  -1.971e-05  !(-2.305e-05, -1.638e-05)
        real(dl) :: vvvp9 =   1.038e-07  !(8.992e-08, 1.178e-07)
        real(dl) :: vvvp10 =  -2.233e-10  !(-2.449e-10, -2.018e-10)

        !Goodness of fit:
        ! SSE: 5.498e-18
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 2.387e-11
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        !!!!!!!!!!!!!!!!!!!!!
        if (y==-2) then
            if (x .le. 0.5) then !n=-2.005
                MagInt_Vector_PiPi = (x**(1.99_dl))*(va*exp(vb*x) + vc*exp(vd*x))
                !write(*,*) "non-he,t, MagInt_Vector_PiPi  = ", x, MagInt_Vector_PiPi   ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.0008) then
                MagInt_Vector_PiPi = (x**(2._dl*y+6._dl))*(vvp1*x6 + vvp2*x5 + vvp3*x4 +&
                vvp4*x3 + vvp5*x2 + vvp6*x + vvp7)
                !write(*,*) "non-he,t, MagInt_Vector_PiPi  = ", x, MagInt_Vector_PiPi   !good
            else if ( x .GE. 0.0008 .and. x .le. 0.5) then
                MagInt_Vector_PiPi = (x**(2._dl*y+6._dl))*(vv2p1*x9 + vv2p2*x8 + vv2p3*x7 + vv2p4*x6 + &
                vv2p5*x5 + vv2p6*x4 + vv2p7*x3 + vv2p8*x2 + vv2p9*x + vv2p10)
            !write(*,*) "non-he,2, MagInt_Vector_PiPi  = ", x, MagInt_Vector_PiPi   !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.5) then
                MagInt_Vector_PiPi = (x**(2._dl*y+6._dl))*((16*x)/105 + 28/(15*x) - 7/4)
                !write(*,*) "non-he,t, MagInt_Vector_PiPi  = ", x, MagInt_Vector_PiPi
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.000756) then   !-0.005000,  5.99_dl
                MagInt_Vector_PiPi = vvv1p1*x6 + vvv1p2*x5 + vvv1p3*x4 + vvv1p4*x3
                !write(*,*) "non-he,t, MagInt_Vector_PiPi  = ", x, MagInt_Vector_PiPi  !good ,check
            else if ( x .GE. 0.000756 .and. x .le. 0.5) then
                MagInt_Vector_PiPi = vvvp1*x9 + vvvp2*x8 + vvvp3*x7 + vvvp4*x6 +&
                vvvp5*x5 + vvvp6*x4 + vvvp7*x3 + vvvp8*x2 + vvvp9*x + vvvp10
                !write(*,*) "non-he,2, MagInt_Vector_PiPi  = ", x, MagInt_Vector_PiPi  !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Vector_PiPi = (x**(2._dl*y+6._dl))*(4/(35*x3) - 8/(315*x) &
                - 5/(12*x4) + 28/(75*x5) + 1/50)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Vector_PiPi = (x**(2._dl*y+6._dl))*(7/(960*x2) - 1/(12*x4) &
                + 4/(15*x5) - 5/(12*x6) + 4/(15*x7) - 1/1920)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Vector_PiPi = (x**(2._dl*y+6._dl))*(92/(1575*x5) - 32/(10395*x3)&
                - 2/(9*x6) + 296/(735*x7) - 5/(12*x8) + 28/(135*x9) + 2/11025)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Vector_PiPi

    !---------------------------------------------------------------------------
    !> this functions computes the vector Pi-Pi integral for non-helical modes
    function MagInt_Vector_Helical_PiPi(x,y)
        implicit none
        !variables
        real(dl) :: MagInt_Vector_Helical_PiPi, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        !t vs tp   ! N=-1.5
        !   f(x) = hp1*x9 + hp2*x8 + hp3*x7 + hp4*x6 +
        !                 hp5*x5 + hp6*x4 + hp7*x3 + hp8*x2 + hp9*x + hp10
        !Coefficients (with 95% confidence bounds):

        !!!!!!!!!!! t<0.00041
        !Linear model Poly9:
        !   f(x) = p1*x^9 + p2*x^8 + p3*x^7 + p4*x^6 +
        !                 p5*x^5 + p6*x^4 + p7*x^3 + p8*x^2 + p9*x + p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: hp1 =    2.41e+20 ! (2.239e+20, 2.582e+20)
        real(dl) :: hp2 =  -4.997e+17  !(-5.309e+17, -4.685e+17)
        real(dl) :: hp3 =   4.468e+14  !(4.232e+14, 4.705e+14)
        real(dl) :: hp4 =  -2.279e+11  !(-2.376e+11, -2.183e+11)
        real(dl) :: hp5 =   7.538e+07  !(7.307e+07, 7.769e+07)
        real(dl) :: hp6 =  -1.919e+04  !(-1.952e+04, -1.887e+04)
        real(dl) :: hp7 =       13.61  !(13.58, 13.63)
        real(dl) :: hp8 =   4.383e-05  !(4.284e-05, 4.482e-05)
        real(dl) :: hp9 =  -2.437e-10  !(-2.596e-10, -2.277e-10)
        real(dl) :: hp10 =   4.308e-16  !(3.654e-16, 4.962e-16)

        !Goodness of fit:
        ! SSE: 6.75e-30
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.677e-16

        !!!!!  0.00041<t<0.004
        !Linear model Poly9:
        !   f(x) = h1p1*x9 + h1p2*x8 + h1p3*x7 + h1p4*x6 +
        !                 h1p5*x5 + h1p6*x4 + h1p7*x3 + h1p8*x2 + h1p9*x + h1p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: h1p1 =   4.669e+13  !(4.411e+13, 4.927e+13)
        real(dl) :: h1p2 =  -1.107e+12  !(-1.16e+12, -1.054e+12)
        real(dl) :: h1p3 =   1.169e+10  !(1.122e+10, 1.215e+10)
        real(dl) :: h1p4 =   -7.35e+07  !(-7.58e+07, -7.121e+07)
        real(dl) :: h1p5 =   3.186e+05  !(3.117e+05, 3.255e+05)
        real(dl) :: h1p6 =       -1163  !(-1177, -1150)
        real(dl) :: h1p7 =        9.74  !(9.724, 9.756)
        real(dl) :: h1p8 =   0.0009195  !(0.0009082, 0.0009308)
        real(dl) :: h1p9 =  -1.721e-07  !(-1.764e-07, -1.677e-07)
        real(dl) :: h1p10 =   1.766e-11  !(1.697e-11, 1.836e-11)

        !Goodness of fit:
        ! SSE: 3.071e-25
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 2.237e-14

        !!!!!t>0.004

        !Linear model Poly9:
        !  f(x) = h2p1*x9 + h2p2*x8 + h2p3*x7 + h2p4*x6 +
        !                 h2p5*x5 + h2p6*x4 + h2p7*x3 + h2p8*x2 + h2p9*x + h2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: h2p1 =   6.728e+06  !(6.676e+06, 6.78e+06)
        real(dl) :: h2p2 =  -2.252e+06  !(-2.267e+06, -2.237e+06)
        real(dl) :: h2p3 =   3.336e+05  !(3.318e+05, 3.355e+05)
        real(dl) :: h2p4 =  -2.923e+04  !(-2.935e+04, -2.91e+04)
        real(dl) :: h2p5 =        1747  !(1742, 1752)
        real(dl) :: h2p6 =      -87.09  !(-87.23, -86.96)
        real(dl) :: h2p7 =       6.293  !(6.291, 6.295)
        real(dl) :: h2p8 =     0.01206  !(0.01204, 0.01208)
        real(dl) :: h2p9 =  -2.885e-05  !(-2.896e-05, -2.875e-05)
        real(dl) :: h2p10 =   3.672e-08  !(3.651e-08, 3.692e-08)

        !Goodness of fit:
        ! SSE: 3.807e-17
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 6.397e-11

        !!!!n=-1.005, new
        !  f(x) = hvp1*x5 + hvp2*x4 + hvp3*x3 + hvp4*x2 + hvp5*x + hvp6
        !!!!T<0.004

        !    f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
        !Coefficients (with 95% confidence bounds):
        real(dl) ::  hvp1 =       2.204 ! (1.633, 2.774)
        real(dl) ::  hvp2 =      -1.566  !(-1.572, -1.56)
        real(dl) ::  hvp3 =       1.347  !(1.347, 1.347)
        real(dl) ::  hvp4 =  -9.115e-09  !(-4.086e-08, 2.263e-08)
        real(dl) ::  hvp5 =   1.824e-11  !(-1.017e-12, 3.75e-11)
        real(dl) ::  hvp6 =   2.628e-13  !(2.594e-13, 2.662e-13)

        !Goodness of fit:
        ! SSE: 8.458e-26
        ! R-square: 1
        ! Adjusted R-square: 1
        ! RMSE: 1.103e-14

        !!!!T>0.004
        !   f(x) = hv2p1*x6 + hv2p2*x5 + hv2p3*x4 + hv2p4*x3 + hv2p5*x2 +
        !                 hv2p6*x + hv2p7
        !Coefficients (with 95% confidence bounds):
        real(dl) ::  hv2p1 =      -0.447 ! (-0.4828, -0.4112)
        real(dl) ::  hv2p2 =       0.442  !(0.4351, 0.4489)
        real(dl) ::  hv2p3 =      -1.531  !(-1.532, -1.531)
        real(dl) ::  hv2p4 =       1.346  !(1.346, 1.346)
        real(dl) ::  hv2p5 =   3.909e-06  !(3.546e-06, 4.273e-06)
        real(dl) ::  hv2p6 =  -2.131e-08  !(-2.456e-08, -1.807e-08)
        real(dl) ::  hv2p7 =   4.679e-11  !(3.645e-11, 5.713e-11)

        !Goodness of fit:
        ! SSE: 2.445e-18
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.62e-11
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8

        if (y==-2) then
            if (x .le. 0.5) then
                MagInt_Vector_Helical_PiPi = (x**(2._dl*y+6._dl))*(2 - x2/4 - (4*x)/3)
                !write(*,*) "t, IntVPP2  = ", x, MagInt_Vector_Helical_PiPi
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.00041) then
                MagInt_Vector_Helical_PiPi =  hp1*x9 + hp2*x8 + hp3*x7 + hp4*x6 + &
                hp5*x5 + hp6*x4 + hp7*x3 + hp8*x2 + hp9*x + hp10
                !write(*,*) "t, MagInt_Vector_Helical_PiPi  = ", x, MagInt_Vector_Helical_PiPi   !good
            else if ( x .GE. 0.00041 .and. x .le. 0.004) then
                MagInt_Vector_Helical_PiPi = h1p1*x9 + h1p2*x8 + h1p3*x7 + h1p4*x6 +&
                h1p5*x5 + h1p6*x4 + h1p7*x3 + h1p8*x2 + h1p9*x + h1p10
                !write(*,*) "1, MagInt_Vector_Helical_PiPi  = ", x, MagInt_Vector_Helical_PiPi !good
            else if ( x .GE. 0.004 .and. x .le. 0.5) then
                MagInt_Vector_Helical_PiPi = h2p1*x9 + h2p2*x8 + h2p3*x7 + h2p4*x6 +&
                h2p5*x5 + h2p6*x4 + h2p7*x3 + h2p8*x2 + h2p9*x + h2p10
                !write(*,*) "2, MagInt_Vector_Helical_PiPi  = ", x, MagInt_Vector_Helical_PiPi   !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if ( y < -0.95 .and. y .GE. -1.1) then !-1.05
            if (x .le. 0.004) then
                MagInt_Vector_Helical_PiPi = hvp1*x5 + hvp2*x4 + hvp3*x3 + hvp4*x2 + hvp5*x + hvp6
                !write(*,*) "t, MagInt_Vector_Helical_PiPi  = ", x, MagInt_Vector_Helical_PiPi   !good
            else if ( x .GE. 0.004 .and. x .le. 0.5) then
                MagInt_Vector_Helical_PiPi = hv2p1*x6 + hv2p2*x5 + hv2p3*x4 + hv2p4*x3 + hv2p5*x2 + &
                hv2p6*x + hv2p7
                !write(*,*) "t, MagInt_Vector_Helical_PiPi  = ", x, MagInt_Vector_Helical_PiPi   !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.5) then
                MagInt_Vector_Helical_PiPi =(x**(2._dl*y+6._dl))*(4/(9*x3) - 1/(4*x2) - 4/(15*x) + 1/6)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Vector_Helical_PiPi = (x**(2._dl*y+6._dl))*(1/(24*x2) - 1/(4*x4) + 4/(15*x5) - 1/320)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Vector_Helical_PiPi = (x**(2._dl*y+6._dl))*(8/(75*x5) - 4/(315*x3) - 1/(4*x6) + 4/(21*x7) + 1/1050)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Vector_Helical_PiPi = (x**(2._dl*y+6._dl))*(1/(192*x4) - 5/(72*x6) + 4/(21*x7) - 1/(4*x8) + 4/(27*x9) - 1/24192)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Vector_Helical_PiPi

    !---------------------------------------------------------------------------
    !> this functions computes the tensor Pi-Pi integral for helical modes
    function MagInt_Tensor_PiPi(x,y)

        implicit none

        real(dl) :: MagInt_Tensor_PiPi, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9
        real(dl) :: tw2, tw3, tw4, tw5, tw6, tw7, tw8, tw9
        real(dl) :: lx,lx2, lx3, lx4, lx5, lx6, lx7, lx8, lx9

        !!!!! t vs p, n=-2.05
        !   f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
        !Coefficients (with 95% confidence bounds):
        real(dl) :: tp1 =   1.251e+04  !(-1.572e+05, 1.822e+05)
        real(dl) :: tp2 =       -2488  !(-2.811e+04, 2.313e+04)
        real(dl) :: tp3 =       205.1  !(-1193, 1603)
        real(dl) :: tp4 =      -11.01  !(-44.2, 22.17)
        real(dl) :: tp5 =      -1.016  !(-1.339, -0.693)
        real(dl) :: tp6 =       5.661  !(5.66, 5.662)

        !Goodness of fit:
        ! SSE: 0.05733
        !R-square: 0.965
        !Adjusted R-square: 0.9649
        !RMSE: 0.004379

        !!!!!!!!!!!!! n=-1.55,  !!!!!!!!!!!!!!!!!!!!!!T<0.0028  ! t vs p
        !General model Exp2:
        !    f(x) = a*exp(b*x) + c*exp(d*x)
        !Coefficients (with 95% confidence bounds):
        real(dl) :: tta =       10.79 ! (1.504, 20.07)
        real(dl) :: ttb =      -140.1  !(-371.9, 91.74)
        real(dl) :: ttc =       1.609  !(-7.802, 11.02)
        real(dl) :: ttd =         236  !(-595.3, 1067)

        !Goodness of fit:
        ! SSE: 0.7727
        !R-square: 0.9724
        !Adjusted R-square: 0.9717
        !RMSE: 0.07991

        !!! T> 0.0028 , logt vs p

        !Linear model Poly9:
        !   f(x) = p1*x^9 + p2*x^8 + p3*x^7 + p4*x^6 +
        !                 p5*x^5 + p6*x^4 + p7*x^3 + p8*x^2 + p9*x + p10
        !     tt2p1*lx9 + tt2p2*lx8 + tt2p3*lx7 + tt2p4*lx6 +
        !                 tt2p5*lx5 + tt2p6*lx4 + tt2p7*lx3 + tt2p8*lx2 + tt2p9*lx + tt2p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: tt2p1 =  -6.081e-07  !(-6.61e-07, -5.552e-07)
        real(dl) :: tt2p2 =   -2.27e-05  !(-2.443e-05, -2.097e-05)
        real(dl) :: tt2p3 =  -0.0003774  !(-0.0004021, -0.0003527)
        real(dl) :: tt2p4 =   -0.003692  !(-0.003893, -0.003491)
        real(dl) :: tt2p5 =    -0.02366  !(-0.02469, -0.02263)
        real(dl) :: tt2p6 =     -0.1047  !(-0.1081, -0.1013)
        real(dl) :: tt2p7 =     -0.3302  !(-0.3376, -0.3227)
        real(dl) :: tt2p8 =      -0.813  !(-0.8231, -0.8029)
        real(dl) :: tt2p9 =      -2.916  !(-2.923, -2.908)
        real(dl) :: tt2p10 =       1.304  !(1.301, 1.306)

        !Goodness of fit:
        ! SSE: 4.017e-08
        ! R-square: 1
        !Adjusted R-square: 1
        !RMSE: 3.543e-06

        !!!!!!!  N=-0.005
        !!!!!!!!!!!!!!T<0.004
        !    f(x) = tttp1*x9 + tttp2*x8 + tttp3*x7 + tttp4*x6 +
        !                 tttp5*x5 + tttp6*x4 + tttp7*x3 + tttp8*x2 + tttp9*x + tttp10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: tttp1 =   8.318e+15  !(-7.842e+15, 2.448e+16)
        real(dl) :: tttp2 =  -1.637e+14  !(-4.643e+14, 1.369e+14)
        real(dl) :: tttp3 =   1.372e+12  !(-9.813e+11, 3.725e+12)
        real(dl) :: tttp4 =  -6.388e+09  !(-1.646e+10, 3.684e+09)
        real(dl) :: tttp5 =    1.81e+07  !(-7.517e+06, 4.373e+07)
        real(dl) :: tttp6 =  -3.214e+04  !(-7.154e+04, 7262)
        real(dl) :: tttp7 =       35.89  !(0.3286, 71.45)
        real(dl) :: tttp8 =    -0.02228  !(-0.03956, -0.004992)
        real(dl) :: tttp9 =   6.341e-06  !(2.618e-06, 1.006e-05)
        real(dl) :: tttp10 =   2.888e-11  !(-2.505e-10, 3.082e-10)

        !Goodness of fit:
        ! SSE: 3.428e-18
        !R-square: 0.9999
        !Adjusted R-square: 0.9999
        !RMSE: 1.42e-10

        !! t vs tp,  t<0.0011
        !General model:
        !    f(x) = p1*x^3
        !Coefficients (with 95% confidence bounds):
        real(dl) :: ttt1p1 =      0.6948 ! (0.6887, 0.701)

        !Goodness of fit:
        ! SSE: 2.483e-12
        !R-square: 0.9709
        !Adjusted R-square: 0.9709
        !RMSE: 6.406e-08

        !!!!!!! T>0.004
        ! f(x) = p1*x^5 + p2*x^4 + p3*x^3 + p4*x^2 + p5*x + p6
        !Coefficients (with 95% confidence bounds):
        real(dl) :: ttt2p1 =     0.08104 ! (0.08055, 0.08153)
        real(dl) :: ttt2p2 =     -0.5865  !(-0.5866, -0.5864)
        real(dl) :: ttt2p3 =      0.6244  !(0.6244, 0.6244)
        real(dl) :: ttt2p4 =   -1.92e-06  !(-2.046e-06, -1.794e-06)
        real(dl) :: ttt2p5 =   1.511e-08  !(1.362e-08, 1.66e-08)
        real(dl) :: ttt2p6 =  -3.963e-11  !(-4.559e-11, -3.367e-11)
        !Goodness of fit:
        ! SSE: 2.152e-19
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 8.774e-12
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8
        lx = log(x)
        lx2 = lx*lx
        lx3 = lx*lx2
        lx4 = lx*lx3
        lx5 = lx*lx4
        lx6 = lx*lx5
        lx7 = lx*lx6
        lx8 = lx*lx7
        lx9 = lx*lx8

        if (y==-2) then
            if (x .le. 0.5) then  !y=-2.05, 5.66 > y=-2,  5.55, t^1.9> t^2 check
                !MagInt_Tensor_PiPi = (x**(1.9_dl))*(tp1*x5 + tp2*x4 + tp3*x3 + tp4*x2 + tp5*x + tp6)
                MagInt_Tensor_PiPi = (x**(2._dl*y+6._dl))*(tp1*x5 + tp2*x4 + tp3*x3 + tp4*x2 + tp5*x + tp6)
                !write(*,*) "t, MagInt_Tensor_PiPi = ", x, MagInt_Tensor_PiPi!, x**(2._dl*y+6._dl)   !,good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then !-1.55
            if (x .le. 0.0029) then
                MagInt_Tensor_PiPi = (x**(2.9_dl))*(tta*exp(ttb*x) + ttc*exp(ttd*x))
                !write(*,*) "t, MagInt_Tensor_PiPi = ", x, MagInt_Tensor_PiPi!    good
            else if ( x .GE. 0.0029 .and. x .le. 0.5) then
                MagInt_Tensor_PiPi = (x**(2.9_dl))*(tt2p1*lx9 + tt2p2*lx8 + tt2p3*lx7 + tt2p4*lx6 +&
                tt2p5*lx5 + tt2p6*lx4 + tt2p7*lx3 + tt2p8*lx2 + tt2p9*lx + tt2p10)
                !write(*,*) "2, MagInt_Tensor_PiPi = ", x, MagInt_Tensor_PiPi!   good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.5) then
                MagInt_Tensor_PiPi = (x**(2._dl*y+6._dl))*(28/(15*x) - (4*x)/105 + x2/32 - 5/4)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.001) then !0.005
                MagInt_Tensor_PiPi = ttt1p1*x3
                !write(*,*) "t, MagInt_Tensor_PiPi , = ", x, MagInt_Tensor_PiPi!    good
            else if ( x .GE. 0.001 .and. x .le. 0.004) then
                MagInt_Tensor_PiPi = (tttp1*x9 + tttp2*x8 + tttp3*x7 + tttp4*x6 + tttp5*x5&
                + tttp6*x4 + tttp7*x3 + tttp8*x2 + tttp9*x + tttp10)
                !write(*,*) "1, MagInt_Tensor_PiPi = ", x, MagInt_Tensor_PiPi!     good
            else if ( x .GE. 0.004 .and. x .le. 0.5) then
                MagInt_Tensor_PiPi = ttt2p1*x5 + ttt2p2*x4 + ttt2p3*x3 + ttt2p4*x2 + ttt2p5*x + ttt2p6
                !write(*,*) "2, MagInt_Tensor_PiPi = ", x, MagInt_Tensor_PiPi!     good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Tensor_PiPi = (x**(2._dl*y+6._dl))*(4/(63*x) + 1/(32*x2)&
                + 32/(105*x3) - 7/(12*x4) + 28/(75*x5) - 2/25)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Tensor_PiPi = (x**(2._dl*y+6._dl))*(8/(15*x5) - 7/(48*x4)&
                - 13/(960*x2) - 7/(12*x6) + 4/(15*x7) + 11/3840)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Tensor_PiPi = (x**(2._dl*y+6._dl))*(52/(10395*x3) + 148/(1575*x5)&
                - 127/(288*x6) + 556/(735*x7) - 7/(12*x8) + 28/(135*x9) - 29/22050)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Tensor_PiPi

    !---------------------------------------------------------------------------
    !> this functions computes the tensor Pi-Pi integral for helical modes
    function MagInt_Tensor_Helical_PiPi(x,y)

        implicit none

        real(dl) :: MagInt_Tensor_Helical_PiPi, x, y
        real(dl) :: x2, x3, x4, x5, x6, x7, x8, x9

        !n=, T VS TP, x=t, y=tp , N=-1.5
        !!!!t<0.00068
        !Linear model Poly9:
        !    f(x) = p1*x^9 + p2*x^8 + p3*x^7 + p4*x^6 +
        !                  p5*x^5 + p6*x^4 + p7*x^3 + p8*x^2 + p9*x + p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: thp1 =   -6.17e+18 ! (-6.56e+18, -5.781e+18)
        real(dl) :: thp2 =   2.103e+16  !(1.985e+16, 2.221e+16)
        real(dl) :: thp3 =   -3.09e+13  !(-3.238e+13, -2.942e+13)
        real(dl) :: thp4 =   2.589e+10  !(2.488e+10, 2.69e+10)
        real(dl) :: thp5 =  -1.405e+07  !(-1.446e+07, -1.365e+07)
        real(dl) :: thp6 =        5866  !(5771, 5961)
        real(dl) :: thp7 =       -5.81  !(-5.823, -5.797)
        real(dl) :: thp8 =  -3.565e-05  !(-3.654e-05, -3.476e-05)
        real(dl) :: thp9 =   3.154e-10  !(2.895e-10, 3.412e-10)
        real(dl) :: thp10 =  -8.042e-16  !(-9.341e-16, -6.742e-16)

        !Goodness of fit:
        ! SSE: 5.316e-29
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 4.281e-16

        !!!!0.04>t>0.00068
        !Linear model Poly9:
        !    f(x) = th1p1*x9 + th1p2*x8 + th1p3*x7 + th1p4*x6 +
        !                  th1p5*x5 + th1p6*x4 + th1p7*x3 + th1p8*x2 + th1p9*x + th1p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: th1p1 =  -1.221e+13  !(-1.659e+13, -7.819e+12)
        real(dl) :: th1p2 =   3.122e+11  !(2.175e+11, 4.069e+11)
        real(dl) :: th1p3 =  -3.593e+09  !(-4.476e+09, -2.71e+09)
        real(dl) :: th1p4 =   2.493e+07  !(2.027e+07, 2.958e+07)
        real(dl) :: th1p5 =  -1.208e+05  !(-1.361e+05, -1.056e+05)
        real(dl) :: th1p6 =       501.9  !(469.8, 533.9)
        real(dl) :: th1p7 =      -4.097  !(-4.14, -4.054)
        real(dl) :: th1p8 =  -0.0005457  !(-0.0005814, -0.0005099)
        real(dl) :: th1p9 =   1.245e-07  !(1.08e-07, 1.41e-07)
        real(dl) :: th1p10 =  -1.607e-11  !(-1.929e-11, -1.285e-11)
        !Goodness of fit:
        ! SSE: 1.557e-25
        ! R-square: 1
        ! Adjusted R-square: 1
        ! RMSE: 1.664e-14
        !!!!!t>0.004
        !Linear model Poly9:
        !    f(x) = p1*x^9 + p2*x^8 + p3*x^7 + p4*x^6 +
        !                  p5*x^5 + p6*x^4 + p7*x^3 + p8*x^2 + p9*x + p10
        !Coefficients (with 95% confidence bounds):
        real(dl) :: th2p1 =   -3.39e+06  !(-3.418e+06, -3.361e+06)
        real(dl) :: th2p2 =   1.134e+06  !(1.126e+06, 1.142e+06)
        real(dl) :: th2p3 =  -1.678e+05  !(-1.688e+05, -1.668e+05)
        real(dl) :: th2p4 =   1.469e+04  !(1.462e+04, 1.476e+04)
        real(dl) :: th2p5 =      -876.9  !(-879.7, -874)
        real(dl) :: th2p6 =       43.76  !(43.69, 43.84)
        real(dl) :: th2p7 =      -2.481  !(-2.483, -2.48)
        real(dl) :: th2p8 =   -0.006016  !(-0.006027, -0.006005)
        real(dl) :: th2p9 =   1.434e-05  !(1.429e-05, 1.44e-05)
        real(dl) :: th2p10 =  -1.818e-08  !(-1.829e-08, -1.807e-08)

        !Goodness of fit:
        ! SSE: 1.147e-17
        !R-square: 1
        ! Adjusted R-square: 1
        !RMSE: 3.511e-11

        !!!!!!!!!!!!! n=-1.005
        !!!!t<0.004
        !Linear model Poly9:
        !   f(x) = tthp1*x9 + tthp2*x8 + tthp3*x7 + tthp4*x6 +
        !                 tthp5*x5 + tthp6*x4 + tthp7*x3 + tthp8*x2 + tthp9*x + tthp10
        !Coefficients (with 95% confidence bounds):
        !   real(dl) :: tthp1 =   9.325e+11 ! (4.539e+11, 1.411e+12)
        !  real(dl) :: tthp2 =  -1.819e+10  !(-2.704e+10, -9.335e+09)
        ! real(dl) :: tthp3 =   1.497e+08  !(8.094e+07, 2.185e+08)
        ! real(dl) :: tthp4 =  -6.757e+05  !(-9.667e+05, -3.846e+05)
        ! real(dl) :: tthp5 =        1817  !(1088, 2546)
        !real(dl) :: tthp6 =      -1.421  !(-2.521, -0.322)
        !real(dl) :: tthp7 =     -0.6705  !(-0.6714, -0.6695)
        !real(dl) :: tthp8 =  -1.563e-06  !(-2.024e-06, -1.102e-06)
        !real(dl) :: tthp9 =   3.931e-10  !(2.931e-10, 4.931e-10)
        !real(dl) :: tthp10 =   2.368e-13  !(2.298e-13, 2.438e-13)

        !Goodness of fit:
        ! SSE: 7.211e-26
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.032e-14

        ! General model:
        !   f(x) = tthp3*x7 + tthp4*x6 + tthp5*x5 + tthp6*x4 + tthp7*x3
        !Coefficients (with 95% confidence bounds):
        real(dl) :: tthp3 =      -8.835 ! (-3128, 3110)
        real(dl) :: tthp4 =       2.627  !(-533.5, 538.8)
        real(dl) :: tthp5 =     -0.3234  !(-34.05, 33.4)
        real(dl) :: tthp6 =       1.278  !(0.3597, 2.195)
        real(dl) :: tthp7 =      -0.661  !(-0.6701, -0.6519)

        !Goodness of fit:
        ! SSE: 1.324e-10
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 1.152e-07

        !!!!!t>0.004
        !Linear model Poly8:
        !     f(x) = tth2p1*x8 + tth2p2*x7 + tth2p3*x6 + tth2p4*x5 +
        !                   tth2p5*x4 + tth2p6*x3 + tth2p7*x2 + tth2p8*x + tth2p9
        !Coefficients (with 95% confidence bounds):
        real(dl) :: tth2p1 =       46.34  !(-51.23, 143.9)
        real(dl) :: tth2p2 =      -17.38  !(-42.44, 7.679)
        real(dl) :: tth2p3 =       2.904  !(0.2214, 5.587)
        real(dl) :: tth2p4 =     -0.7574  !(-0.9126, -0.6023)
        real(dl) :: tth2p5 =       1.532  !(1.527, 1.538)
        real(dl) :: tth2p6 =     -0.6731  !(-0.6732, -0.673)
        real(dl) :: tth2p7 =  -1.508e-06  !(-2.715e-06, -3.022e-07)
        real(dl) :: tth2p8 =   5.845e-09  !(-1.237e-09, 1.293e-08)
        real(dl) :: tth2p9 =  -9.516e-12  !(-2.555e-11, 6.523e-12)

        !Goodness of fit:
        ! SSE: 6.927e-19
        !R-square: 1
        !Adjusted R-square: 1
        !RMSE: 8.627e-12
        x2 = x*x
        x3 = x*x2
        x4 = x*x3
        x5 = x*x4
        x6 = x*x5
        x7 = x*x6
        x8 = x*x7
        x9 = x*x8

        if (y==-2) then
            if (x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = (x**(2._dl*y+6._dl))*((2*x)/3 + x2/4)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==(-1.5_dl)) then
            if (x .le. 0.00068) then
                MagInt_Tensor_Helical_PiPi = thp1*x9 + thp2*x8 + thp3*x7 + thp4*x6 + &
                thp5*x5 + thp6*x4 + thp7*x3 + thp8*x2 + thp9*x + thp10
                !write(*,*) "t, MagInt_Tensor_Helical_PiPi  = ", x, MagInt_Tensor_Helical_PiPi!    good
            else if ( x .GE. 0.00068 .and. x .le. 0.004) then
                MagInt_Tensor_Helical_PiPi = th1p1*x9 + th1p2*x8 + th1p3*x7 + th1p4*x6 + th1p5*x5&
                + th1p6*x4 + th1p7*x3 + th1p8*x2 + th1p9*x + th1p10
                !write(*,*) "1, MagInt_Tensor_Helical_PiPi  = ", x, MagInt_Tensor_Helical_PiPi!    good
            else if ( x .GE. 0.004 .and. x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = th2p1*x9 + th2p2*x8 + th2p3*x7 + th2p4*x6 +&
                th2p5*x5 + th2p6*x4 + th2p7*x3 + th2p8*x2 + th2p9*x + th2p10
                !write(*,*) "2, MagInt_Tensor_Helical_PiPi = ", x, MagInt_Tensor_Helical_PiPi!    !good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y == -1) then
            if (x .le. 0.004) then
                MagInt_Tensor_Helical_PiPi = tthp3*x7 + tthp4*x6 + tthp5*x5 + tthp6*x4 + tthp7*x3
                ! tthp1*x9 + tthp2*x8 + tthp3*x7 + tthp4*x6 + tthp5*x5&
                !          + tthp6*x4 + tthp7*x3 + tthp8*x2 + tthp9*x + tthp10
                !write(*,*) "t, MagInt_Tensor_Helical_PiPi  = ", x, MagInt_Tensor_Helical_PiPi!   good
            else if ( x .GE. 0.004 .and. x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = tth2p1*x8 + tth2p2*x7 + tth2p3*x6 + tth2p4*x5 +&
                tth2p5*x4 + tth2p6*x3 + tth2p7*x2 + tth2p8*x + tth2p9
                !write(*,*) "2, MagInt_Tensor_Helical_PiPi = ", x, MagInt_Tensor_Helical_PiPi ! good
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==0) then
            if (x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = (x**(2._dl*y+6._dl))*(2/(5*x) + 1/(4*x2) - 2/(9*x3) - 1/3)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==1) then
            if (x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = (x**(2._dl*y+6._dl))*(1/(4*x4) - 1/(12*x2) &
                - 2/(15*x5) + 3._dl/320._dl)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==2) then
            if (x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = (x**(2._dl*y+6._dl))*(2/(63*x3) - 4/(25*x5) &
                + 1/(4*x6) - 2/(21*x7) - 2._dl/525._dl)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else if (y==3) then
            if (x .le. 0.5) then
                MagInt_Tensor_Helical_PiPi = (x**(2._dl*y+6._dl))*(5/(36*x6) - 1/(64*x4) &
                - 2/(7*x7) + 1/(4*x8) - 2/(27*x9) + 5/24192)
            else !t>0.5
                write(*,*) "Fitting functions not computed for t>0.5"
                stop
            end if
        else
            write(*,*) "Fitting functions not computed for this spectral index"
            stop
        end if

    end function MagInt_Tensor_Helical_PiPi


end module MagCAMB_fitting_functions
