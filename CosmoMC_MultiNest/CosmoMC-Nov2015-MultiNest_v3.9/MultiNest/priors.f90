module priors

    use utils1
    use symmetricbetadistribution

    contains

    !=======================================================================
    ! Prior distribution functions: r is a uniform deviate from the unit
    !
    !

    ! Uniform[0:1]  ->  Delta[x1]

    function DeltaFunctionPrior(r,x1,x2)

        implicit none

        double precision r,x1,x2,DeltaFunctionPrior

        DeltaFunctionPrior=x1

    end function DeltaFunctionPrior

    !=======================================================================
    ! Uniform[0:1]  ->  Uniform[x1:x2]

    function UniformPrior(r,x1,x2)

        implicit none

        double precision r,x1,x2,UniformPrior

        UniformPrior=x1+r*(x2-x1)

    end function UniformPrior

    !=======================================================================
    ! Uniform[0:1]^d  -> Sorted Uniform[x1:x2]^d

    subroutine UniformSortPrior(r,r_min,r_max)
        ! This transforms the unit hypercube to a "forced identifiablity" prior.
        ! This means that the (r1,r2,...,rn) variables are uniformly distributed in the
        ! physical prior space between r_min and r_max, but have been sorted so
        ! that r1<r2<...<rn. This amounts to choosing a non-separable prior such that:
        !
        ! P(r1)       =  n/(r_max-r_min) * ((r_max-r_min)-r1)^(n-1)  r_min<r1<r_max
        ! P(r2|r1)    = (n-1)/(r_max-r1) * ( (r_max-r1) - r1)^(n-2)     r1<r2<r_max
        ! P(r2|r1,r2) = (n-2)/(r_max-r2) * ( (r_max-r2) - r1)^(n-3)     r2<r3<r_max
        ! ...
        ! P(rn|r1,...,rn-1) =  1/(r_max-rn-1)                         rn-1<rn<r_max
        !
        ! The first of these is simply the probability density for the smallest
        ! of n points in the range [r_min,r_max]
        !
        ! For the next highest point it is the smallest of n-1 points in the
        ! range [r1,r_max]. And so on.
        !
        !
        ! To perform this transformation, it is cleanest to do this in two
        ! steps. First transform the variable r(i) on the interval [0,1] to the
        ! minimum of n-(i-1) randomly drawn points in [0,1]. This is done via the inverse of
        ! the CDF:
        ! CDF(x) = 1 - (1-x)^(n-(i-1))  =>   r(i) -> 1-(1-r(i))^( 1/(n-(i-1)) )
        !
        ! Second, this is then transform uniformly from [0,1] to [r(i-1),r_max] 

        implicit none

        double precision, intent(inout) :: r(:) ! input:  multinest Hypercube parameters
        ! output: physical cosmological parameters
        double precision r_min ! The lower limit on the range
        double precision r_max ! The upper limit on the range

        ! working variables

        integer n_prior ! the dimension
        integer i       ! iteration counter
        double precision r_end ! saves the varible r(i-1) 


        ! Get the size of the array
        n_prior = size(r)

        ! Initialise the counter variables
        r_end = r_min
        i=1


        ! generate the smallest of n variables
        do while (n_prior>0)

            ! Transform from a uniform random variable in [0,1]
            ! to the smallest of d=n_prior variables in [0,1]^d
            r(i)  = 1 - ( 1 - r(i) )**(1.d0/n_prior)

            ! Transform from [0,1] to [r_end,r_max]
            r(i)  = UniformPrior(r(i),r_end,r_max)

            ! Save this variable for the next iteration
            r_end = r(i)

            ! decrease the number of variables
            n_prior = n_prior-1
            ! move on to the next loop
            i=i+1
        end do








    end subroutine UniformSortPrior

    !=======================================================================
    ! Uniform[0:1]  ->  LogUniform[x1:x2]

    function LogPrior(r,x1,x2)

        implicit none

        double precision r,x1,x2,LogPrior
        double precision lx1,lx2

        if (r.le.0.0d0) then
            LogPrior=-1.0d32
        else
            lx1=dlog10(x1)
            lx2=dlog10(x2)
            LogPrior=10.d0**(lx1+r*(lx2-lx1))
        endif

    end function LogPrior

    !=======================================================================
    ! Uniform[0:1]  ->  Sin[x1:x2]  (angles in degrees):

    function SinPrior(r,x1,x2)

        implicit none

        double precision r,x1,x2,SinPrior
        real cx1,cx2,deg2rad
        parameter(deg2rad=0.017453292)

        cx1=cos(x1*deg2rad)
        cx2=cos(x2*deg2rad)
        SinPrior=1.d0*acos(cx1+r*(cx2-cx1))

    end function SinPrior

    !=======================================================================
    ! Uniform[0:1]  ->  Cauchy[mean=x0,FWHM=2*gamma]

    function CauchyPrior(r,x0,gamma)

        implicit none

        double precision r,x0,gamma,CauchyPrior
        real Pi
        parameter(Pi=3.141592654)

        CauchyPrior=x0+gamma*tan(Pi*(r-0.5))

    end function CauchyPrior

    !=======================================================================
    ! Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]

    function GaussianPrior(r,mu,sigma)

        implicit none

        double precision r,mu,sigma,GaussianPrior
        double precision SqrtTwo
        parameter(SqrtTwo=1.414213562d0)

        if (r.le.1.0d-16.or.(1.0d0-r).le.1.0d-16) then
            GaussianPrior=-1.0d32
        else
            GaussianPrior=mu+sigma*SqrtTwo*dierfc(2.d0*(1.d0-r))
        endif

    end function GaussianPrior

    !=======================================================================
    ! Uniform[0:1]  ->  Ellipsoid[minimum=x1,maximum=x2]^n

    function EllipsoidPrior(r,x1,x2,n)

        implicit none

        double precision r,x1,x2,EllipsoidPrior
        integer n
        integer ifault
        integer u

        EllipsoidPrior = invsymmetricbeta(r,1d0+n/2d0)
        EllipsoidPrior = x1 + (x2-x1) * EllipsoidPrior

    end function EllipsoidPrior
    !=======================================================================
    ! Uniform[0:1]^d  ->  Gaussian[mean=mean_vec_i,variance=cov_mat_ij]

    subroutine MultivariateGaussianPrior(r,mean_vec,cov_mat)
        ! This function converts a vector of uniformly distrited hypercube
        ! parameters r(:) from multinest into a multivariate gaussian
        ! distribution of theta_i:
        !
        !  P(theta_i) ~ exp(-1/2 * (theta_i - mu_i) (Sigma^-1_ij) (theta_j - mu_j))
        !
        ! where: 
        !   theta_i are cosmological parameters, 
        !   mu_i are their means,
        !   sigma_ij are their covariance matrix.
        !
        ! The subroutine works by transforming the vector of numbers it is
        ! given in the variable 'r(:)' according to the covariance matrix and
        ! mean it is given.
        !
        ! The means and covariance matrices are implicitly used via the
        ! 'marginalBm' variable.
        !
        ! The trick to converting is to use the fact that a non-separable prior
        ! can be separated using:
        !
        ! pi(x) = pi(x1) * pi(x2|x1) * pi(x3|x1,x2) * ... * pi(xn|x1,...,xn-1)
        !
        ! where the general term pi(xj|x1,...,xj-1) is found by marginalising
        ! over the final j+1 variables, and then conditioning on the first j-1.
        ! For gaussian distributions, the process of marginalising and
        ! conditioning is very simple (analytically speaking. For a very nice
        ! demonstration see:
        !
        !   http://gbhqed.wordpress.com/2010/02/21/conditional-and-marginal-distributions-of-a-multivariate-gaussian/
        !
        ! and for a more formal discussion see:
        !
        !   http://orbit.dtu.dk/fedora/objects/orbit:78104/datastreams/file_2691604/content
        !
        ! The upshot of this is that to marginalise over the final j+1
        ! variables, you just exclude the lower three 'quarters' of the
        ! covariance matrix:
        !
        !   [ . . . . . . x x x ]  
        !   [ . . . . . . x x x ]  
        !   [ . . . . . . x x x ]   
        !   [ . . . . . . x x x ]   
        !   [ . . . . . . x x x ]   
        !   [ . . . . .jj x x x ]   
        !   [ x x x x x x x x x ]   
        !   [ x x x x x x x x x ]   
        !   [ x x x x x x x x x ]   
        !
        ! In order to condition on the first j-1 variables, we split the
        ! remaining matrix into 4 parts:
        !
        !   [ . . . . . | x x x ]  
        !   [ . . . . . | x x x ]  
        !   [ . . B . . C x x x ]   
        !   [ . . . . . | x x x ]   
        !   [ . . . . . | x x x ]   
        !   [ ---C^T---jj x x x ]   
        !   [ x x x x x x x x x ]   
        !   [ x x x x x x x x x ]   
        !   [ x x x x x x x x x ]   
        !
        ! A (j-1)x(j-1) matrix B, two identical vectors (C and its transpose),
        ! and the jj diagonal element. 
        ! It is now a standard result that the distribution pi(xj|x1,...xj-1) is
        ! gaussian with 
        !   mean     = mu_j + C^T B^-1 (x1-mu1, ... , xj-1 - muj-1)^T
        !   variance = sigma_jj - C^T B^-1 C
        !
        ! Since we take the covariance matrix as an input to the program, the
        ! only non-trivial things we need to calculate are the B^-1 matrices.


        implicit none
        double precision, intent(inout) :: r(:) ! input:  multinest Hypercube parameters
        ! output: physical cosmological parameters

        double precision, intent(in) :: mean_vec(:)  ! The means of the gaussian
        double precision, intent(in) :: cov_mat(:,:) ! The covariance matrix of the gaussian


        ! working variables
        double precision,allocatable :: theta(:)  ! The physical parameter working variable
        integer n_prior ! the dimension of the gaussian
        double precision mean     
        double precision variance 
        double precision sigma     ! Standard deviation
        integer i

        double precision,allocatable :: margeBm(:,:)
        ! margeBm(j,:,:) corresponds to an (j-1)x(j-1)
        ! matrix which is the inverse of the matrix B detailed above.


        ! Get the dimension of the input array
        n_prior = size(cov_mat,DIM=1)
        allocate(theta(n_prior))

        ! Check for any issues with array shape
        if (     size(r,DIM=1)        /= n_prior   &
            .or. size(mean_vec,DIM=1) /= n_prior   &
            .or. size(cov_mat,DIM=2)  /= n_prior ) &
            write(*,*)"Warning:  MultivariateGaussianPrior: input arrays do not match in size"




        ! For each of the variables in r(:)
        do i=1,n_prior

            ! a) All means and variances are these with a modifier
            mean = mean_vec(i)
            variance = cov_mat(i,i)

            if(i>1) then     ! (if i=1, then there is no matrix B)

                ! b) Invert the matrix B (as explained above)
                allocate(margeBm(i-1,i-1))

                ! invert the cov_mat and store in margeBm
                margeBm = cov_mat(:i-1,:i-1)
                call invert_symmetric_matrix(margeBm)

                ! c) Modify the mean and variance according to:
                !      mean     = mu_j + C^T B^-1 (x1-mu1, ... , xj-1 - muj-1)^T
                !      variance = sigam_jj - C^T B^-1 C
                mean = mean + dot_product( cov_mat(i,:i-1) ,matmul(margeBm(:i-1,:i-1),theta(:i-1) - mean_vec(:i-1) ) )
                variance = variance + dot_product( cov_mat(i,:i-1) ,matmul(margeBm(:i-1,:i-1), cov_mat(i,:i-1) ) )

                deallocate(margeBm)

            endif

            ! d) Calculate the variable theta_i from the variable r(:)
            sigma = sqrt(variance)
            theta(i) = GaussianPrior(r(i),mean,sigma)

        enddo

        ! Pass theta back to r(:) for outputting (theta is local)
        r = theta

    end subroutine MultivariateGaussianPrior

    !=======================================================================
    ! Uniform[0:1]^d  ->  UniformEllipsoid[mean=mean_vec_i,variance=cov_mat_ij]

    subroutine EllipsoidalUniformPrior(r,mean_vec,cov_mat)
        ! This function converts a vector of uniformly distrited
        ! hypercube parameters r(:) from multinest into a prior which
        ! is uniform on the ellipsoidal bound defined by the
        ! covariance matrix distribution of theta_i:
        !
        !  P(theta_i) ~ 1/Volume(Cov) : (theta_^T-mu_) Cov^-1 (theta_-mu_) < 1
        !               0             : otherwise
        !
        ! where: 
        !   theta_ are cosmological parameters, 
        !   mu_ are their means,
        !   Cov is their covariance matrix.
        !
        ! The subroutine works by transforming the vector of numbers
        ! it is given in the variable 'r(:)' according to the
        ! covariance matrix and mean it is given.
        !
        ! The means and covariance matrices are implicitly used via the
        ! 'marginalBm' variable.
        !
        ! The trick to converting is to use the fact that a
        ! non-separable prior can be separated using:
        !
        ! pi(x) = pi(x1) * pi(x2|x1) * pi(x3|x1,x2) * ... * pi(xn|x1,...,xn-1)
        !
        ! where the general term pi(xj|x1,...,xj-1) is found by
        ! marginalising over the final j+1 variables, and then
        ! conditioning on the first j-1.
        !
        ! We aim to calculate p( theta_i | theta_1,...,theta_i-1 )
        !
        ! We begin by splitting up the vector of variables into two
        ! vectors y,z.
        !
        !               [  y  ]
        ! theta - mu =  [ --- ]
        !               [  z  ]
        !                      
        !        [ theta_1 - mu_1 ]          [ theta_i+1 - mu_i+1 ]  
        !    y = [ theta_2 - mu_2 ],     z = [ theta_i+2 - mu_i+2 ]  
        !        [       ||       ]          [         ||         ]  
        !        [ theta_i - mu_i ]          [  theta_n  -  mu_n  ]  
        !
        ! The inverse of the covariance matrix will split into four
        ! accordingly:
        !
        !          [  A  | B^T  ]  
        ! Cov^-1 = [ ----+----- ]  
        !          [  B  | C^-1 ]  
        !
        ! and the condition for the edge of the ellipsoidal
        ! distribution:
        !
        ! theta^T Cov^-1 theta = 1
        !
        ! becomes:
        !
        ! ( z + CBy )^T C^-1 ( z + CBy ) = 1 - y^T (A- B^T C B ) y 
        !
        ! This shows that at fixed y, the vector z describes an
        ! off-center ellipsoid with a covariance matrix:
        ! 
        ! (1 - y^T D y ) C   , where  D = A - B^T C B 
        !
        ! When you marginalise over the z, you are left with a
        ! distribution which is proportional (removing factors
        ! independent of y) to
        !
        ! (1 - y^T D y )^((n-i)/2) 
        !
        ! since the matrix C has dimension (n-i)x(n-i)
        !
        ! We now separate theta_i from the rest of y in order to
        ! condition on the others
        !
        !     [   theta_1 - mu_1   ]    [  w  ]
        ! y = [         ||         ] =  [ --- ]
        !     [ theta_i-1 - mu_i-1 ]    [  x  ]
        !     [ ------------------ ]
        !     [   theta_i - mu_i   ]
        !
        ! And D separates accordingly
        !
        !     [  M  | v^T ]
        ! D = [ ----+---- ]
        !     [  v  |  p  ]
        !
        ! So we find that:
        !   (1 - y^T D y )^ (n-i)/2 
        ! = ( 1 - w^T M w - (v^t w + w^t w) x - p x^2 )^ (n-i)/2
        !
        ! If we fix the vector w, and just allow x to vary, we find
        ! that this is a scaled, symmetric 'Beta distribution' with
        ! parameter alpha = 1 + (n-i)/2, so that it is non-zero in the
        ! region x_- < x < x_+, with
        !
        ! x_+/-  =  (-v^T w +/- sqrt(p + (v^T w )^2 - p w^T M w ) ) /p
        !
        ! We don't need to worry about normalisation as re-scaled beta
        ! distributionis defined by its zeroes. The aim of this
        ! algorithm is therefore to calculate the positions of these
        ! zeroes




        implicit none
        double precision, intent(inout) :: r(:) ! input:  multinest Hypercube parameters
        ! output: physical cosmological parameters

        double precision, intent(in) :: mean_vec(:)  ! The means of the gaussian
        double precision, intent(in) :: cov_mat(:,:) ! The covariance matrix of the gaussian


        ! working variables
        double precision,allocatable :: theta(:)  ! The physical parameter working variable
        integer n_prior ! the dimension of the gaussian
        integer i
        double precision disc,x1,x2

        double precision,allocatable,dimension(:,:) :: covinv,A,B,Cinv,C,D
        double precision,allocatable,dimension(:)   :: w

        double precision,allocatable, save :: M(:,:,:), v(:,:), p(:)
        logical :: calculate

        double precision,allocatable,dimension(:) :: mean_vec_saved  ! The means of the gaussian
        double precision,allocatable,dimension(:,:) :: cov_mat_saved ! The covariance matrix of the gaussian



        ! Get the dimension of the input array
        n_prior = size(cov_mat,DIM=1)
        allocate(theta(n_prior))
        allocate(w(n_prior))

        ! Check for any issues with array shape
        if (     size(r,DIM=1)        /= n_prior   &
            .or. size(mean_vec,DIM=1) /= n_prior   &
            .or. size(cov_mat,DIM=2)  /= n_prior ) &
            write(*,*)"Warning:  EllipsoidalUniformPrior: input arrays do not match in size"


        if ( (.not. allocated(mean_vec_saved) ) .or. (.not. allocated(cov_mat_saved) ) ) then
            calculate = .true.
            allocate(mean_vec_saved(n_prior))
            allocate(cov_mat_saved(n_prior,n_prior))

        else if ( all(cov_mat_saved .eq. cov_mat ) .and. all(mean_vec_saved .eq. mean_vec ) ) then
            calculate = .false.

        else
            calculate = .true.
        endif

        ! Save the covariance matrix and mean vector for checking next time this
        ! function is called
        cov_mat_saved = cov_mat
        mean_vec_saved = mean_vec


        if (calculate) then
            ! Allocate arrays as necessary
            allocate(covinv(n_prior,n_prior))
            allocate(A(n_prior,n_prior))
            allocate(B(n_prior,n_prior))
            allocate(Cinv(n_prior,n_prior))
            allocate(C(n_prior,n_prior))
            allocate(D(n_prior,n_prior))

            if (allocated(M)) deallocate(M)
            allocate(M(n_prior,n_prior,n_prior))

            if (allocated(v)) deallocate(v)
            allocate(v(n_prior,n_prior))

            if (allocated(p)) deallocate(p)
            allocate(p(n_prior))



            ! initialise the inverse covariance matrix
            covinv = cov_mat
            ! Invert the covariance matrix
            call invert_symmetric_matrix(covinv)
        endif


        ! For each of the variables in r(:)
        do i=1,n_prior

            if ( calculate ) then
                ! Split up the covariance matrix into A,B,C^-1
                A(:i,:i) = covinv(:i,:i)
                D(:i,:i) = A(:i,:i)

                if (i/=n_prior) then
                    B(:n_prior-i,:i)            = covinv(i+1:,:i)
                    Cinv(:n_prior-i,:n_prior-i) = covinv(i+1:,i+1:)

                    ! invert Cinv
                    C(:n_prior-i,:n_prior-i) = Cinv(:n_prior-i,:n_prior-i)
                    call invert_symmetric_matrix(C(:n_prior-i,:n_prior-i))

                    ! Calculate D =  A - B^T C B 
                    ! --------------------------

                    D(:i,:i) = D(:i,:i) -                 &
                        matmul(matmul(                        &
                        transpose( B(:n_prior-i,:i) ),      &
                        C(:n_prior-i,:n_prior-i)         ), &
                        B(:n_prior-i,:i)                 )
                endif


                ! separate D into M,v,p
                if (i/=1) then
                    M(:i-1,:i-1,i) = D(:i-1,:i-1)
                    v(:i-1,i)      = D(i,:i-1) 
                endif
                p(i) = D(i,i)


            endif

            ! calculate the x_+/- values
            !
            ! x_+/-  =  (-v^T w +/- sqrt( p + (v^T w )^2 -  p w^T M w ) ) / p
            !
            ! Note that we have saved M,v and p from previous calcultions

            x2 = mean_vec(i)

            disc = p(i)     ! initialise the discriminant

            if (i/=1) then
                w(:i-1) = theta(:i-1) - mean_vec(:i-1)
                disc = disc + dot_product(v(:i-1,i),w(:i-1)) **2  
                disc = disc - p(i) * dot_product( w(:i-1) , matmul( M(:i-1,:i-1,i) , w(:i-1) ) )
                x2= x2 - dot_product(v(:i-1,i),w(:i-1)) / p(i)
            endif

            disc = sqrt(disc) / p(i)

            x1= x2 - disc
            x2= x2 + disc

            ! Convert r_i to theta_i
            theta(i)=EllipsoidPrior(r(i),x1,x2,n_prior-i)

        enddo

        ! Pass theta back to r(:) for outputting (theta is local)
        r = theta

    end subroutine EllipsoidalUniformPrior

    !=======================================================================
    subroutine invert_symmetric_matrix(matrix)
        ! using Lapack routines: http://www.netlib.org/lapack/double/        
        implicit none

        double precision, intent(inout) :: matrix(:,:)
        integer j,k,n,information

        n = size(matrix,DIM=1) 
        if (n /= size(matrix,DIM=2) ) write(*,*) 'invert_matrix: n/=m'
        do j=1,n
            do k=1,n
                if ( matrix(k,j) /= matrix(j,k) ) write(*,*) 'Matrix not symmetric!'
            enddo
        enddo

        ! Convert matrix to Cholesky factorisation of a symmetric matrix
        ! http://www.netlib.org/lapack/double/dpotrf.f
        call DPOTRF('L',n,matrix,n,information)

        ! Convert Cholesky factorisation to the inverse matrix
        ! http://www.netlib.org/lapack/double/dpotri.f
        call DPOTRI('L',n,matrix,n,information)

        ! Convert from lower triangle to symmetric matrix
        do j=1,n
            do k=1,j-1
                matrix(k,j) = matrix(j,k)
            enddo
        enddo

    end subroutine invert_symmetric_matrix

    !=======================================================================
    ! Uniform[0:1]  ->  LogNormal[mode=a,width parameter=sigma]

    function LogNormalPrior(r,a,sigma)

        implicit none

        double precision r,a,sigma,LogNormalPrior
        double precision SqrtTwo,bracket
        parameter(SqrtTwo=1.414213562d0)

        bracket=sigma*sigma+sigma*SqrtTwo*dierfc(2.d0*r)
        LogNormalPrior=a*dexp(bracket)

    end function LogNormalPrior

    !=======================================================================       
    ! Inverse of complimentary error function in double precision

    function dierfc(y)

        implicit none
        double precision y,dierfc
        double precision qa,qb,qc,qd,q0,q1,q2,q3,q4,pa,pb,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18
        double precision p19,p20,p21,p22,x,z,w,u,s,t
        double precision infinity
        parameter (infinity=5.0d0)
        parameter (qa=9.16461398268964d-01, &
            qb=2.31729200323405d-01, &
            qc=4.88826640273108d-01, &
            qd=1.24610454613712d-01, &
            q0=4.99999303439796d-01, &
            q1=1.16065025341614d-01, &
            q2=1.50689047360223d-01, &
            q3=2.69999308670029d-01, &
            q4=-7.28846765585675d-02)
        parameter (pa=3.97886080735226000d+00, &
            pb=1.20782237635245222d-01, &
            p0=2.44044510593190935d-01, &
            p1=4.34397492331430115d-01, &
            p2=6.86265948274097816d-01, &
            p3=9.56464974744799006d-01, &
            p4=1.16374581931560831d+00, &
            p5=1.21448730779995237d+00, &
            p6=1.05375024970847138d+00, &
            p7=7.13657635868730364d-01, &
            p8=3.16847638520135944d-01, &
            p9=1.47297938331485121d-02, &
            p10=-1.05872177941595488d-01, &
            p11=-7.43424357241784861d-02)
        parameter (p12=2.20995927012179067d-03, &
            p13=3.46494207789099922d-02, &
            p14=1.42961988697898018d-02, &
            p15=-1.18598117047771104d-02, &
            p16=-1.12749169332504870d-02, &
            p17=3.39721910367775861d-03, &
            p18=6.85649426074558612d-03, &
            p19=-7.71708358954120939d-04, &
            p20=-3.51287146129100025d-03, &
            p21=1.05739299623423047d-04, &
            p22=1.12648096188977922d-03)
        if (y==0.0) then
            dierfc=infinity
            return
        endif  
        z=y
        if (y .gt. 1) z=2-y
        w=qa-log(z)
        u=sqrt(w)
        s=(qc+log(u))/w
        t=1/(u+qb)
        x=u*(1-s*(0.5d0+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t
        t=pa/(pa+x)
        u=t-0.5d0
        s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12
        s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2) &
            *u+p1)*u+p0)*t-z*exp(x*x-pb)
        x=x+s*(1+x*s)
        if (y .gt. 1) x=-x
        dierfc=x

    end function dierfc

    !=======================================================================
end module priors
