module symmetricbetadistribution

    ! This module is an implementation of the paper:
    ! "Inverting the Symmetrical Beta Distribution" 
    !    --  PIERRE L’ECUYER and RICHARD SIMARD
    ! http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.9005&rep=rep1&type=pdf
    !
    ! The main subroutine of the module is the inversion of the cdf
    ! (cumulative distribution function) of the beta distribution:
    !
    ! f(x) = (x(1-x))^(alpha-1)  / B(alpha,alpha)
    !
    ! F(x) = int( (t(1-t))^(alpha-1)  , {t,0,x} ) / B(alpha,alpha)
    !      = B(x;alpha,alpha)/B(alpha,alpha)
    ! Note that we focus only on the algorithm for 1<alpha<200,
    ! although this module could easily be extended to lower/higher
    ! values of alpha.

    ! we denote u = F(x), and
    ! v = 1/2 - u
    ! y = 1/2 - x

    ! All series used to approximate F are computed with relative
    ! tolerance:
    double precision eps
    parameter( eps = 1d-15 )
    ! which means that we neglect all terms smaller than eps times the
    ! current sum

    ! It is useful to have the parameter sqrt(pi) precomputed
    double precision sqrtpi
    parameter(sqrtpi = 1.7724538509055160272981674833411452d0) 


    contains


    function invsymmetricbeta(u,alpha)
        ! This inverts the symmetric beta distribution
        ! F: [0,1] -> [0,1]
        ! F:   x   -> B(x,alpha,alpha)/B(alpha,alpha)
        ! using a Newton-Raphson (NR) root finder

        implicit none
        double precision             :: u,alpha
        double precision             :: invsymmetricbeta
        logical upper_half
        double precision x0,x1,y0,y1,v
        integer i,max_its
        parameter (max_its = 11)

        ! Start out with some basic checks.
        ! We actually only need to work in the domain (0,0.5), since we
        ! know the function values at 0 and 0.5, and can calculate the
        ! upper half of the range by symmetry
        if (u<0d0 .or. u>1d0) then
            write(*,*) 'invsymmetricbeta: out of bounds'
            invsymmetricbeta = 0d0
            return
        else if (u==0d0) then
            invsymmetricbeta = 0d0
            return
        else if (u==1d0) then
            invsymmetricbeta = 1d0
            return
        else if (u==0.5d0) then
            invsymmetricbeta = 0.5d0
            return
        else if (u>0.5d0) then
            ! If it lies in the upper half, use symmetry to map it to
            ! the lower half
            u = 1-u
            ! and set a flag so that we remap the answer
            upper_half= .true.
        endif

        ! The u in (0,5) can be split into two halfs, each with a
        ! preferred series.
        if ( u<u_splitter(alpha) ) then
            ! In the lower half, this is the best starting point
            x0 = (u*alpha*SymmetricBeta(alpha))**(1d0/alpha)

            do i=1,max_its
                ! iterate NR
                ! We stop after 11 iterations (which empirically never
                ! occurs)
                x1 = x0 - ( F9(x0,alpha) - u )/f(x0,alpha)

                ! if we're within tolerence, exit
                if ( .not. abovetol(x1-x0,x0) ) exit

                ! if x1 has moved outside the range, then choose a
                ! 'safer' start point just below 1/2.
                if (x1>0.5d0) then
                    x0 = 0.5d0-eps
                else
                    ! set x0 to x1 for the next iteration
                    x0 = x1
                endif
            enddo
        else
            ! In the upper half, define new variables 
            ! v= 1/2 - u
            ! y= 1/2 - x
            v = 0.5d0 - u

            ! Due to the flatness of the function around this region,
            ! this is always the best start point
            y0 = v*FouraB(alpha)

            do i=1,max_its 
                ! iterate NR
                ! We stop after 11 iterations (which empirically never
                ! occurs)
                y1 = y0 - ( H10(y0,alpha) - v )/g(y0,alpha)

                ! if we're within tolerence, exit
                if ( .not. abovetol(y1-y0,y0) ) exit

                ! set y0 to y1 for the next iteration
                y0 = y1
            enddo

            ! transform from y to x
            x1 = 0.5d0 - y1

        endif

        if (upper_half) then
            ! if we were working in the upper half map the result onto
            ! (0.5,1)
            invsymmetricbeta = 1 - x1
            ! and re-set u back to its original value
            u=1-u
        else
            invsymmetricbeta = x1
        endif
            
    end function invsymmetricbeta 



    function u_splitter(alpha)
        ! It turns out that this empirical curve is a very good
        ! separator for the regions where the two series are most
        ! useful for alpha>1
        implicit none
        double precision, intent(in) :: alpha
        double precision u_splitter

        u_splitter = 1d0 / (2.5 + 2.25*sqrt(alpha) )
    end function



    function F9 (x,alpha)
        ! This denotes the function in equation (9) of the paper. 
        ! It is an (optimal) power series expansion of the inverse
        ! beta distribution u
        implicit none
        double precision, intent(in) :: x,alpha
        double precision             :: F9

        F9 = x**alpha * (1-x)**(alpha-1) / (alpha * SymmetricBeta(alpha) ) &
             * Hypergeometric2F1(1d0, 1-alpha, 1+alpha, -x/(1-x) )

    end function F9

    function f (x,alpha)
        ! This is simply the (unintegrated) probability distribution
        implicit none
        double precision, intent(in) :: x,alpha
        double precision             :: f

        f = (x*(1-x))**(alpha-1) / SymmetricBeta(alpha)
    end function f

    function H10 (y,alpha)
        ! This denotes the function in equation (10) of the paper. 
        ! It is an (optimal) power series expansion of the inverse
        ! beta distribution v = 1/2 - u in terms of y=1/2 - x
        implicit none
        double precision, intent(in) :: alpha,y
        double precision             :: H10

        H10 = y * (1-4*y**2)**alpha / FouraB(alpha) &
              * Hypergeometric2F1(0.5 + alpha,1d0,1.5d0,4*y*y)

    end function H10

    function g (y,alpha)
        ! This is simply the (unintegrated) probability distribution,
        ! with y=1/2 - x as a variable
        implicit none
        double precision, intent(in) :: y,alpha
        double precision             :: g

        g = (0.25d0 - y**2 )**(alpha-1) / SymmetricBeta(alpha)
    end function g


    function SymmetricBeta(alpha)
        ! compute B(alpha,alpha) where B is the beta function:
        ! B(a,b) = gamma(a)gamma(b)/gamma(a+b)
        ! This is a normalising factor 
        implicit none
        double precision, intent(in) :: alpha
        double precision             :: SymmetricBeta

        SymmetricBeta = gamma(alpha)**2/gamma(2*alpha)

    end function SymmetricBeta

    function FouraB(alpha)
        ! Care must be taken in computing the expression:
        ! 4^(alpha−1)*B(alpha,alpha)
        ! A naive calculaton of the two factors separately will give
        ! rise to catastrophic loss of precision for large alpha
        implicit none
        double precision,intent(in) :: alpha
        double precision            :: FouraB

        FouraB = sqrtpi / 2d0

        if (alpha < 10) then
            FouraB = FouraB *  gamma(alpha) / gamma(alpha + 0.5d0)
        else
            FouraB = FouraB / sqrt( (alpha-0.5d0) * Hypergeometric2F1(-0.5d0,-0.5d0,alpha-0.5d0,1d0) ) 
        endif

    end function FouraB


    function Hypergeometric2F1(a,b,c,z)
        implicit none
        double precision, intent(in)  :: a,b,c,z
        double precision              :: Hypergeometric2F1
        integer n
        double precision change

        ! This computes the hypergeometric 2F1 function using a
        ! truncated power series, stopping when the relative change
        ! due to higher terms is less than epsilon

        ! http://en.wikipedia.org/wiki/Hypergeometric_function#The_hypergeometric_series

        Hypergeometric2F1=0d0
        n=0

        do while( abovetol(change,Hypergeometric2F1) )
           change = Pochhammer(a,n) * Pochhammer(b,n) / Pochhammer(c,n) * z**n / gamma(1d0+n)

           Hypergeometric2F1 = Hypergeometric2F1 + change
           n=n+1
        enddo

    end function Hypergeometric2F1




    recursive function Pochhammer (x,n) result (xn)
        ! This function computes the rising factorial x^(n):
        ! for a non-negative integer n and real x
        !
        ! http://en.wikipedia.org/wiki/Pochhammer_symbol

        implicit none
        double precision, intent(in)  :: x
        integer,          intent(in)  :: n
        double precision              :: xn

        if (n<=0) then
            xn = 1
        else
            xn = Pochhammer(x,n-1)*(x+n-1)
        endif

    end function Pochhammer

    function abovetol (change,current_sum)
        ! This function outputs true if change is outside of the
        ! tolerance epsilon from current_sum
        implicit none
        double precision, intent(in) :: change,current_sum
        logical abovetol

        if (current_sum == 0d0) then
            ! this check is useful for entering loops
            abovetol = .true.
        else
            abovetol = abs( change/current_sum ) > eps
        endif

    end function abovetol

end module 
