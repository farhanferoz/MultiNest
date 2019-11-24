module like
	
use params
use utils1
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood
	double precision a,b,pi
	integer i,j
	
	pi = 4d0 * atan(1d0)
	
	do i = 1, sdim
		Cube(i) = spriorran(i, 1) + ( spriorran(i, 2) - spriorran(i, 1) ) * Cube(i)
	enddo
	
	a = 0d0
	b = 0d0
	do i = 1, sdim
		a = a + Cube(i)**2
		b = b + cos( 2 * pi * Cube(i) )
	enddo
	a = -0.2 * sqrt( a / sdim )
	b = b / sdim
	
	slhood = 20d0 * exp(a) + exp(b) - 20d0 - exp(1d0)

end subroutine slikelihood
      
!=======================================================================

end module like
