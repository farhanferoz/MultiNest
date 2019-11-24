module like

use params
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood
	double precision TwoPi
	integer i
         
	TwoPi = 6.2831853d0

	slhood = - sdim / 2d0 * log( TwoPi )
	do i = 1, sdim
		slhood = slhood - log( sigma(i) )
	enddo
	
	slhood = slhood - sum( ( ( Cube( 1:sdim ) - center ) / sigma( 1:sdim ) ) ** 2d0 ) / 2d0

end subroutine slikelihood
      
!=======================================================================

end module like
