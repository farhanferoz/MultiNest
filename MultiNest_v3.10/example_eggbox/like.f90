module like
	
use params
use utils1
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood
	integer i
	
	slhood = 1d0
	
	do i = 1, sdim
		Cube(i) = spriorran(i, 1) + ( spriorran(i, 2) - spriorran(i, 1) ) * Cube(i)
		slhood = slhood * cos( Cube(i)/2d0 )
	enddo
	
	slhood = (2d0 + slhood)**5d0

end subroutine slikelihood
      
!=======================================================================

end module like
