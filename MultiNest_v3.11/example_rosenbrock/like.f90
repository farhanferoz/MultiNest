module like
	
use params
use utils1
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
	double precision Cube(nest_nPar),slhood
	integer i,j
	
	slhood = 0d0
	
	do i = 1, sdim
		Cube(i) = spriorran(i, 1) + ( spriorran(i, 2) - spriorran(i, 1) ) * Cube(i)
	enddo
	
	do i = 1, sdim - 1
		slhood = slhood - ( ( ( 1d0 - Cube(i) ) ** 2 ) + 100d0 * ( ( Cube(i+1) - ( Cube(i) ** 2 ) ) ** 2 ) )
	enddo

end subroutine slikelihood
      
!=======================================================================

end module like
