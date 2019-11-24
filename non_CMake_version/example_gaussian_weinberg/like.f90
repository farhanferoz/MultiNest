module like

use params
use utils1
use priors
implicit none
      
contains
      
      
!=======================================================================

subroutine slikelihood(Cube, slhood)
         
	implicit none
      
	double precision Cube(nest_nPar), slhood, loglike
	double precision TwoPi
	integer i, j
         
	twoPi = 8d0 * atan(1d0)
	
	slhood = 0.0
	do i = 1, num_components
		loglike = log(alpha(i)) - (ndim_components / 2.0) * log(twoPi * variance)
		loglike = loglike - sum(((Cube(:) - centers(i, :))**2) / (2*variance))
		
		if (i == 1) then
			slhood = loglike
		else
			slhood = LogSumExp(loglike, slhood)
		endif
	enddo

end subroutine slikelihood
      
!=======================================================================

end module like
