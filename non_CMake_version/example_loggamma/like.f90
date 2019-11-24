module like

use params
use utils1
use priors
implicit none
      
contains
      
      
!=======================================================================

double precision function log_of_log_gamma(x, nu, lambda)

	implicit none
	
	double precision x, nu, lambda
	double precision z
	
	z = (x - nu) / lambda
	log_of_log_gamma = z - exp(z) - log(abs(lambda))
	
end function log_of_log_gamma
      
      
!=======================================================================

double precision function log_of_normal(x, mu, sigma)

	implicit none
	
	double precision x, mu, sigma
	double precision pi
	
	pi = 4d0 * atan(1d0)
	
	log_of_normal = -log(sigma) - log(2d0 * pi) / 2d0
	log_of_normal = log_of_normal - (((x - mu) / sigma)**2) / 2d0
	
end function log_of_normal
      
      
!=======================================================================

subroutine slikelihood(Cube, slhood)
         
	implicit none
      
	double precision Cube(nest_nPar), slhood, loglike, y, z
	integer i
	
	slhood = 0.0
	do i = 1, sdim
		Cube(i) = UniformPrior(Cube(i), center_min, center_max)
	
		if (i == 1) then
			y = log(0.5) + log_of_log_gamma(Cube(i), 10d0, 1d0)
			z = log(0.5) + log_of_log_gamma(Cube(i), -10d0, 1d0)
			loglike = LogSumExp(y, z)
		else if (i == 2) then
			y = log(0.5) + log_of_normal(Cube(i), 10d0, 1d0)
			z = log(0.5) + log_of_normal(Cube(i), -10d0, 1d0)
			loglike = LogSumExp(y, z)
		else if (i .le. (sdim + 2) / 2) then
			loglike = log_of_log_gamma(Cube(i), 10d0, 1d0)
		else
			loglike = log_of_normal(Cube(i), 10d0, 1d0)
		endif
		
		slhood = slhood + loglike
	enddo

end subroutine slikelihood
      
!=======================================================================

end module like
