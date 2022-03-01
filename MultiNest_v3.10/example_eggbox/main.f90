program main

	use params
      	use nestwrapper
      
	implicit none
      
      
	integer i
	double precision pi
	
	pi = 4d0*atan(1d0)
      
	!setting priors
	spriorran(1:sdim,1)=0d0
      	spriorran(1:sdim,2)=10d0*pi
      
      	!no parameters to wrap around
      	nest_pWrap=0

      	call nest_Sample
      	stop
end
