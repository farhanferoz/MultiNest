program main

	use params
      	use nestwrapper
      
	implicit none
      
      
	integer i
      
	!setting priors
	spriorran(1:sdim,1)=-32.768
      	spriorran(1:sdim,2)=32.768
      
      	!no parameters to wrap around
      	nest_pWrap=0

      	call nest_Sample
      	stop
end
