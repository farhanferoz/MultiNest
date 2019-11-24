program main

	use params
	use nestwrapper
      
	implicit none
	
	integer i
      
      	!no parameters to wrap around
      	nest_pWrap = 0
	
	do i = 1, sdim
		sigma(i) = 0.001
	enddo

      	call nest_Sample
end
