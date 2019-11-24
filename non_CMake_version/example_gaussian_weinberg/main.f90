program main

	use params
	use nestwrapper
	use RandomNS
      
	implicit none
	
	integer i, j
	
	double precision urv, alpha_sum
	
	!initialise the random number generator
	if(data_seed < 0) then
        	!take the seed from system clock
        	call InitRandomNS(1)
        else
        	call InitRandomNS(1, data_seed)
	end if
	
	!set the weights & centers of the Gaussian components in the mixture model
	do i = 1, num_components
		do j = 1, ndim_components
			urv = ranmarns(0)
			centers(i, j) = UniformPrior(urv, center_min, center_max)
		enddo
		
		urv = ranmarns(0)
		alpha(i) = ExponentialPrior(urv, 1d0)
	enddo
	alpha_sum = sum(alpha)
	alpha = alpha / alpha_sum
	
	
	call killRandomNS()

      	call nest_Sample
end
