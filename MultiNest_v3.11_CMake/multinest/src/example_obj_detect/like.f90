module like

use params
use RandomNS
implicit none
      
double precision GLLhood0
      
contains

!----------------------------------------------------------------------
 
subroutine sCreateData(var,sdata,spredict)
         
	implicit none

	double precision var(4),sdata(snpix,snpix),tmp
	logical spredict
	integer i,j,k
	double precision gap!pixel gap in each direction between clusters
	character(len=1000) datafile

	sdata=0.
	   
	!position cluster uniformly if autopos flag is set
	if(.not.spredict) then
		if(autopos) then
			gap=snpix*0.8/(sqrt(1.*snclstr)-1.)
			do i=1,int(sqrt(1.*snclstr))
				do j=1,int(sqrt(1.*snclstr))
					spos(j*2-1+int(sqrt(1.*snclstr))*(i-1)*2)=snpix*0.1+(i-1)*gap
					spos(j*2+int(sqrt(1.*snclstr))*(i-1)*2)=snpix*0.1+(j-1)*gap
				end do
			end do
		end if
	end if
         
         
	do i=1,snpix
		do j=1,snpix
			if(.not.spredict) then
				do k=1,snclstr
					sdata(i,j)=sdata(i,j)+samp(k)*exp(-((i-0.5-spos(k*2-1))**2.+(j-0.5-spos(k*2))**2.)/(2.*(ssig(k)**2.)))
				end do
                  	else
                  		sdata(i,j)=var(3)*exp(-((i-0.5-var(1))**2.+(j-0.5-var(2))**2.)/(2.*(var(4)**2.)))
                  	end if
              	end do	
	end do
         
	if (.not.spredict) then
!		datafile='data.dat'
!            	open(2,file=datafile,form='formatted',status='replace')
            	if(dseed<0) then
        		!take the seed from system clock
        		call InitRandomNS(1)
        	else
        		call InitRandomNS(1,dseed)
	  	end if
            
            	do i=1,snpix
               		do j=1,snpix
               			tmp=Gaussian1NS(0)
                  		sdata(i,j)=sdata(i,j)+tmp*snoise
               		end do
            	end do
            
!            	do j=1,snpix
!               	do i=1,snpix
!               		write(2,*)sdata(i,j)
!               	end do
!               	write(2,*)""
!            	end do
!            
!            	close(2)
		call killRandomNS()
	end if

end subroutine sCreateData
      
!=======================================================================

subroutine sInitialize(snullev)
         
	implicit none
         
	integer i,j
	double precision Chisq
	double precision snullev
	double precision TwoPi
         
	TwoPi=6.2831853
	GLLhood0 = 0.
	GLLhood0=-dble(snpix*snpix)*log(snoise)-dble(snpix*snpix)*log(TwoPi)/2.
         
	Chisq=0.

	do i=1,snpix
		do j=1,snpix
			Chisq=Chisq+(sodata(i,j))*(1./(snoise**2))*(sodata(i,j))
            	end do
	end do

	snullev=GLLhood0-Chisq/2.
	write(*,*)"Null Evidence=",snullev

end subroutine sInitialize
      
!=======================================================================

subroutine slikelihood(Cube,slhood)
         
	implicit none
      
        double precision Cube(sdim),slhood
	double precision spdata(snpix,snpix),Chisq,temp(4)
	integer i,j

	j=0
	do i=1,4
		if(spriorran(i,1)==spriorran(i,2)) then
			temp(i)=spriorran(i,1)
            	else
            		j=j+1
            		temp(i)=(spriorran(i,2)-spriorran(i,1))*Cube(j)+spriorran(i,1)
		end if
	end do
	Cube(1:sdim)=temp(1:sdim)
	 
	call sCreateData(temp,spdata,.true.)
         
	slhood=0.
	Chisq=0.

	do i=1,snpix
		do j=1,snpix
			Chisq=Chisq+((sodata(i,j)-spdata(i,j))/snoise)**2.
            	end do
	end do
	slhood=GLLhood0-Chisq/2.d0

end subroutine slikelihood
      
!=======================================================================

end module like
