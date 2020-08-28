program main

	use initialize
	implicit none
	integer :: i,j,k,monteCarloSteps = 10000
	real*8 :: pe, TE, ii,jj,kk
	real :: speed
	
	call posInit()
	call initForce()
	call initVelocities()
	
	open(81, file = 'speed3.dat')
	
	thermoCounter = 500
	
	do j = 1, monteCarloSteps
		call verletAlgorithm()
		do i = 1, particles
			speed = sqrt(velocity_x(i)**2 + velocity_y(i)**2 + velocity_z(i)**2)		
			write(81, *) speed
		end do 
	end do
	
	close(81)
	
end program main 


