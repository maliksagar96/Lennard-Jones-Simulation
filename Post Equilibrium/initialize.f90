module initialize

	  real*8 :: length = 15.0d0                                         !Box is cubical box of Boxsize = boxsize*boxsize*boxsize
    integer, parameter :: particles = 600
    real, parameter :: dt = 0.0050d0
    real*8, parameter :: sigma = 1.0d0, eps = 1.0d0, mass = 1.0d0, temp = 1.0d0
    real*8, parameter :: rc = 2.50d0*sigma                            !Cutoff Distance
    real*8 :: velocity_x(particles), velocity_y(particles), velocity_z(particles)
    real*8 :: x(particles), y(particles), z(particles)
    integer, parameter :: seed_x = 2550, seed_y = 456, seed_z = 588 !To check random number calling for differernt seeds
    real*8 :: dt2by2 = (dt**2.0d0)/(mass*2.0d0), dtby2 = dt/(2.0d0*mass)
    real*8 :: potential = 0.0d0, kineticEnergy = 0.0d0
    real*8 :: force_x(particles), force_y(particles), force_z(particles)
    real*8 :: newForce_x(particles), newForce_y(particles), newForce_z(particles)
    real*8 :: f_rc = eps*((12.0d0*(sigma**12)/(rc**13)) -(6.0d0*(sigma**6)/rc**7))
    real*8 :: u_rc = eps*((sigma**12)/(rc**12) - (sigma**6)/(rc**6))
    real*8 :: KE_initial
    integer :: updateCounter, thermoCounter

		contains
		
		    !Minimum image convention
    function minDist(x1, x2) result(dist)
      implicit none
      real*8, intent(in) :: x1, x2
      real*8 :: dist

        if((x1 - x2) > length /2.0d0) then
          dist = x1 -x2 - length

        else if ((x1 - x2)< (-length /2.0d0)) then
          dist = x1 - x2 + length

        else
          dist = x1 - x2
        end if
    end function minDist

    ! Always provide square root of r^2
    function calcForce(xi, xj, r) result(f)
      implicit none
      real*8, intent(in) :: xi, xj,r
      real*8::f, sig6, xDiff
      xDiff = minDist(xi, xj)
      sig6 = (sigma**6)/(r**6)
      f = 6.0d0*sig6*eps*(2.0d0*sig6-1.0d0)*(xDiff)/(r**2) - f_rc*(xDiff)/abs(r)
    end function calcForce

    function calcKE(vx, vy, vz) result(KE)
      implicit none
      real*8, intent(in):: vx, vy, vz
      real*8 ::KE
      KE = 0.50d0 * mass * (vx**2 + vy**2 + vz**2)
    end function calcKE

    function calcPE(r) result(PE)
      implicit none
      real*8, intent(in) :: r
      real*8 :: PE, sig6
      sig6 = (sigma**6)/(r**6)
      PE = eps*sig6*(sig6 - 1.0d0) - u_rc + r * f_rc -rc * f_rc
    end function calcPE

		
		subroutine posInit() 
			implicit none
			integer :: i
			open(21, file = 'x.dat')
			open(22, file = 'y.dat')
			open(23, file = 'z.dat')
		
			do i = 1, particles
				read(21,*) x(i)
			end do 
		
			do i = 1, particles
				read(22,*) y(i)
			end do 
		
			do i = 1, particles
				read(23,*) z(i)
			end do 
		
		
			close(21)
			close(22)
			close(23)
		end subroutine posInit

		subroutine updatePos()
      implicit none
      integer :: i,j
      real*8::r

        do i = 1, particles

            x(i) = x(i) + velocity_x(i) * dt + force_x(i) * dt2by2
            y(i) = y(i) + velocity_y(i) * dt + force_y(i) * dt2by2
            z(i) = z(i) + velocity_z(i) * dt + force_z(i) * dt2by2
            !! Boundary cross condition
            x(i) = MODULO(x(i), length)
            y(i) = MODULO(y(i), length)
            z(i) = MODULO(z(i), length)
        end do
    end subroutine updatePos


		subroutine initForce()

			implicit none
			integer :: i
			open(41, file = 'fx.dat')
			open(42, file = 'fy.dat')
			open(43, file = 'fz.dat')
			
      do i = 0, particles
        newForce_x(i) = 0.0d0
        newForce_y(i) = 0.0d0
        newForce_z(i) = 0.0d0
      end do

		
			do i = 1, particles
				read(41,*) force_x(i)
			end do 
		
			do i = 1, particles
				read(42,*) force_y(i)
			end do 
		
			do i = 1, particles
				read(43,*) force_z(i)
			end do 
		
			close(41)
			close(42)
			close(43)
		end subroutine initForce
	
		subroutine updateForce()
      implicit none
      integer :: i,j
      real*8 :: r, fx, fy, fz
      integer :: counter = 0

      potential = 0.0d0

      do i = 0, particles
        newForce_x(i) = 0.0d0
        newForce_y(i) = 0.0d0
        newForce_z(i) = 0.0d0
      end do
      potential = 0.0d0
      do i = 1, (particles-1)
        do j = (i+1), particles
          r = dsqrt(minDist(x(i),x(j))**2 + minDist(y(i), y(j))**2 + minDist(z(i), z(j))**2)

          if(r < rc) then
            !force on ith particle due to all other particles

            fx = calcForce(x(i), x(j), r)
            fy = calcForce(y(i), y(j), r)
            fz = calcForce(z(i), z(j), r)

            newForce_x(i) = newForce_x(i) + fx
            newForce_y(i) = newForce_y(i) + fy
            newForce_z(i) = newForce_z(i) + fz

            !force on jth particle due to ith
            newForce_x(j) = newForce_x(j) - fx
            newForce_y(j) = newForce_y(j) - fy
            newForce_z(j) = newForce_z(j) - fz

            potential = potential + calcPE(r)

          end if
        end do
      end do
    end subroutine

	
		subroutine initVelocities()
			implicit none
			integer :: i
			open(31, file = 'vx.dat')
			open(32, file = 'vy.dat')
			open(33, file = 'vz.dat')
			open(71, file = 'KE_initial.dat')
		
			do i = 1, particles
				read(31,*) velocity_x(i)
			end do 
		
			do i = 1, particles
				read(32,*) velocity_y(i)
			end do 
		
			do i = 1, particles
				read(33,*) velocity_z(i)
			end do 
		
			read(71, *) KE_initial
		
		
			close(31)
			close(32)
			close(33)
			close(71)
		end subroutine initVelocities

    subroutine updateVelocity()
      implicit none
      integer :: i,j
      real*8 :: scalingFactor
      kineticEnergy = 0.0d0

      do i = 1, particles
        velocity_x(i) = velocity_x(i) + (force_x(i) + newForce_x(i)) * dtby2         !Heart of the subroutine. Mass is incorporated in dt2by2
        velocity_y(i) = velocity_y(i) + (force_y(i) + newForce_y(i)) * dtby2
        velocity_z(i) = velocity_z(i) + (force_z(i) + newForce_z(i)) * dtby2
        kineticEnergy = kineticEnergy + calcKE(velocity_x(i), velocity_y(i), velocity_z(i))
      end do
    !!Scaling velocities according to thermoStat
		
      if(mod(updateCounter, thermoCounter) == 0) then
      scalingFactor = KE_initial/kineticEnergy
      scalingFactor = dsqrt(scalingFactor)
        do i = 1, particles
          velocity_x(i) = velocity_x(i) * scalingFactor
          velocity_y(i) = velocity_y(i) * scalingFactor
          velocity_z(i) = velocity_z(i) * scalingFactor
        end do
      end if

      do i = 1, particles
        force_x(i) = newForce_x(i)
        force_y(i) = newForce_y(i)
        force_z(i) = newForce_z(i)
      end do
    end subroutine updateVelocity

		subroutine verletAlgorithm()
      implicit none
      call updatePos()
      call updateForce()
      call updateVelocity()
    end subroutine verletAlgorithm

end module

