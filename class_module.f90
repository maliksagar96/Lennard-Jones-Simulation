module class_module
  implicit none

    real*8 :: length = 15.0d0                                         !Box is cubical box of Boxsize = boxsize*boxsize*boxsize
    integer, parameter :: particles = 600
    real, parameter :: dt = 0.0050d0/2.0d0
    real*8, parameter :: sigma = 1.0d0, eps = 1.0d0, mass = 1.0d0, temp = 1.0d0
    real*8, parameter :: rc = 2.50d0*sigma                            !Cutoff Distance
    real*8 :: velocity_x(particles), velocity_y(particles), velocity_z(particles)
    real*8 :: x(particles), y(particles), z(particles)
    integer, parameter :: seed_x = 2550, seed_y = 456, seed_z = 588 !To check random number calling for differernt seeds
    real*8 :: dt2by2 = (dt**2.0d0)/(mass*2.0d0), dtby2 = dt/(2.0d0*mass)
    real*8 :: potential = 0.0d0, kineticEnergy = 0.0d0
    real*8 :: force_x(particles), force_y(particles), force_z(particles)
    real*8 :: newForce_x(particles), newForce_y(particles), newForce_z(particles)
    real*8:: f_rc = eps*((12.0d0*(sigma**12)/(rc**13)) -(6.0d0*(sigma**6)/rc**7))
    real*8 :: u_rc = eps*((sigma**12)/(rc**12) - (sigma**6)/(rc**6))


contains
  !Defining all the funtions first then defining all the subroutines

    function velInit(arrayDim, seed) result(v)
      implicit none
      integer, intent(in) :: arrayDim, seed
      real*8, dimension(arrayDim) :: v
      integer :: i, j = 0
      real*8 :: vc = 1.0d0   !velocity constant
      call srand(seed)
      do i = 1, size(v)
        v(i) = vc*(rand() - 0.05d0)
      end do
    end function velInit

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

    ! Initializing position with a equally spaced particles with a constant separation between 2 adjacent molecules in any direction
    subroutine posInit()

      implicit none
      integer :: i, j, k, pcounter, dummyCounter, counter
      real :: sepFactor = 1.20d0
      integer :: maxCount
      sepFactor = sigma*sepFactor
      maxCount = floor(length/sepFactor)
      x(1) = sepFactor
      y(1) = sepFactor
      z(1) = sepFactor
      counter = 0
      do i = 1, particles-1
        counter = counter + 1
        if(counter == maxCount) then
          counter = 0
          x(i+1) = sepFactor
        else
          x(i+1) = x(i) + sepFactor
        end if

      end do

      counter = 0
      do i = 1, particles-1
        counter = counter + 1

        if(counter == maxCount) then
          counter = 0
          y(i+1) = y(i) + sepFactor
        else
          y(i+1) = y(i)
        end if

        if(y(i+1) > sepFactor*maxCount) then
          y(i+1) = sepFactor
        end if
      end do

      counter = 0
      do i = 1, particles-1
        counter = counter + 1

        if(counter == (maxCount**2)) then
          counter = 0
          z(i+1) = z(i) + 2.0d0*sepFactor
        else
          z(i+1) = z(i)
        end if
      end do
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
      !Calculating force on ith element
      integer :: i,j

      do i = 1, particles
        newForce_x(i) = 0.0d0
        newForce_y(i) = 0.0d0
        newForce_z(i) = 0.0d0
        force_x(i) = 0.0d0
        force_y(i) = 0.0d0
        force_z(i) = 0.0d0
      end do
      call updateForce()

      do i = 1, particles
        force_x(i) = newForce_x(i)
        force_y(i) = newForce_y(i)
        force_z(i) = newForce_z(i)
      end do

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
      real*8 :: avgVx, avgVy, avgVz, total_x
      !initializing without normalizing and Making Velocity of COM zero
      velocity_x = velInit(size(velocity_x), seed_x)
      velocity_y = velInit(size(velocity_y), seed_y)
      velocity_z = velInit(size(velocity_z), seed_z)
      !averaging velocities
      do i = 1, size(velocity_x)
        avgVx = velocity_x(i)+avgVx
        avgVy = velocity_y(i)+avgVy
        avgVz = velocity_z(i)+avgVz
      end do

      avgVx = avgVx/dfloat(particles)
      avgVy = avgVy/dfloat(particles)
      avgVz = avgVz/dfloat(particles)

      !Making Velocity of COM of zero
      do i = 1, particles
        velocity_x(i) = avgVx - velocity_x(i)
        velocity_y(i) = avgVy - velocity_y(i)
        velocity_z(i) = avgVz - velocity_z(i)
      end do
    end subroutine initVelocities

    subroutine updateVelocity()
      implicit none
      integer :: i,j

      kineticEnergy = 0.0d0

      do i = 1, particles
        velocity_x(i) = velocity_x(i) + (force_x(i) + newForce_x(i)) * dtby2         !Heart of the subroutine. Mass is incorporated in dt2by2
        velocity_y(i) = velocity_y(i) + (force_y(i) + newForce_y(i)) * dtby2
        velocity_z(i) = velocity_z(i) + (force_z(i) + newForce_z(i)) * dtby2
        kineticEnergy = kineticEnergy + calcKE(velocity_x(i), velocity_y(i), velocity_z(i))
      end do

      do i = 1, particles
        force_x(i) = newForce_x(i)
        force_y(i) = newForce_y(i)
        force_z(i) = newForce_z(i)
      end do
    end subroutine updateVelocity


    subroutine verletAlgorithm()
      call updatePos()
      call updateForce()
      call updateVelocity()
    end subroutine verletAlgorithm


end module class_module
