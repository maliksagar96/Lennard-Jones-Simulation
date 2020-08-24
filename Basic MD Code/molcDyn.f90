program main

    use class_module

    implicit none
    integer :: i,j,k,monteCarloSteps = 5000
    real*8 :: pe, TE, ii,jj,kk
    !Initializing everyhing
    call posInit()
    call initForce()
    call initVelocities()

    !open(20, file = 'KE.dat')
    !open(30, file = 'PE.dat')
    !open(40, file = 'TE.dat')
    open(50, file = 'diffT.dat')

    do j = 1, monteCarloSteps
      call verletAlgorithm()

      TE = kineticEnergy + potential
      !if(mod(j,50) == 0) then
        !do i = 1, particles
      !    ii = ii + velocity_x(i)
      !    jj = jj + velocity_y(i)
      !!  end do
        !print *,ii
      !  print *,jj
      !  print *,kk
        !write(20, *), j,kineticEnergy/dfloat(particles)
        !write(30, *), j,potential/dfloat(particles)
        !write(40, *), j, TE/dfloat(particles)
        write(50, *), j, TE/dfloat(particles)
      !end if
    end do

    !close(20)
    !close(30)
    !close(40)
    close(50)
end program main
