module extModule

	implicit none
	
	integer :: inDim = 6000000
	integer, parameter :: inverseInterval = 1000
	integer, dimension(:), allocatable :: speedArray 
	real, dimension(:), allocatable :: dataArray, v
	integer, dimension(:), allocatable ::freqV
	integer :: dataNum = 71
	
	contains
	
	subroutine dataInput() 
		implicit none
		integer :: i
		open(dataNum, file = 'speed3.dat')
		allocate(dataArray(inDim))	
		allocate(speedArray(inDim))
		do i = 1, inDim
			read(dataNum, *) dataArray(i)
		end do 
		speedArray = int(inverseInterval * dataArray)			
		close(dataNum)
	end subroutine dataInput
	
	subroutine calcFreq() 
		implicit none
		integer :: i
		integer :: maxValue
		maxValue = maxval(speedArray)
		maxValue = maxValue + 1 
		
		open(65, file = 'xAxis3.dat')
		open(69, file = 'yAxis3.dat')
		allocate(v(maxValue))
		allocate(freqV(maxValue))
		freqV = 0
		v = 0
		
		do i = 1, (size(v) - 1)
			v(i+1) = v(i) + (1.0/inverseInterval)
		end do 
		
		do i = 1, size(speedArray)
			freqV(speedArray(i) + 1) = freqV(speedArray(i) + 1) + 1
		end do 
		
		do i = 1, size(freqV)
			write(65, *) v(i)
			write(69, *) freqV(i)
		end do 
		
		close(69)
		close(65)
		
	end subroutine calcFreq
	
	end module