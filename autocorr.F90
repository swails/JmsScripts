! Calculation of the autocorrelation function written by Adrian Roitberg

program autocorr

   integer, parameter :: MAXVALS = 1000000
   double precision, dimension(MAXVALS) :: energy
   double precision, dimension(MAXVALS) :: acoef
   double precision :: mu, denom, junk
   integer :: ntot

   acoef(:) = 0.d0
   energy(:) = 0.d0

   do i=1, MAXVALS
      read(*,*,end=10) time, energy(i)
   end do

   read(*,*,end=10) time, junk
   write(0,*) 'Warning: More than ', MAXVALS, ' points detected!'
   write(0,*) '         Adjust MAXVALS in autocorr.F90 and recompile.'
   i = i + 1

10 ntot = i - 1


   mu = 0.0
   do i = 1, ntot
      mu = mu + energy(i)
   end do

   mu = mu / dble(ntot)


   denom = 0.0

   do i = 1, ntot
      denom = denom + (energy(i) - mu) * (energy(i) - mu)
   end do

   do i = 0, ntot
      do k = i + 1, ntot
         acoef(i+1) = acoef(i + 1) + (energy(k - i) - mu) * (energy(k) - mu)
      end do
      acoef(i + 1) = acoef(i + 1) / denom
      write(*,*) i, acoef(i + 1)
   end do

end program autocorr
