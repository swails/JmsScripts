      program onedhis 

      implicit none
      real :: frames, distance, mindis, maxdis
      real :: lowerl, upperl, binsize, bincenter
      real, dimension(:), allocatable :: Prob
      integer :: n, i, numstate, inputstatus, binnum, numargs
      character *100 :: inputfile, outputfile

      numargs=iargc()
      if (numargs .eq. 0 .or. numargs > 2) then
        print *, "Usage: 1dhistogram.x inputfile outputfile"
        stop
      end if
      call getarg(1, inputfile)
      call getarg(2, outputfile)
      open (unit=10, file=inputfile)
      open (unit=20, file=outputfile)

      do n=1, 1
        read(10,*) frames, distance
        mindis=distance
        maxdis=distance
      end do
      rewind(10)
      
      numstate=0
      do
        read (10, *, IOSTAT=inputstatus) frames, distance
        if (inputstatus < 0) EXIT
        numstate=numstate+1
        mindis=min(mindis, distance)
        maxdis=max(maxdis, distance)
      end do
      print *, "number of lines=", numstate
      print *, "min dist=",mindis, "max dist=",maxdis
      rewind(10)
      
      lowerl=mindis-0.5
      upperl=maxdis+0.5
      binsize=0.1
      binnum=(upperl-lowerl)/binsize+1
      allocate(prob(binnum))
      rewind(10)
      
      do i=1, binnum
        Prob(i)=0.0
      enddo

      do n=1, numstate
        read(10,*) frames, distance
        i=(distance-lowerl)/binsize+1
        Prob(i)=Prob(i)+1.0
      enddo

      do i=1, binnum
        bincenter=lowerl+i*binsize-0.5d0*binsize
        Prob(i)=Prob(i)/numstate
        write(20, *) bincenter, Prob(i)
      enddo
      
      deallocate(prob)
      end
