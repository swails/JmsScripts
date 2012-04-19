      program twodhis 

      implicit none
      real :: x, y, xmin, xmax, ymin, ymax
      real :: xbinsize, ybinsize, xcenter, ycenter
      real, dimension(:,:), allocatable :: Prob
      integer :: n, i, j, numstate, inputstatus, xbinnum, ybinnum

      open (unit=10, file="hist.in")
      open (unit=20, file="hist.out")

      numstate=0
      do n=1, 1
        read(10, *) x, y
        xmin=x
        xmax=x
        ymin=y
        ymax=y
      end do
      rewind(10)
      
      do
        read (10, *, IOSTAT=inputstatus) x, y
        if (inputstatus < 0) EXIT
        numstate=numstate+1
        xmin=min(xmin, x)
        xmax=max(xmax, x)
        ymin=min(ymin, y)
        ymax=max(ymax, y)
      end do
      print *, "number of lines=", numstate
      print *, "x minimum=", xmin
      print *, "x maximum=", xmax
      print *, "y minimum=", ymin
      print *, "y maximum=", ymax
      rewind(10)
      
      xbinsize=0.1
      ybinsize=0.1
      xbinnum=(xmax-xmin)/xbinsize+1
      ybinnum=(ymax-ymin)/ybinsize+1
      print *, "number of bins in x dimension:", xbinnum
      print *, "number of bins in y dimension:", ybinnum
      allocate(prob(xbinnum, ybinnum))
      rewind(10)
      
      do i=1, xbinnum
        do j=1, ybinnum
          Prob(i,j)=0.0
        end do
      end do

      do n=1, numstate
        read(10, *) x, y
        i=(x-xmin)/xbinsize+1
        j=(y-ymin)/ybinsize+1
        Prob(i,j)=Prob(i,j)+1.0
      end do

      do i=1, xbinnum
        do j=1, ybinnum
          xcenter=xmin+i*xbinsize-0.5d0*xbinsize
          ycenter=ymin+j*ybinsize-0.5d0*ybinsize
          Prob(i,j)=Prob(i,j)/numstate
          write(20, *) xcenter, ycenter, Prob(i,j)
        end do
      end do
      
      deallocate(prob)
      close(10)
      close(20)
      end
