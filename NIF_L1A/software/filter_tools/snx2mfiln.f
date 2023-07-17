      subroutine snx2mfiln(fsinex, yr, x, cov,osite)
      
* Convert SINEX format to MATLAB mfile format

      integer luin,luout
      data luin,luout/16,17/
      integer i,j
      integer nsites
      parameter (nsites=50)
      character text*80
      character fsinex*80
      integer iyr,doy,sec
      real*8 yr(3*nsites), x(3*nsites), cov(3*nsites,3*nsites)      
      real*8 osite

* Input files
      open(luin,file=fsinex)

* Solution epochs

      do while (text(1:18) .ne. '+SOLUTION/ESTIMATE') 
       read(luin,'(a)',end=99) text
      enddo

      read(luin,'(a)',end=99) text
      read(luin,'(a)',end=99) text

      n=0
      do while (text(1:18) .ne. '-SOLUTION/ESTIMATE')
       n=n+1
       read(text,'(27x,i2,1x,i3,1x,i5)') iyr,doy,sec
       yr(n) = iyr + (doy-1 + sec/86400.0)/366.0

       read(luin,'(a)',end=99) text
      enddo
      
      rewind(luin)

* Solution estimates

      do while (text(1:18) .ne. '+SOLUTION/ESTIMATE') 
       read(luin,'(a)',end=99) text
      enddo

      read(luin,'(a)',end=99) text
      read(luin,'(a)',end=99) text

      n=0
      do while (text(1:18) .ne. '-SOLUTION/ESTIMATE')
       n=n+1 
       read(text,'(47x,e21.15)') x(n)
       read(luin,'(a)',end=99) text
      enddo
      rewind(luin)

* Solution covariance

      do while (text(1:25) .ne. '+SOLUTION/MATRIX_ESTIMATE') 
       read(luin,'(a)',end=99) text
      enddo

      if (text(27:32) .ne. 'L COVA') then
       write(*,*) 'Unknown matrix_estimate'
       stop
      endif

      read(luin,'(a)',end=99) text
      read(luin,'(a)',end=99) text

      do while (text(1:25) .ne. '-SOLUTION/MATRIX_ESTIMATE')
       read(text,'(1x,i5,1x,i5,1x,e21.14)') i,j,cov(i,j)
       if(i.ne.j) then
        cov(j,i)=cov(i,j)
       endif
       if (j.lt.i) then
        read(text,'(35x,e21.14)') cov(i,j+1)
        if(i.ne.j+1) then
        cov(j+1,i)=cov(i,j+1)
        endif
       endif
       if (j+1.lt.i) then
        read(text,'(57x,e21.14)') cov(i,j+2)
        if(i.ne.j+2) then
        cov(j+2,i)=cov(i,j+2)
        endif
       endif
       read(luin,'(a)',end=99) text
      enddo
      rewind(luin)
c
c  Copy the value of nsite to the double precision so to pass to
c  Matlab
c
      osite=n/3.0
c 
  99  close(luin)
      return
      end
