      subroutine matrdvel(fname, vel, outcov, outnvec, llh)
      
      implicit none
c
c     Read from a standard displacement file
c
      integer       maxvec
      parameter (maxvec=100)

      integer       nvec
      character*150 fname
      character*8   outstanam(maxvec), stanam(maxvec), refsta
      character*3   crdsys
      real*8        xyz(3,maxvec), vel(3,maxvec-1)
      real*8        outxyz(3,maxvec)
      real*8        llhtmp(3), llh(3,maxvec), outnvec
      real*8        cov(3*maxvec,3*maxvec)
      real*8        outcov(3*(maxvec-1),3*(maxvec-1))

c     ... Local vars
      integer      uin
      integer      i, j, k, m, n, refind
      real*8       tmp
      real*8       tmpvel(3)
c
c     ... Initialization
c
      refsta = ' '
      uin = 10
c
c     ... Open velocity file
c
      open(unit=uin, file=fname, form='formatted', status='old', 
     +     err=900)

c     ... Number of parameters
c
      read(unit=uin, fmt='(X,I4,A8,X,A3)') n, refsta, crdsys
      nvec = n/3
      if (nvec.gt.maxvec) then
	 print *,' ** Number of Vectors exceeds Maxvec ** '
         goto 900
      endif
c
c     ... Station names and coordinates
c      if (refsta.ne.' ' .and. refsta.ne.'*')
c      might need some code for more cases where refsta is undefined
c
      do i = 1, nvec+1
         read(unit=uin,fmt='(A8,3(F19.7))')
     +           stanam(i), xyz(1,i), xyz(2,i), xyz(3,i)
	 if ( stanam(i) .eq. refsta ) then
		refind = i
		print *, 'match', refind
        endif
      enddo

c     ... Parameter names and values
c
      m=1
      do i = 1, nvec
         read(unit=uin,fmt='(X,I4,10X,A8,2X,E22.15,6X,E21.15)')
     +       k, stanam(i), tmpvel(1), tmp
         cov(k,k) = tmp**2
c         read(unit=uin,fmt='(X,I4,10X,A8,2X,E22.15,6X,2X,E21.15)')
         read(unit=uin,fmt='(X,I4,10X,A8,2X,E22.15,6X,E21.15)')
     +       k, stanam(i), tmpvel(2), tmp
         cov(k,k) = tmp**2
c         read(unit=uin,fmt='(X,I4,10X,A8,2X,E22.15,6X,2X,E21.15)')
         read(unit=uin,fmt='(X,I4,10X,A8,2X,E22.15,6X,E21.15)')
     +       k, stanam(i), tmpvel(3), tmp
         cov(k,k) = tmp**2
         if(stanam(i).ne.refsta) then
            do j=1,3  
               vel(j,m) = tmpvel(j)
            enddo
            m=m+1
         endif 
      enddo

c
c     ... Now the covariance elements
c
      do i = 1, 3*nvec
         do j = 1, i-1
            read(unit=uin,fmt='(X,I4,2X,I4,2X,E22.15)') i, j, tmp
            cov(i,j) = tmp*dsqrt(cov(i,i)*cov(j,j))
            cov(j,i) = cov(i,j)
         enddo
      enddo
c
c     ...get rid of zeros from covariance matrix
c
 
      do i = 1, 3*nvec
         if (i .lt. ((refind*3)-2) ) then
            do j= 1,i
               outcov(i,j) = cov(i,j)
               outcov(j,i) = outcov(i,j)
            enddo
          else if(i .gt. (refind*3)) then
            m=i-3
            do j= 1,i
               if(j .lt. ((refind*3)-2)) then
                  n=j 
                  outcov(m,n) = cov(i,j)
                  outcov(n,m) = outcov(m,n)
                else if(j .gt. (refind*3)) then
                  n=j-3
                  outcov(m,n) = cov(i,j)
                  outcov(n,m) = outcov(m,n)
                endif
             enddo
           endif 
        enddo


c
c    ... Reshuffle stanam and xyz so that refsta is last
c

      do i = 1, (nvec+1)
         if (i .lt. refind) then
             do j=1,3
                outxyz(j,i) = xyz(j,i)
             enddo
             outstanam(i) = stanam(i)
         else if(i. gt. refind) then
             do j=1,3
                outxyz(j,i-1) = xyz(j,i)
             enddo
             outstanam(i-1) = stanam(i)
         endif
      enddo
      outstanam(nvec+1) = refsta
      do j=1,3
         outxyz(j,nvec+1) = xyz(j, refind)
      enddo

c
c     ... Copying nvec to a double precision # to pass to Matlab
c     
      outnvec = nvec

c
c    ... Convert xyz station coords to lat & lon
c

      do i = 1, (nvec+1)
         call xyz2llh(outxyz(1,i), llhtmp)
         do k =1,3
            llh(k,i) = llhtmp(k)
         enddo 
      enddo
      

      close(unit = uin)

      return

  900 print *,' ** Error in opening file ** '

      end

      subroutine xyz2llh(xyz, llh)
      implicit none
c
      real*8         xyz(3), llh(3)
c
c     ... Local variables
c
      real*8         lat, lon, h
      real*8         a, f, e2, p, r, num, den, n,
     +               deg2rad, oldh
      integer        iter
c
c     ... Code
c
      deg2rad = datan(1.d0)/45.d0
c
c     ... Transform from x,y,z to lat,lon,h
c
      a = 6378137.d0
      f = 1.d0/298.257223563d0
      e2 = 2.d0*f - f*f
      p = dsqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
      r = dsqrt(p*p + xyz(3)*xyz(3))
      lon = datan2(xyz(2),xyz(1))
c
c     ... First iteration on lat and h
c            - assumes h = 0
c
      lat = datan2(xyz(3)/p,(1.d0-e2))
      n = a/sqrt((1.d0 - e2*sin(lat)**2))
      h = p/cos(lat) - n
c
c        ... Iterate until h converges (should be quick since h << n)
c
      oldh = -1.d9
      iter = 0
      num = xyz(3)/p
      do while (abs(h - oldh) .gt. 0.0001)
         iter = iter + 1
         oldh = h
         den = 1.d0 - e2*n/(n+h)
         lat = datan2(num,den)
         n = a/sqrt((1.d0 - e2*sin(lat)**2))
         h = p/cos(lat) - n
      enddo
c      write(stderr,*) 'Converged in ', iter, ' iterations'
c
c     ... Unit conversions
c
      llh(1) = lat/deg2rad
      llh(2) = lon/deg2rad
      llh(3) = h
c
      return
      end
