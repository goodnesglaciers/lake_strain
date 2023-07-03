      subroutine errlps(b,m,mdim,pa,pb,al)
c
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
c  this subroutine computes error ellipses from the variance-covariance matrix
c  taken from main12.for
c  program by WH Prescott
c  modified slightly be P Segall 2/22/88
c
c  modified by P Segall 1/10/88 so that the axes pa and pb are:
c
c                   x^2                 y^2
c           _______________     +    ________________   =  1
c            Chi^2_p*(s_x)^2         Chi^2_p*(s_y)^2
c  and
c      s_y^2 and s_x^2 are the prinicpal components of the covariance
c      matrix and x, y are associated components of displacement.
c      Chi^2_p is the Chi-squared value corresponding to probability 
c      p (with two degrees of freedom).  The semi axes of the 
c      p error ellipse are thus:
c    
c                       a or b = sqrt(2.0 * Chi^2_p) * s_x or s_y
c      with
c            Probability   50%   75%    90%       95%        99%
c     sqrt(2.0 * Chi^2_p) 1.177 1.665  2.146     2.448      3.035
c
c     The routine currently returns 95% error ellipses
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
      real*8 b(mdim,mdim),pa(mdim),pb(mdim),al(mdim)
      do 2 l=2,m,2
      i=l/2
      arg1 = 2. * b(l,l-1)
      arg2 = b(l-1,l-1) - b(l,l)
      if((arg1.ne.0.0).or.(arg2.ne.0.0))  go to 1
      if(b(l-1,l-1).eq.b(l,l)) then
c.......isotropic covariance matrix
c        pa(i) = 2.0*sqrt(b(l,l))
        pa(i) = 2.448*sqrt(b(l,l))
        pb(i) = pa(i)
        g     = 1.570796
        go to 2
      else
c.......zero covariance matrix
        pa(i) = 0.0
        pb(i) = 0.0
        g     = 1.570796
        go to 2
      endif
    1 continue
      g = 0.5 * atan2(arg1,arg2)
      cg = cos(g)**2
      sg = sin(g)**2
      scg = sin(2.*g)
c
      arg3 = b(l-1,l-1) * cg
     +      +b(l,l-1)   * scg
     +      +b(l,l)     * sg
      arg4 = b(l-1,l-1) * sg
     +      -b(l,l-1)   * scg
     +      +b(l,l)     * cg
c
      pa(i) = 0.0
c      if (arg3.ge.0.0) pa(i) = 2. * sqrt(arg3)
      if (arg3.ge.0.0) pa(i) = 2.448 * sqrt(arg3)
      pb(i) = 0.0
c      if (arg4.ge.0.0) pb(i) = 2. * sqrt(arg4)
      if (arg4.ge.0.0) pb(i) = 2.448 * sqrt(arg4)
      if(pa(i).ge.pb(i)) go to 2
      paa = pa(i)
      pa(i) = pb(i)
      pb(i) = paa
      g = g + 1.5707963
    2 al(i) = 90. -g/0.01745329d0


      return
      end
