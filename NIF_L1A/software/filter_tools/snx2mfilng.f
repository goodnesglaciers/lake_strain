c   snx2mfilng.f
c     Gareway function for snx2mfiln.f
c
c
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      integer plhs(*), prhs(*)
      integer nlhs, nrhs
c
      integer nsites
      parameter(nsites=50)
c
      integer mxCreateFull,mxGetPr
      integer yr, x, cov, osite
      character*80 fsinex
c
c-----
c
c   Check for proper number of arguments
c
      if (nrhs.ne.1) then
        call mexErrMsgTxt('Only one input arguments allowed')
      endif
      if (nlhs.gt.4) then
        call mexErrMsgTxt('Too much arguments')
      endif
c
c   Check that input is a string
c
      if (mxIsString(prhs(1)).ne.1) then
        call mexErrMsgTxt('argument must be a string!')
      end if
      ilen=mxGetN(prhs(1))
      i=mxGetString(prhs(1),fsinex,ilen)
c
c   Create matrices for the return arguments
c
c     epoch
      plhs(1)=mxCreateFull(3,nsites,0)
c     
c     coordinates(x,y,z)
      plhs(2)=mxCreateFull(3,nsites,0)
c     
c     covariance matrix
      plhs(3)=mxCreateFull(3*nsites,3*nsites,0)
c     
c     number of sites
      plhs(4)=mxCreateFull(1,1,0)
c
c   Dereference arguments to get array pointers
c
      yr=mxGetPr(plhs(1))
      x=mxGetPr(plhs(2))
      cov=mxGetPr(plhs(3))
      osite=mxGetPr(plhs(4))
c
c   Call the computational subroutine
c
      call snx2mfiln(fsinex,%val(yr),%val(x),
     1             %val(cov),%val(osite))
c
      return
      end
