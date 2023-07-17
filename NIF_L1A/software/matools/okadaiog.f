      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      integer plhs(*), prhs(*)
      integer nlhs, nrhs
c
      integer mxGetM, mxGetN, mxGetPr, mxCreateFull
      integer m, n
      integer alp, xrec, source, uout
c
c check for proper number of arguments
c
      if (nrhs .ne. 3) then
        call mexErrMsgTxt('three input arguments required')
      elseif (nlhs .gt. 1) then
        call mexErrMsgTxt('one output arguments required')
      endif
c
c Check the dimensions of XREC. 
c
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
c
      if ( (m.ne.2) .or. (n.ne.1) ) then
        call mexErrMsgTxt('OKADAIO requires XREC to be 2 x 1')
      endif
c
c Check the dimensions of SOURCE. 
c
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
c
      if ( (m.ne.10) .or. (n.ne.1) ) then
        call mexErrMsgTxt('OKADAIO requires SOURCE to be 10 x 1')
      endif
c
c Create a matrix for return argument
c
      plhs(1) = mxCreateFull(9,1,0)
c
c Dereference arguments to get array pointers 
      uout  = mxGetPr(plhs(1))
      alp   = mxGetPr(prhs(1))
      xrec  = mxGetPr(prhs(2))
      source = mxGetPr(prhs(3))
c
c Call the Computational Routine
      call okadaio(%VAL(alp), %VAL(xrec), %VAL(source), %VAL(uout))
c
      return
      end

