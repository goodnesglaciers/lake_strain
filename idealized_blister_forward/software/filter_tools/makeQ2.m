  function Q = makeQ2(statedim, Nbasis, delt, alpha, tau)
%%
%% Q = makeQ2(statedim, Nbasis, delt, alpha, tau)
%%
%%  subroutine to make "Process" Covariance for Kalman filter
%%  Q is (statedim x statedim)
%%
%%  Input:
%%	statedim = dimension of state vector
%%	Nbasis	 = number of basis functions
%%	alpha    = scale parameter for integrated random walk
%%	tau      = scale parameter for  random walk
%%	delt     = t_k - t_{k-1} at present epoch
%%	
%%  Paul Segall, July 1997

	Q = zeros(statedim);
	
	for i = 1:Nbasis
		index = 3*(i-1)+2;
		Q(index,index) = alpha^2*delt^3/3;
		Q(index,index+1) = alpha^2*delt^2/2;
		Q(index+1,index) = Q(index,index+1);		
		Q(index+1,index+1) = alpha^2*delt;		
	end
	
% recall statedim = 3*Nbasis + 2Nsites;
	
	Nsites = (statedim-3*Nbasis)/2;

	for i = 1:Nsites
		index=3*Nbasis+i;
		Q(index, index) = tau^2*delt;
	end
		
		
