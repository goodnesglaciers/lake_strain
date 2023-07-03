  function F = makeF(statedim, Nbasis, delt)
%%
%%  F = makeF(statedim, Nbasis, delt)
%%
%%  subroutine to make state transition matrix for Kalman filter
%%  F is (statedim x statedim)
%%
%%  Input:
%%	statedim = dimension of state vector
%%	Nbasis	 = number of basis functions
%%	deltt    = t_k - t_{k-1} at present iteration
%%
	F = eye(statedim);
	
	for i = 1:Nbasis
		index = 3*(i-1)+2;
		F(index,index+1) = delt;
	end
