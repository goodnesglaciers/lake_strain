    function F = makeF_sparse(statedim, Nbasis, Nframe, delt)
%%
%%  F = makeF_f(statedim, Nbasis, Nframe, delt)
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
		index = 2*(i-1)+1;
		F(index,index+1) = delt;
    end

    % L its the same so leave as eye
    
    % for frame it's zero mean
	for i = statedim-Nframe:statedim   
		F(i,i) = 0;
	end

	F = sparse(F);
