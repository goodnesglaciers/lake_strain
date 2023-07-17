  function H = makeH(Nsites, Nbasis, tk, G)
%%
%%  H = makeH(Nsites, Nbasis, tk, G)
%%
%%  subroutine to make the H matrix for Kalman filter
%%  H is (N x 3Nbasis+Nsites)
%%
%%  Input:
%%	Nsites   = number of observation sites
%%	Nbasis	 = number of basis functions
%%	tk       = t_k time at present iteration
%%	G        = Matrix that maps slip to data
%%


%% first construct the submatrix [tk,1,0]

	sub = zeros(Nbasis, 3*Nbasis);
	for i = 1:Nbasis
		sub(i,3*(i-1)+1:3*(i-1)+3) = [tk,1,0];
	end
	
%% form H
	H = [G*sub, eye(Nsites)];
