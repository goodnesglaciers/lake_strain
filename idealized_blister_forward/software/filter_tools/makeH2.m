    function H = makeH2(Index, Nbasis, tk, G)
%%
%%  H = makeH(Index, Nbasis, tk, G)
%%
%%  subroutine to make the H matrix for Kalman filter
%%  Paul Segall, Stanford University
%%
%%  Modified for irregularly sampled data
%%  7/3/97
%%
%%  H has dimension (sum(Index) x 3*Nbasis+2*Nsites)
%%
%%  Input:
%%	Index    = Vector showing which sites occupied at epoch
%%	Nbasis	 = number of basis functions
%%	tk       = t_k time of present eopch
%%	G        = Matrix that maps slip to data
%%


%% first construct the submatrix [tk,1,0]

	sub = zeros(Nbasis, 3*Nbasis);
	for i = 1:Nbasis
		sub(i,3*(i-1)+1:3*(i-1)+3) = [tk,1,0];
	end
	

%% now form the matrix G appropriate for the present epoch

N = sum(Index);
Nsites = length(Index);

Ir = zeros(N,Nsites);

count = 1;
for i = 1:Nsites
   if Index(i) == 1;
      Gatt(count,:) = G(i,:);
      Ir(count,i) = 1;
      count = count + 1;
   end
end



%% form H
	H = [Gatt*sub, Ir, Ir];
