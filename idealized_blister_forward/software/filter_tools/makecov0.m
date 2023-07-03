  function cov_0g0 = makecov0(statedim, Nbasis, gamma, Lambda)
%%
%%  cov_0g0 = makecov0(statedim, Nbasis, gamma2, Lambda)
%%
%%  subroutine to make apriori state covariance matrix
%%
%%  Input:
%%      statedim = dimension of state vector
%%      Nbasis   = number of basis functions
%%      gamma    = prior standard deviation
%%      Lambda   = singular values (in vector form)


	Small = 1.0e-6*gamma^2*min(Lambda);  
	Int   = 1.0e+6*gamma^2*max(Lambda);

%% most components are set near zero
        cov_0g0 = Small*eye(statedim);

%% components corresponding to slip scale with singular value
        for j=1:Nbasis
                cov_0g0(3*(j-1)+1,3*(j-1)+1) = gamma^2*Lambda(j);
        end

%% components corresponding to positions have sigmas of 1 m.

	Ncomps = 0.5*(statedim - 3*Nbasis);

        for j=1:Ncomps
                cov_0g0(statedim-Ncomps+j, statedim-Ncomps+j) = Int;
        end

