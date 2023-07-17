  function [val, sig2] = netloglik_ir(t, Data, G, theta1, theta2)
%%
%% [val, sig2] = netloglik_ir(t, Data, G, theta1, theta2)
%%
%% Evaluate Loglikelihood for network model
%% Data is      d(t) = G*s(t) + B(t) + white noise
%% slip is: 	s(t) = b*t + W(t)
%%
%% Input :  
%%	t	= vector of observation times 
%%	Data 	= Nsites*Nepochs matrix of observations
%%	G 	= Nsites*Nbasis maps slip to Signal
%% 		Nsites	= number of observation sites
%% 		Nepochs	= number of observation epochs
%%	sigmas; 
%%		theta1 = tau/sigma;
%%		theta2 = alpha/sigma;
%% 		tau	= random walk standard deviation, mm/sqrt(yr)
%% 		sigma	= white noise standard deviation, mm
%% 		alpha 	= scale of Weiner process	
%%
%% Output:      val  = C - 2 * log(L)
%% 		C  (= n - n log n) is a constant 
%% 		L is the profile likelihood  
%%		sig2 = estimated variance
%%
%%  State Vector: x = [b_1, W_1(t), dot{W_1(t)}, b_2, W_2(t), ....
%%			B_1(t), B_2(t),	B_3(t),...B_Nsites(t)]'
%%
%%  Observation eq:  d_k = H_k*x_k + eps_k		eps_k ~ N(0,sig^2*R_k)
%%  Update equation:  x_k+1 = F_k*x_k + delta_k		delta_k ~ N(0,sig^2*Q_k)
%%
%% Paul Segall and Mark Matthews
%% Stanford University
%%
%%  modified 7/9/97 For irregularly sampled data, by P. Segall

%% Determine some dimensions
	[Nsites, Nepochs] = size(Data);
	n = Nepochs;
	[ncheck, Nbasis]  = size(G);
	if ncheck ~= Nsites
		disp('G and Data are inconsistent')
	end
	statedim = 3*Nbasis + 2*Nsites;

	t = t - min(t);
	
%% Some parameters
	Big = 10000;
	Small = 1.0e-4;
	tau = theta1; alpha=theta2; 
	vsum = 0; rsum = 0;

%% Matricies to be used in Square Root Filter
	s = zeros(statedim+Nsites);

%% Starting (prior) values for state and covariance matrix:
        x = [zeros(3*Nbasis + Nsites,1); first_obs(Data)'];
	covx = Small*eye(statedim);
	for j=1:Nbasis
		covx(3*(j-1)+1,3*(j-1)+1) = Big;
	end
        for j=1:Nsites
                covx(statedim-Nsites+j, statedim-Nsites+j) = Big;
        end

%% First update step to get x_1/1 and covariance, using straight Kalman filter
 	[datt, index] = actualdata(Data(:,1));
	N = sum(index);
 
 	H = makeH2(index, Nbasis, t(1), G);
        nu = datt - H*x;           % "Inovation" or prediction error
 	R = makeR(index);
	invH = inv(R + H*covx*H');
 	clear datt

        g = covx*H'*invH;    % Kalman gain
        x = x + g*nu;
        covx = covx - g*H*covx;

	%% Components of Likelihood
	rsum = rsum + nu'*invH*nu;
	vsum = vsum + log(  det(R + H*covx*H')  );

%% Run Square Root Filter
	for k = 2:n
		%% Prediction Step:
		delt = t(k) - t(k-1);
		F = makeF(statedim, Nbasis, delt);  % State Transition Matrix
 		Q = makeQ2(statedim, Nbasis, delt, alpha, tau);
		 				    % "Process" Covariance
		x = F*x;
		covx =  F*covx*F' + Q;
		sqrsigma = chol(covx)';

		%% Update Step:
 		[datt, index] = actualdata(Data(:,k));  
 		Nobs = sum(index);
  		H = makeH2(index, Nbasis, t(k), G);
     		nu = datt - H*x;        % "Inovation" or prediction error
 		clear datt

		u = H*sqrsigma;
		R = makeR(index);	%  Normalized data covariance

		invH = inv(R + u*u');
		g = covx*H'*invH;		% Kalman gain
     	  	x = x + g*nu;
     	   	
 		s = [ chol(R)', zeros(Nobs,statedim); u', sqrsigma'];
                [q1,r1] = qr(s);
                sqrsigma = r1(Nobs+1:Nobs+statedim, ...
 			Nobs+1:Nobs+statedim)';
                covx = sqrsigma * sqrsigma';

		%% Components of Likelihood
		rsum = rsum + nu'*invH*nu;
		vsum = vsum + log(  det(R + u*u')  );
		N = N + Nobs;
	end;
%%
val = N*log(rsum) + vsum;
sig2 = rsum/N;
