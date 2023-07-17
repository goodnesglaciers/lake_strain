  function [x_kgk,sigma_kgk,x_kgn,sigma_kgn,x_kp1gk, sigma_kp1gk, rsum, N, vsum] = ...
	netfilt_l(sites, files, G, sigmas, x0, var0, u_pre, origin, smooth)
%%
%%  [x_kgk,sigma_kgk,x_kgn,sigma_kgn,x_kp1gk, sigma_kp1gk, rsum, N, vsum] = ...
%%		netfilt_l(sites, files, G, sigmas, x0, var0, u_pre, smooth)
%%
%% Square Root filter for estimating network model
%% Data is:     d(t) = G*s(t) + B(t) + white noise
%% slip is: 	s(t) = b*t + W(t)
%%
%% Input :  
%%	t	= vector of observation times 
%%	Data 	= Nsites*Nepochs matrix of observations
%%	G 	= Nsites*Nbasis maps slip to Signal
%% 		Nsites	= number of observation sites
%% 		Nepochs	= number of observation epochs
%%	sigmas = [tau, sigma, alpha]
%% 		tau	= random walk standard deviation, mm/sqrt(yr)
%% 		sigma	= white noise standard deviation, mm
%% 		alpha 	= scale of Weiner process
%%	x0 	= state vector at initial epoch (length statedim)
%%	var0 	= diagonal of state covariance at initial epoch 
%%			(statedim X statedim); statedim = 3*Nbasis + Nsites;
%%	smooth 	= Run smoother if smooth = 1;
%%	
%% Output : 
%%		x_kp1gk 	= predicted state vector, x_{k+1|k}
%%		sigma_kp1gk 	= corresponding covariance
%%		x_kgk		= filtered state vector, x_{k|k}
%%		sigma_kgk	= corresponding covariance
%%		x_kgn		= filtered state vector, x_{k|N}
%%		sigma_kgn	= corresponding covariance
%%
%%  State Vector: x = [b_1, W_1(t), dot{W_1(t)}, b_2, W_2(t), ....
%%			B_1(t), B_2(t),	B_3(t),...B_Nsites(t)]'
%%
%%  Observation eq:  d_k = H_k*x_k + eps_k,          eps_k ~ N(0,sig^2*R_k)
%%  Update equation:  x_k+1 = F_k*x_k + delta_k,   delta_k ~ N(0,sig^2*Q_k)
%%
%%  P. Segall and M. Matthews
%%  Stanford University
%%  modified 10/8/96 by P. Segall
%%  modified 7/3/97 For irregularly sampled data, by P. Segall
%%
%% Determine some dimensions
	[Nepochs, rl] = size(files);
	n = Nepochs;
	[Nsites, Nbasis]  = size(G);
	statedim = 3*Nbasis + 2*Nsites;
	
%% Variance parameters
%%	tau = sigmas(1)/sigmas(2); 
%%	alpha=sigmas(3)/sigmas(2);

%% Variance parameters
	sigma = sigmas(1);
	tau = sigmas(2); 
	alpha = sigmas(3);

%% Matricies to be used in smoothing operation & for Square Root Filter
	s = zeros(statedim+Nsites);
 	x_kp1gk = zeros(n,statedim); sigma_kp1gk = zeros(n,statedim^2);  
	x_kgk = zeros(n,statedim); sigma_kgk = zeros(n,statedim^2);  
	x_kgn = zeros(n,statedim); sigma_kgn = zeros(n,statedim^2);

%% Starting (prior) values for state and covariance matrix:
        x = x0;
	covx = var0;

        x_kp1gk(1,:) = x';              % x_{k+1|k} for smoothing
        sigma_kp1gk(1,:) = covx(:)'; % Sigma_{k+1|k} for smoothing
	vsum = 0; rsum = 0;


%% First update step to get x_1/1 and covariance, using straight Kalman filter

%% READ IN SOME DATA
	[pos, pos_cov, t(1), sitecodes] = read_sinex(files(1,:));
        disp(['read file:   ', files(1,:)])
	[d_t, cov_t, H] = makedandH(sitecodes, sites, pos, pos_cov, ...
		0, G, u_pre, origin);

	N = length(d_t);
        nu = d_t - H*x;           % "Inovation" or prediction error

	if t < 91.4603		%scaling based on JGR results
		R = sigma^2*6*cov_t;
	else
		R = sigma^2*10*cov_t;
	end 

        invH = inv(R + H*covx*H');
        g = covx*H'*invH;    % Kalman gain
        x = x + g*nu;
        covx = covx - g*H*covx;

	clear d_t  cov_t pos  pos_cov   sitecodes
        x_kgk(1,:) = x';                % k_{1|1} for smoothing operation
        sigma_kgk(1,:) = covx(:)';   % Sigma_{1|1} for smoothing operation

	%% Components of Likelihood
	rsum = rsum + nu'*invH*nu;
	vsum = vsum + log(  det(R + H*covx*H')  );

%% Run Square Root Filter
	for k = 2:n
		%% Prediction Step:
		[pos, pos_cov, t(k), sitecodes] = read_sinex(files(k,:));
		disp(['read file:   ', files(k,:)])

		delt = t(k) - t(k-1);
		F = makeF(statedim, Nbasis, delt);  % State Transition Matrix
		Q = makeQ2(statedim, Nbasis, delt, alpha, tau);
		 				    % "Process" Covariance
		x = F*x;
		covx =  F*covx*F' + Q;
		sqrsigma = chol(covx)';
		x_kp1gk(k,:) = x';		% k_{k+1|k} for smoothing 
		sigma_kp1gk(k,:) = covx(:)'; % Sigma_{k+1|k} for smoothing


		%% Update Step:
		[d_t, cov_t, H] = makedandH(sitecodes, sites, pos, pos_cov, ...
			t(k)-t(1), G, u_pre, origin);
		Nobs = length(d_t);
    		nu = d_t - H*x;        % "Inovation" or prediction error
		u = H*sqrsigma;

		if t < 91.4603		%scaling based on JGR results
			R = sigma^2*6*cov_t;
		else
			R = sigma^2*10*cov_t;
		end 
    
                invH = inv(R + u*u');
                g = covx*H'*invH;               % Kalman gain
     	  	x = x + g*nu;
     	   	
		s = [ chol(R)', zeros(Nobs,statedim); u', sqrsigma'];
                [q1,r1] = qr(s);
                sqrsigma = r1(Nobs+1:Nobs+statedim, ...
			Nobs+1:Nobs+statedim)';

                covx = sqrsigma * sqrsigma';

		clear d_t  cov_t pos  pos_cov   sitecodes
		x_kgk(k,:) = x';		% k_{k|k} for smoothing 
		sigma_kgk(k,:) = covx(:)'; 	% Sigma_{k|k} for smoothing

		%% Components of Likelihood
		rsum = rsum + nu'*invH*nu;
		vsum = vsum + log(  det(R + u*u')  );
		N = N + Nobs;


	end;
if smooth ==1;

%% Run Smoother 
%%
%%  x_k|N = x_k|k + s_k*(x_k+1|N - x_k+1|k)
%%  sig_k|N = sig_k|k + s_k*(sig_k+1|N - sig_k+1|k)*s^T
%%
%%   where:
%%   	s_k = sig_k|k  *  F_k^T * inv(sig_k+1|k)
%%      
%%   note that 
%%		x_k+1|k = F_k* x_k|k 
%%	& 	sig_k+1|k = F_k*sig_k|k*F_k^T  + Q_k
%%
%%
 	x_kgn(n,:) = x'; 
 	sigma_kgn(n,:) = sigma_kgk(n,:);

 	clear s, s = zeros(statedim,statedim);
	sigkgk = eye(statedim);
	sigkp1gN = eye(statedim);
 
for k = (n-1):-1:1
	delt = t(k+1) - t(k);
	F = makeF(statedim, Nbasis, delt);  % State Transition Matrix
	Q = makeQ2(statedim, Nbasis, delt, alpha, tau);
		 				    % "Process" Covariance
	sigkgk(:) = sigma_kgk(k,:);
	sigkp1gk =  F*sigkgk*F' + Q;
	s = sigkgk*F'*inv(sigkp1gk);

	x_k1gk = F*x_kgk(k,:)';
 	x = x_kgk(k,:)' + s*(x_kgn(k+1,:)' -  x_k1gk);
	x_kgn(k,:) = x';

	sigkp1gN(:) = sigma_kgn(k+1,:);

 	sigkgN     = sigkgk + s*(sigkp1gN -  sigkp1gk)*s';
	sigma_kgn(k,:) = sigkgN(:)';
end;


end;

	
