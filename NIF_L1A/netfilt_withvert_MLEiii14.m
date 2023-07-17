function [val, sig2, x_s, xcov] = ...
	  netfilt_withvert_MLEiii14(ts_enu, Kern, x0, var0, sitesV, X0_enu, V0_enu, ...
      V0_cov_enu, ref_epoch, origin, R1, R2, smooth, c, sig_f,...
      theta2C, theta3C, theta2B, theta3B, theta1)


%       function [x_s, xcov, outrsum, N, outvsum, val, sig2] = ...
% 	  netfilt_withvert_MLE(ts_enu, Kern, x0, var0, sitesV, X0_enu, V0_enu, ...
%       V0_cov_enu, ref_epoch, origin, R1, R2, smooth, c, sig_f,...
%       tau, sigma, alphaC, alphaB, gammaC, gammaB)
 
	  
%%
%%    function [x_kgk, x_kgn, rsum, N, vsum] = ...
%%        netfilt_sparse(sites1, sites2, files1, files2, G, G2, ...
%%               sigmas, sig_f, X0, X0_enu_apr, V0_enu_apr, V0_enu_cov, ...
%%		 x0, var0, origin, R1, R2, ref_epoch, eq_epoch, smooth)
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
%%  modified 1/14/98 to implement smoothing at each step of iteration
%%  modified 6/23/98 to change vector storage of covariance for version 5
%%  modified 10/15/98 to make first basisfunction uniform slip
%%  9/2013   attempt to undo lots of things back to original nif but keep gpsts and state objects JM

%-(1). Determine some dimensions ------------------------------
        Nepochs = size(ts_enu.epochs,1);
        n = Nepochs;
        [Nsites, Nbasis]  = size(Kern);
	    Nframe = length(sig_f);
        statedim = 2*Nbasis + Nsites + Nframe;

%-(2). Variance parameters ------------------------------------
%        hypers=[tau,sigma,alpha,gamma];
    sigma=2;
    tau=theta1*sigma;
    alpha1=theta2C*sigma;
    alpha2=theta2B*sigma;
    alpha3=alpha2;
    gamma1=theta3C*sigma;
    gamma2=theta3B*sigma;
    gamma3=gamma2;
    
%(3).	PETER's state-variable object
   	x_s=state([zeros(statedim,1)]);
    xcov=state(zeros(statedim,statedim),'swap14/');

%-(4). aprior values for state and covariance matrix ----------

    x_s(1,0) = x0';
	xcov(1,0)=var0;
	x=x0';
    covx = sparse(var0);
	sqrsigma = chol(covx)';
    vsum = 0; rsum = 0;


%-(5). First update step to get x_1|1 and covariance,
%      using straight Kalman filter ---------------------------

    %-(5.2). create data vector "d_t", its covariance "cov_t",
    %        and design matrix "Hd"
    %  note all H matricies should include the whole statevector RW and IRW
        
	pos=ts_enu(:,1).d;
	pos_cov=ts_enu(:,1).dcov;
	
	%%%[pos2, pos_cov2]=rdim(pos,pos_cov);   %for ditching vertical comps
	sitecodes=char(ts_enu(:,1).sites);
	sites=char(ts_enu(:,:).sites);
	disp('Hello 1')
    %[d_t, cov_t, Hd, hd] = ...
    %   makedandH_sparse5_ext_novert(pos2',pos_cov2, sitecodes, sites, G, X0_enu, V0_enu, V0_cov_enu, ...
	%	origin, ref_epoch, ts_enu.epochs(1),  alpha, sigma, tau, x, c2);
    [d_t, cov_t, Hd ] = ...
        makedandH_withvert(pos',pos_cov, sitecodes, sites, sitesV, Kern,  X0_enu, V0_enu, V0_cov_enu, ...
		origin, ref_epoch, ts_enu.epochs(1));
    
    N = length(d_t);

    %whos d_t; keyboard;
    %-(5.3). add pseduo data for spatial smoothing
	R3=R2;   % both are basal
    Nsubf1=size(R1,1); Nsubf2=size(R2,1); Nsubf3=size(R3,1);
	[d_pseudo, Hs, Rsmooth]  = smoothH_vel_sparse(R1, R2,R3, Nbasis, statedim, gamma1, gamma2, gamma3);
    Nsmooth = length(d_pseudo);
    %Nsmooth=0;
	
    % NOW DO POSITIVITY PSEUDO DATA
    % [d_pseudop, Hp, hp]= positive_sparse_ext(x,alpha,rho,Nbasis,statedim);
    %[d_pseudop, Hp, hp]= positive_thrust_neg_ss(x,alpha,rho,Nbasis,statedim);
    %[d_pseudop, Hp, hp]= slip_angle_ekf(x,alpha,rho,Nbasis,statedim,-90,45);
    %Nsmooth = Nsmooth + length(d_pseudop);
    %-(5.4). construct full data vector and its covariance
	%whos Hd Hs Hp hd hs hp d_t d_pseudo d_pseudop
	%disp('Hello3')
	%keyboard
	%d = [d_t; d_pseudo; d_pseudop];
	%H = [Hd; Hs; Hp];
	%h = [hd; hs; hp];
    
    % concatinate real and pseudo data
    d= [d_t; d_pseudo];
    %H= [Hd;];  % not needed for straight filter
    H= [Hd; Hs];
    %-(5.3). "Inovation" or prediction error 
    %whos d H x
    %keyboard
    nu = d - H*x; 

    %-(5.4). scaling of data covariance "cov_t" by the hyperparameter "sigma"
    

        R = sparse([sigma*cov_t,     zeros(N,Nsmooth); ...
             zeros(Nsmooth, N), Rsmooth]);
         whos R
     %    keyboard
	
    %-(5.5). kalman gain  %NOT USING SMOOTHING DATA FOR g??????

	ud = Hd*sqrsigma;
	Vd = R(1:N,1:N) + ud*ud';
	invV = inv(Vd);
    u = H*sqrsigma;
    V = R + u*u';
    invH = inv(V);
    g = covx*H'*invH;    % Kalman gain
    
    %-(5.6). updated state x_1|1 
        x = x + g*nu;

    %-(5.7). updated state state covariance C_1|1 
	[A, p] = chol(R);
	if p~=0
		Nd = length(d)-p+1;
		CR = [A, zeros(p-1,Nd);zeros(Nd, p-1), zeros(Nd,Nd)];
	else
		CR = A;
    end
	s = [ CR, zeros(N+Nsmooth,statedim); u', sqrsigma'];
	r1 = qr(s);
	sqrsigma = r1(N+Nsmooth+1:N+Nsmooth+statedim, ...
                        N+Nsmooth+1:N+Nsmooth+statedim)';
	covx = sqrsigma * sqrsigma';

    %-(5.8). store x_1|1 and its covariance into
    %        "x_kgk(1,:)" and "sigma_kgk(1,:)".
	x_s(1,1) = x;
	xcov(1,1)=full(covx);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%MLE NEW STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %-(5.8). components of likelihood
    %        TAKE OUT ONLY COMPONENTS ASSOCIATED WITH REAL DATA;

        nu1 = nu(1:N);
        rsum = rsum + nu1'*invV*nu1;
% change to avoid ZERO determinant
        EV = eig(full(Vd));
        EVsum = sum(log(EV));
        vsum = vsum + EVsum;
%        vsum = vsum + log(  det( V(1:N,1:N) )  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       
   %keyboard
    %-(5.10). clear some variable for the next use
    clear d_t  cov_t pos  pos_cov   sitecodes


%-(6). Run Square Root Filter ----------------------------------

	%n=160;
	%keyboard
  for k = 2:n
	k
    %-(6.1). Prediction step
       %-(6.1.1.5). store epoch in t(k) for later use
       t(k) = ts_enu.epochs(k);   % store epoch in t(k) for later use

       %-(6.1.1.6). observation at this epoch

       %-(6.1.2). define the epoch interval
       delt = ts_enu.epochs(k) - ts_enu.epochs(k-1);
       %delt=delt*(1/365);  % time units days to years to match everything else?

      %-(6.1.3). construct state transition matrix "F" and
      %          process noise matrix "Q"
      % keyboard
        F = makeF_sparse(statedim, Nbasis, Nframe, delt);
        Q = makeQ_sparse(statedim, Nbasis, Nframe, Nsubf1, Nsubf2, Nsubf3, delt, alpha1, alpha2, alpha3, tau, sig_f);
		
	  %ANTENNA CHANGES
	  if(ts_enu.epochs(k) >=2006.020 & ts_enu.epochs(k) <= 2006.026)
	   disp('ANTENNA CHANGE!!!!!')
	   kk = GetIndex(char(ts_enu.sites), 'P506');
       if(kk>=1)
           %keyboard
	    ie=Nbasis*3+kk*2-1;
	    in=Nbasis*3+kk*2;
	    Q(ie,ie)=5e6;;
	    Q(in,in)=5e6;
       end
	  end


      %-(6.1.4). predicted state "x_k+1|k" and its covariance "C_k+1|k"

        x = F*x;
        covx =  F*covx*F' + Q;	
	    %x_s(k,k-1) = x;      %just for writing out to disk
	    %xcov(k,k-1)=covx;    % ditto
       
       %whos Q
       %keyboard
      %-(6.1.5). Cholesky decomposition of "C_k+1|k"


      %-(6.2). Updated step
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%MLE NEW STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   	 tau = x(statedim-4);
% 		 sigma=x(statedim-3); 
% 	         alpha=x(statedim-2);
% 	         gamma=x(statedim-1);
% 	         rho=  x(statedim);
% 		 tau
% 		 sigma
% 		 alpha
% 		 gamma
% 		 rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %-(6.2.1). create data vector "d_t", its covariance "cov_t",
      %          and design matrix "H"
    	pos=ts_enu(:,k).d;
    	pos_cov=ts_enu(:,k).dcov;
        pos_cov=.05^2*eye(length(pos));
        
	    %[pos2, pos_cov2]=rdim(pos,pos_cov);
    	sitecodes=char(ts_enu(:,k).sites);
	    sites=char(ts_enu(:,:).sites);
        [d_t, cov_t, Hd ] = ...
        makedandH_withvert(pos',pos_cov, sitecodes, sites, sitesV, Kern,  X0_enu, V0_enu, V0_cov_enu, ...
		origin, ref_epoch, t(k));
     %keyboard
        
      %-(6.2.2). dimension of data vector "d_t"
        Nobs = length(d_t);

      %-(6.2.3). add pseduo data for smoothing
      [d_pseudo, Hs, Rsmooth]  = smoothH_vel_sparse(R1, R2,R3, Nbasis, statedim, gamma1, gamma2, gamma3);
      Nsmooth = length(d_pseudo);
 
      %-(6.2.4). dimension of data vector "d_t"
	  %  Nsmooth = length(d_pseudo);
      %Nsmooth=0;
      
 
    % NOW DO POSITIVITY PSEUDO DATA
    %   [d_pseudop, Hp, hp]= positive_sparse_ext(x,alpha,rho,Nbasis,statedim);
    %	Note for Cascadia: Average pure thrust direction is ~ -125 to -150, 
	%	average subduction direction is ~ -125. So strike-slip comp should
	%	be > 0, ie right lateral.
       %[d_pseudop, Hp, hp]= positive_thrust_neg_ss(x,alpha,rho,Nbasis,statedim);
       %[d_pseudop, Hp, hp]= slip_angle_ekf(x,alpha,rho,Nbasis,statedim,-90,45);
    %   Nsmooth = Nsmooth + length(d_pseudop);
    %-(5.4). construct full data vector and its covariance
	%d = [d_t; d_pseudo; d_pseudop];
	%H = [Hd; Hs; Hp];
	%h = [hd; hs; hp];
	%cov_t=cov_t*1e1;
    
      % concatinate real and pseudo data
      d= [d_t; d_pseudo];
      %H= [Hd;];  % not needed for straight filter
       H= [Hd; Hs];
       %-(5.3). "Inovation" or prediction error 
       %whos d H x
       %keyboard
      nu = d - H*x; 
      %-(5.4). scaling of data covariance "cov_t" by the hyperparameter "sigma"
        R = sparse([sigma*cov_t,     zeros(Nobs,Nsmooth); ...
             zeros(Nsmooth, Nobs), Rsmooth]);

        nu=d-H*x;
	    [W,p]=chol(covx);
    	if p~=0
            warning('covx no longer positive definite.')
        end
	    v=W*H';
        K=covx*H'*inv(v'*v + R);
         %Update the state and covariance
	     %keyboard
          x=x + K*nu;
          r=chol([R + v'*v,v'*W ;W*v,W'*W]);
          C=r(end-statedim+1:end,end-statedim+1:end);
          covx=C'*C;
          [dsave(k),d2save(k)]=max(d);
          x_s(k,k) = x;
	      xcov(k,k)=full(covx); 
	      clear sigma_kgk; sigma_kgk = zeros(1,statedim^2);  
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%MLE NEW STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
%     -(6.2.9). Components of Likelihood

% % %         nu1 = nu(1:Nobs);
% % %         rsum = rsum + nu1'*invV*nu1;
% % %         EV = eig(full(Vd));
% % %         EVsum = sum(log(EV));
% % %         vsum = vsum + EVsum;
% % %         N = N + Nobs;
        
        
%        TAKE OUT ONLY COMPONENTS ASSOCIATED WITH REAL DATA;

        nu1 = nu(1:N);
        rsum = rsum + nu1'*invV*nu1;
% change to avoid ZERO determinant
        EV = eig(full(Vd));
        EVsum = sum(log(EV));
        vsum = vsum + EVsum;
%        vsum = vsum + log(  det( V(1:N,1:N) )  );
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %-(6.2.10). clear some variables for the next use

        clear d_ts cov_t pos  pos_cov   sitecodes
  end;
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%MLE NEW STUFF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %-(17). MLE value and its index

        outrsum = N*log(rsum);
        outvsum = vsum;
        val = N*log(rsum) + vsum;
        sig2 = rsum/N;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keyboard
%-(7). condition whether we run the smoother
%      if "1", continue, else, terminate this
%      routine here.
disp('smoothing')
%keyboard
  if smooth ==1;
  disp('Running Smoother')

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
	%x still is x_NgN
	%sigma_kgn=xcov(end,end);
	

        clear sigma_kgk; sigma_kgk = zeros(1,statedim^2);

        clear sigma_kgn; sigma_kgn = zeros(1,statedim^2);

  %-(7.2). initialize state covariance matrix

 	clear s, s = zeros(statedim,statedim);
	sigkgk = eye(statedim);
	sigkp1gN = eye(statedim);

  %-(7.3). backward smoothing

    for k = (n-1):-1:1
	k
    %-(7.3.1). define the epoch interval

        delt = ts_enu.epochs(k+1) - ts_enu.epochs(k);

    %-(7.3.2). construct the state transition matrix "F" and
    %          process noise covariance "Q"


        F = makeF_sparse(statedim, Nbasis, Nframe, delt);
        %%Q = makeQ_sparse(statedim, Nbasis, Nframe, ...
		%%delt, alpha, tau, sig_f);
       
        Q = makeQ_sparse(statedim, Nbasis, Nframe, Nsubf1, Nsubf2, Nsubf3,  delt, alpha1, alpha2, alpha3, tau, sig_f);

    
	if(ts_enu.epochs(k) >=2006.018 & ts_enu.epochs(k) <= 2006.026)
	  kk = GetIndex(char(ts_enu.sites), 'P506');
          if(kk>=1)
	   ie=Nbasis*3+kk*2-1;
	   in=Nbasis*3+kk*2;
	   Q(ie,ie)=5e6;;
	   Q(in,in)=5e6;
          end
	end   
       %saveQb(k)=Q(Nbasis*3+11*2-1,Nbasis*3+11*2-1);



	%-(7.3.2). load state covariance

         %x_s(k,k) = x;
	sigkgk=xcov(k,k); 
    clear sigma_kgk; sigma_kgk = zeros(1,statedim^2);

    %-(7.3.3). smoothed state vector "x_k|n" and covariance "C_k|n"

	sigkp1gk =  F*sigkgk*F' + Q;
    s = sigkgk*F'*inv(sigkp1gk);
	x_k1gk = F*x_s(k,k); 
    %keyboard
 	x = x_s(k,k) + s*(x_s(k+1,end) -  x_k1gk);
    %keyboard

        sigma_kgn=xcov(k+1,end);  % load sigma_kgn(k+1,:) as sigma_kgn
        sigkp1gN(:) = sigma_kgn;
        clear sigma_kgn; sigma_kgn = zeros(1,statedim^2);

 	sigkgN     = sigkgk + s*(sigkp1gN -  sigkp1gk)*s';
	sigma_kgn = full(sigkgN(:)');

        
	xcov(k,end)=sigma_kgn;
	x_s(k,end)=x;
	
        clear sigma_kgn; sigma_kgn = zeros(1,statedim^2);

    end;
  end;

clearvars -except val sig2 x_s xcov

%save Qs saveQf saveQb
