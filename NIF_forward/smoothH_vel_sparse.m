	function [d_pseudo, Hs, Rs]  = smoothH_vel_sparse(R1, R2, R3, Nbasis, statedim, gamma1, gamma2, gamma3)

% 		 [d_pseudo, Hs, hs]  = smoothH_vel_sparse_ext(t, R1, R2, Nbasis, statedim, gamma, alpha, x);
%
%  function to implement spatial smoothing as pseduo-observations
%  Note there is some question as to whether we should be smoothing
%  velocity or slip, or whether it in fact matters.
%  For now smooth slip
%
% 	P.Segall 12/1/1997
%	corrected 3/1/1999, PSegall


	A = zeros(Nbasis, statedim);	% Use this to smooth velocity

	d_pseudo = zeros(Nbasis,1);
    for i = 1:Nbasis
      A(i,2*(i-1)+1:2*i) = [0,1];
    end

% Because the basis functions are orthonormal, there is no scaling

	%whos R1 A R2
    Nsubfa=size(R1,2);
    Nsubfb=size(R2,2);
    Nsubfc=size(R3,2);
	pre = [R1, zeros(Nsubfa,Nsubfb), zeros(Nsubfa,Nsubfc); ...
          zeros(Nsubfb,Nsubfa), R2,  zeros(Nsubfb,Nsubfc); ...
          zeros(Nsubfc,Nsubfa), zeros(Nsubfc,Nsubfb), R3;];
        
	Hs = pre*A;
	Hs = sparse(Hs);
	
    Rs = [gamma1^2*eye(Nsubfa), zeros(Nsubfa,Nsubfb), zeros(Nsubfa,Nsubfc); ...
          zeros(Nsubfb,Nsubfa), gamma2^2*eye(Nsubfb), zeros(Nsubfb,Nsubfc); ...
          zeros(Nsubfc,Nsubfa), zeros(Nsubfc,Nsubfb), gamma3^2*eye(Nsubfc); ];
