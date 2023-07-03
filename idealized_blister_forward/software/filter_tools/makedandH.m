    function [d_t, cov_t, H] = makedandH(sitecodes, sites, pos, pos_cov, ...
		tk, G, v_pre, origin)
%    function [d_t, cov_t, H] = makedandH(sitecodes, sites, pos, pos_cov, ...
%		 tk, G, v_pre, origin)
%
% input: 
%		sitecodes   = list of stations observed this epoch
%		sites	    = list of all stations
%		pos	    = positions observed at this epoch
%		pos_cov	    = position covariance at this epoch
%		tk	    = time of this epoch
%		G	    = matrix relating displacement to slip
%		v_pre	    = predicted secular velocity
% output: 
%		d_t 	    = data vector for this epoch
%		cov_t	    = data covariance for this epoch
%		H	    = matrix relating data to state vector at
%				curent epoch


[Nstsob, four] = size(sitecodes);
[Ntotal, fur]  = size(sites);
[np, Nbasis] = size(G);


% for each station find matching station and extract 
% appropriate part of G into Gatt. also form matrix Ir

	Gatt = [];
	Ir = zeros(3*Nstsob,3*Ntotal);

for i = 1:Nstsob
	k = GetIndex(sites, sitecodes(i,:));
	Gatt = [Gatt; G(3*(k-1)+1:3*(k-1)+3,:)];
	Ir(3*(i-1)+1:3*(i-1)+3, 3*(k-1)+1:3*(k-1)+3) = eye(3);
	sec_rate(3*(i-1)+1:3*(i-1)+3) = v_pre(:,k);
end



%% rotate the positions and covariance into ENU

	[ENU, covENU] = xyz2enu(pos(:), pos_cov, origin);

%% correct for secular deformation
%% NOTE we are not correcting the covariance for secular deformation

	ENU = ENU - sec_rate'*tk;

%% Now form data & covariance, and kernel as relative positions

	cv = [];
	for i=1:Nstsob-1
		cv = [cv; -eye(3)];
	end
	T = [cv, eye(3*(Nstsob-1))];

	d_t = T*ENU;
	cov_t = T*covENU*T';
	G_t = T*Gatt;
	Ir = T*Ir;


%%  construct the submatrix [tk,1,0]

        sub = zeros(Nbasis, 3*Nbasis);
        for i = 1:Nbasis
                sub(i,3*(i-1)+1:3*(i-1)+3) = [tk,1,0];
        end


%% Finally form H
        H = [G_t*sub, Ir, Ir];

