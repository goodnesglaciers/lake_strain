    function    [d_t, cov_t, Hd ] = makedandH_withvert(pos,pos_cov,...
        sitecodes, sites, sitesV, Kern,  X0_enu_apr, V0_enu_apr, V0_cov_enu, ...
		origin, ref_epoch, cur_epoch);
%
% input: 
%		sitecodes   = list of stations observed this epoch
%		sites	    = list of all stations
%		pos	    = positions observed at this epoch
%		pos_cov	    = position covariance at this epoch
%		epoch	    = this epoch
%		G	    = matrix relating displacement to slip
% output: 
%		d_t 	    = data vector for this epoch
%		cov_t	    = data covariance for this epoch
%		H	    = matrix relating data to state vector at
%				curent epoch
% dimension of state vector:
%               position    = m
%               velocoty    = mm/yr
%               translation = m
%               rotation    = mas
%               scale       = 10^(-8)

	Nstsob = size(sitecodes,1);
	Ntotal = size(sites,1);
	[Ncomps, Nbasis] = size(Kern);

% for each station find matching station and extract 
% appropriate part of G into Gatt. also form matrix Ir

	Gatt = []; G=zeros(Ncomps,Nbasis*2);
	mu = []; xapr = []; xapr_use = [];
    Rtr1 = []; Rro1 = []; Rsc1 = [];
	X_apr = zeros(3*Nstsob,1);
	Ir = zeros(3*Nstsob,3*Ntotal);
	tmpCv = []; Cv = [];
	%Full Frame Terms we don't need
	%subR = rotmatrix(origin(2), origin(1));
    %subR=subR(1:2,1:3);  % ONLY NEED HORIZONTALS BUT STILL HAVE 3+3+1 frame 
    %for i = 1:Ntotal
	%  mu = [mu; 0, -X0(3,i), X0(2,i); ...
    %              X0(3,i), 0, -X0(1,i); ...
    %              -X0(2,i), X0(1,i), 0];
    % end
    
    % already in ENU
    subR=eye(3);
    % turn Kern into G by adding dummy rows for RW terms
    for i=1:Ncomps
     for j=1:Nbasis
      G(i,j*2-1)=Kern(i,j);  % only measure IRW(j*2-1) not RW(j*2)
     end
    end
    
    for i = 1:Nstsob
	 k = GetIndex(sites, sitecodes(i,:));
	 if(k ~= 0)
	  tmpCv(3*i-1:3*i,:) = V0_cov_enu(3*k-1:3*k,:);
     else
	  %disp('problem2')
	  %i
     end
    end

for i = 1:Nstsob
	k = GetIndex(sites, sitecodes(i,:));
    if(k ~= 0)
	 Cv(:,3*i-2:3*i) = tmpCv(:,3*k-2:3*k);
	 Gatt = [Gatt; G(3*k-2:3*k,:)];
	 Ir(3*i-2:3*i, 3*k-2:3*k) = eye(3);
     %keyboard
	 X_apr(3*i-2:3*i) = X0_enu_apr(k*3-2:3*k)	+ V0_enu_apr(3*k-2:3*k)*(cur_epoch-ref_epoch);
	 Rtr1 = [Rtr1; subR];
	 %Rro1 = [Rro1; subR*mu(3*(k-2)+1:3*(k-2)+3,:)];
	 %Rsc1 = [Rsc1; subR*X0(:,k)];
     else
        %disp('problem')
	  %i
     end
end

%if(Nstsob<=Ntotal-1)
%    keyboard
%end
%if(cur_epoch>=168)
%    cur_epoch
%    keyboard
%end

%% transform the unit of rotation (rad -> mas)
%%        Rro1 = (pi/18/36)*10^(-8)*Rro1;

%% Transform the unit of scale ( 1 -> 10^(-8) )
%%        Rsc1 = 10^(-8)*Rsc1;

%% rotate the positions and covariance into ENU
	%[ENU, covENU] = xyz2enu(pos(:), pos_cov, origin);

    ENU=pos;
	covENU=pos_cov;

%% subtract the apriori value of the coordinate at each epoch 
%% from data.  removes secular.
%% Should comment out this part if we aren't using them because it modifies the covariance.

	%whos ENU X_apr
	Data = ENU - X_apr;
	Datacov = covENU + Cv*(cur_epoch-ref_epoch)^2;

%% Now form data & covariance, and kernel as relative positions

    d_t = Data;
    cov_t = Datacov;
    G_t = Gatt;

%%  Heavyside function for discontinuties?
%%
%5	step = Heaviside(cur_epoch, Miyake_epoch);
%       whos G_t
%%  construct the submatrix [1,0]

%    sub = zeros(Nbasis, 3*Nbasis);
%    for i = 1:Nbasis; 
%      sub(i,3*(i-1)+1:3*i) = [1, 0, 0];
%    end

	%keyboard
%% Finally form H
%% state vector is statedim=2*Nbasis +3*Nsites(1) +Nframe;
   % keyboard
    H=[Gatt,Ir,Rtr1];
	
	d_t=Data;
	Hd = sparse(H);
	cov_t = sparse(Datacov);





