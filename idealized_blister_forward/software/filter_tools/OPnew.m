function [RData, E, N, U, RHat, Ep, Np, Up, kk] = OPnew(t, Data, Dhat, X0_enu, sites, rel_site)

%% [RData, E, N, U, RHat, Ep, Np, Up, kk] = OPnew(t, Data, Dhat, X0_enu, sites, rel_site)
%% 
%% Function to pull out Observed "O" and predicted "P" data for plotting.
%% Used with plotOP

%% INPUT:
%%	t 	= vector of observation times
%%	Data	= matrix of observed data
%%	Dhat	= matrix of predicted data
%%	XO_enu	= Apriori coordinates in ENU coordinates
%%	sites   = vector of site names
%%	rel_site= fixed station to which all baselines are formed
%% OUTPUT:
%%	RData 	= Matrix of relative observed data
%%	E,N,U	= East observed, N...
%%	RHat	= Matrix of relative predicted data
%%	Ep,Np,Up= East predicted, N...
%%	kk	= index of rel_site

        [Ncomps,Nepochs] = size(Data);
	Nsites = Ncomps/3;

%% Determine the index of the fixed station
	kk = GetIndex(sites, rel_site);

%% Subtract the apriori coords from observed and predicted
	DData = Data - X0_enu*ones(1,Nepochs);
	DDhat = Dhat - X0_enu*ones(1,Nepochs);


%% Data Series for reference site
	
	RefO = DData(3*(kk-1)+1:3*(kk-1)+3, :);
	RefP = DDhat(3*(kk-1)+1:3*(kk-1)+3, :);

	ObsRef = []; PredRef = [];
for j = 1:Nsites
	ObsRef = [ObsRef;RefO];
	PredRef = [PredRef;RefP];
end

%% Subtract Reference Site motions for Data

	RData = DData - ObsRef;
	RHat  = DDhat - PredRef;


	
% pull out the east north and up components for plotting

for i = 1:Nsites
	E(i ,:) = RData(3*(i-1)+1,:);
	N(i ,:) = RData(3*(i-1)+2,:);
	U(i ,:) = RData(3*(i-1)+3,:);

	Ep(i ,:) = RHat(3*(i-1)+1,:);
	Np(i ,:) = RHat(3*(i-1)+2,:);
	Up(i ,:) = RHat(3*(i-1)+3,:);
end



