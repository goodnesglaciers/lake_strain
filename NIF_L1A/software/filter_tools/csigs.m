   function [E_Sig, N_Sig, U_Sig] = csigs(files, sites, Ncomps, Nepochs, ...
								origin, rel_site)
%% function [E_Sig, N_Sig, U_Sig] = csigs(files, sites, Ncomps, Nepochs, ...
%%								origin, rel_site)
%%
%% create array of sigmas (observed standard deviations for each
%% component in ENU)  for the time series relative to "rel_site".
%% If a station is not observed (including the rel_site) leave a NaN
%% This is to be used in plotting observed and predicted data after
%% running the filter


%% INPUT:
%%	files	=  list of data file names
%%	sites	=  list of site names
%%	Ncomps  =  Number of components 3*Number of sites
%%	Nepochs =  Number of observation epochs
%%	origin  =  origin for rotating coordinates to ENU
%%	rel_site=  fixed station to form all baselines relative to

%% Initialize a Data arrray with NaN's

Rel_Sigs  = NaN*ones(Ncomps, Nepochs);
[nfiles, junk] = size(files);

for k=1:nfiles

	%% Read in the data file
        [pos, pos_cov, T, sitecodes] = read_sinex(files(k,:));
                if length(T) > 1
                        t(k) = T(1);
                else
                        t(k) = T;
                end
        [Nstsob, four] = size(sitecodes);

	%% rotate the positions and covariance into ENU
	[ENU, covENU] = xyz2enu(pos(:), pos_cov, origin);

	%% get index of rel_site
	ii = GetIndex(sitecodes, rel_site, 1);
	if ii ~= 0 

	%% variance of fixed site
	fxd_st_sgs = diag(  covENU(3*(ii-1)+1:3*(ii-1)+3,3*(ii-1)+1:3*(ii-1)+3));

      	  for i = 1:Nstsob
                kk = GetIndex(sites, sitecodes(i,:),1);
		jj = 3*(kk-1)+1;

		%% variance of second site
		snd_st_sgs = diag(  covENU(3*(i-1)+1: 	...
			3*(i-1)+3,3*(i-1)+1:3*(i-1)+3));

		%% Get covariance terms
		cov_term =  diag(   covENU(3*(ii-1)+1:3*(ii-1)+3, ...
			3*(i-1)+1:3*(i-1)+3));

		%% covariance of the difference
		
		covdiff = fxd_st_sgs + snd_st_sgs - 2*cov_term;

		%% where to put the result in the array

		Rel_Sigs(jj:jj+2 , k) = sqrt(covdiff);

      	  end

	end

end



E_Sig = zeros(Ncomps/3,Nepochs);
N_Sig = zeros(Ncomps/3,Nepochs);
U_Sig = zeros(Ncomps/3,Nepochs);

%% Extract the sigmas and associate with each component

for k = 1:Ncomps/3
	E_Sig(k,:) = Rel_Sigs(3*(k-1)+1,:);
	N_Sig(k,:) = Rel_Sigs(3*(k-1)+2,:);
	U_Sig(k,:) = Rel_Sigs(3*(k-1)+3,:);
end
