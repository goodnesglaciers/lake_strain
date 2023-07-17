   function [t, Data] = cda_part(files, sites, Ncomps, Nepochs, origin)
%% function [t, Data] = cda_part(files, sites, Ncomps, Nepochs, origin)
%%
%% create data array for plotting
%% ignores stations for which there are no matches in the site list

%% Initialize a Data arrray with NaN's

Data  = NaN*ones(Ncomps, Nepochs);
[nfiles, junk] = size(files);

for k=1:nfiles

	%% Read in the data file
        disp(['file:   ', files(k,:)])
        [pos, pos_cov, T, sitecodes] = read_sinex(files(k,:));
		if length(T) > 1
			t(k) = T(1);
		else
			t(k) = T;
		end
        [Nstsob, four] = size(sitecodes);


	%% rotate the positions and covariance into ENU
	[ENU, covENU] = xyz2enu(pos(:), pos_cov, origin);

        for i = 1:Nstsob
                kk = GetIndex(sites, sitecodes(i,:),1);
		if kk ~= 0
			jj = 3*(kk-1)+1;
			Data(jj:jj+2 , k) = ENU(3*(i-1)+1:3*(i-1)+3);
		end
        end

end
