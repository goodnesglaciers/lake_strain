function [sites, X0] = read_apcoords(obsfil)
%
% function [sites, X0] = read_apcoords(obsfil)
%
% Subroutine to read in station codes from SINEX files
%
% Paul Segall
% July, 1997
%%
        obs = fopen(obsfil,'r');

% warn user if file was not successfully opened

        if obs  == -1
                disp('Unable to open observation file!')
                 return
        end

% rewind file to begining
        frewind(obs)

%% read  station ids
	sites = fscanf(obs, '%s %*g %*g %*g', [4, inf]);
	sites = sites';


% rewind file and read initial coords
        frewind(obs)
	X0 = fscanf(obs, '%*s %g %g %g', [3, inf]);

	fclose(obs);
