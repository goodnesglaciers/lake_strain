function [sites, sec_vel, sec_sig] = read_apvel(obsfil)
%
% function [sites, sec_vel, sec_sig] = read_apvel(obsfil)
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
	sites = fscanf(obs, '%s %*g %*g %*g %*g %*g %*g', [4, inf]);
	sites = sites';


% rewind file and read a priori velocities
% Format is N; N_sig;  E; E_sig; V; V_sig in m/yr

        frewind(obs)
	r = fscanf(obs, '%*s %g %g %g %g %g %g', [6, inf]);
	r = r';
	sec_vel = [r(:,1), r(:,3), r(:,5)];
	sec_sig = [r(:,2), r(:,4), r(:,6)];

	fclose(obs);
