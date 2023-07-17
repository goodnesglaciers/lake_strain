function sitecodes = read_sinex_code(obsfil)
%	sitecodes = read_sinex_code(obsfil)
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


%% read down to start of '+SOLUTION/EPOCHS'
	match=0;

	while match ~= 1
	    line = fgets(obs);
	    [st, count, err] = sscanf(line,'%s');
	    if length(st) >=16
		if st(1:16) == '+SOLUTION/EPOCHS'
			match = 1;
		end
	    end
	end

% read one extra header line
	line = fgets(obs);

%% read station codes until line  '-SOLUTION/EPOCHS' 

	match=0;
	sitecodes = [];
	while match ~= 1
	    line = fgets(obs);
	    [st, count, err] = sscanf(line,'%s');
		if st(1:16) == '-SOLUTION/EPOCHS'
			match = 1;
		else
			sitecodes = [sitecodes; st(1:4)];	
		end
	end


%% close file

	fclose(obs);
