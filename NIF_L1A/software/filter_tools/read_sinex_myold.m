function [x, cov, obs_date, sitecodes] = read_sinex(data_file)
%	 [x, cov, obs_date, sitecodes] = read_sinex(data_file)
%
% Subroutine to read data from SINEX files
% using Yosuke Aoki's snx2mfiln
% Paul Segall
% July, 1997



	[Obs_date,X,COV,ndata]=snx2mfiln(data_file);

	[Obs_date,x,cov]=FormMatrix(Obs_date,X,COV,ndata);

	obs_date = Obs_date(1,1);

	sitecodes = read_sinex_code(data_file);

	[nsites, four] = size(sitecodes);

	if nsites ~= ndata
		disp('warning: site codes and data vector not same length')
		return
	end

