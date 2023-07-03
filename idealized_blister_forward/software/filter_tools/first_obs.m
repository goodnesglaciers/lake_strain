function [d0] = first_obs(d)


%% Input Data Matrix
%% Return a vector with the first observation for each site
%%

[n, nepochs] = size(d);

for i = 1:n
	j = 1;
	while isnan(d(i,j)) == 1
	j = j+1;
	end
	
	d0(i) = d(i,j);
end



