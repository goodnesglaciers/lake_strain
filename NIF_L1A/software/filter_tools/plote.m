  function plote(t,d, sigs)
%% function plote(t,d, sigs)
%% 
%%	Plot data so lines don't overlap
%%	Modified to work with incomplete data
%%
%%	Input:
%%			t = vector of times
%%			d = matrix of obsrvations
%%			sigs = matrix of standard deviations

	[Nsites,Nepochs] = size(d);
	for i = 1:Nsites
		j =1;
		while isnan(d(i,j)) == 1;
			j = j+1;
		end
		d(i,:) = d(i,:) - d(i,j)*ones(1,Nepochs) ;
	end


	dp = zeros(size(d));
	dp(1,:) = d(1,:);
	errorbar(t, dp(1,:), sigs(1,:)); hold on
	
	for k=2:Nsites
	dp(k,:) = d(k,:) + max( stripnan(dp(k-1,:))) - min( stripnan(d(k,:)));
	errorbar(t, dp(k,:), sigs(k,:)); hold on

	end
	

	xlabel('Time (years)')
	ylabel('Displalecement (mm)')
