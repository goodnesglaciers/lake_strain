	function [slip] = integrate_rate(rate);

%% 	function [slip] = integrate_rate(rate);
%%
%% integrate the slip-rate to determine slip

	Nepochs = size(rate,2);

	slip = zeros(size(rate));

	for k = 1:Nepochs
		temp = rate(:,1:k);
		temp2 = sum(temp')/365.25;
		slip(:,k) = temp2';
	end


