	function ind = date2ind(dte, t)

% Given dates return the indicies of the time vector
% Input:
%	t = vector of observations times
%	dte = vector of dates to find indicies for

diff =  ones(length(dte), length(t)); diff = NaN*diff;

	
for i = 1:length(t)
	diff(:,i) = (dte - t(i)*ones(size(dte)))';
end

[y,ind] = min(abs(diff'));

