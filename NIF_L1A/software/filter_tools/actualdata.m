function [datt, index] = actualdata(d)

n = length(d);
index = ones(size(d));

count = 1;
for i = 1:n
	if isnan(d(i))
		index(i) = 0;
	else
		datt(count) = d(i);
		count = count + 1;
	end
end


datt = datt';

