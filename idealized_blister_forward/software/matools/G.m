	function g = G(t,T,n);

%
% t is time normalized by relaxation time
% T is recurrence interval normalized by relaxation time
% n is order

% initialize terms
	g = zeros(1,length(t));
	eps = 10; k = 0;

while eps > 0.001
	add = exp(-k*T).*(t + k*T).^(n-1);
	g = g + add;
	k = k+1;
	eps = max(abs(add));
end
	g = T.*g.*exp(-t)./fact(n-1);



