	function [v, gamma] = savpres(x,t,H,tR,T,N);

%	 [v, gamma] = savpres(x,t,H,tR,T);
%
% Postseismic displacements due to slip on fault of
% depth H, in an elastic layer of thickness H, over
% a Maxwell viscoelastic half-space. Savage (1990)
%
%	tR = Maxwell relaxation time in same units as time
%	T = recurrence interval

v = zeros(length(x),length(t));
if nargout == 2
	gamma = zeros(length(x),length(t));	
end

% main loop
if nargin ==6
	for n=1:N
		g = G(t/tR, T/tR, n);
		Fn = atan(2*H*x./( (4*n^2-1)*H^2  + x.^2  ))/pi;
		v = v + Fn'*g;

		if nargout == 2
			d1 = (2*n+1)*H;    d2 = (2*n-1)*H;
			Gn = (-d1./(x.^2 + d1^2)  + d2./(x.^2 + d2^2))/pi;
		        gamma = gamma + Gn'*g;	
		end

	end
	n = n+1; 
end

if nargin ==5
	tol = 1e-3; n = 1; vnew = v +100;

	while norm(vnew,1) > tol*norm(v,1) & n < 50
		g = G(t/tR, T/tR, n);
		Fn = atan(2*H*x./( (4*n^2-1)*H^2  + x.^2  ))/pi;
		vnew = Fn'*g;
		v = v + vnew;
		if nargout == 2
			d1 = (2*n+1)*H;    d2 = (2*n-1)*H;
			Gn = (-d1./(x.^2 + d1^2)  + d2./(x.^2 + d2^2))/pi;
		        gamma = gamma + Gn'*g;	
		end

		n = n+1;
	end
end



%% deep dislocation below
  	Fn =  + atan( x/((2*n-1)*H))/pi;
 	g = ones(size(g));
 	v = v + Fn'*g;


	if nargout == 2
		d2 = (2*n-1)*H;
		Gn = (d2./(x.^2 + d2^2))/pi;
		gamma = gamma + Gn'*g;	
	end





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



