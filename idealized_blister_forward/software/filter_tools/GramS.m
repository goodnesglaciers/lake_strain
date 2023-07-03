function  [B_new] = GramS(B);
 
%  [B_new] = GramS(B);
%
%  From new basis functions by taking uniform slip as first basis
%  and using Gram-Schmidt to orthogonalize the other bases relative
%  to it.


	[Nels, Nbasis] = size(B);


%  Form uniform slip basis;
	v = zeros(Nels, Nbasis+1);
	v(:,1)  = ones(Nels, 1)/sqrt(Nels);

%  Do the Gram-Schmidt orthogonalization

	for i = 1:Nbasis

		s = zeros(Nels,1);
		for j = 1:i
		s = s + (  (v(:,j)'*B(:,i) )/(v(:,j)'*v(:,j))  )*v(:,j); 
		end

		v(:,i+1) = B(:,i) - s;

	end

%  Now normalize the v

	B_new(:,1) = v(:,1);
	for i = 2:Nbasis+1
	
		B_new(:,i) = v(:,i)/sqrt(v(:,i)'*v(:,i));

	end
