function C = mmwnans(A,B)



% first test dimensions
	[N,K1] = size(A);
	[K2,M] = size(B);
		if K1 ~= K2
			disp('dimensions are invalid')
			stop
		end


% dimension output
	C = zeros(N,M);
	
% computation
	
for i = 1:N
	for j = 1:M
		sum = 0;
		for k = 1:K1
			if (A(i,k)==0&isnan(B(k,j))) | (B(k,j)==0&isnan(A(i,k)))
				sum = sum;
			else
				sum = sum + A(i,k)*B(k,j);
			end
		end
		C(i,j) = sum;
	end
end


