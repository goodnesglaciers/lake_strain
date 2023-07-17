
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=twoloops(A,B)

%Matrix Multiplication With NaN's

% first test dimensions
        [N,K1] = size(A);
        [K2,M] = size(B);
                if K1 ~= K2
                        disp('dimensions are invalid')
                        stop
                end

% dimension output
        C = zeros(N,M);

% Loop over elements of C performing each inner product
for i = 1:N
        for j = 1:M

		k = 1:K1;
		arow=A(i,k);
		bcol=B(k,j);
   		logic = (arow==0 & isnan(bcol') | isnan(arow) & bcol'==0);
		arow(logic)=0;
		bcol(logic)=0;

                C(i,j) = arow*bcol;	% inner product
        end
end

