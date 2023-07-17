%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=oneloop(A,B)

% Matrix multiplication with nans


% Get dimensions
[m,k1]=size(A);
[k2,n]=size(B);

if k1 ~= k2
   error('Inner matrix dimensions must agree.')
end

% Set output dimensions
C=zeros(m,n);

% Loop over column vectors of B and C
for i=1:n

   % Transpose column vector and repeat downward "m" times to form a 
   % matrix that can be compared and element multiplied with A.
   bcol=B(:,i)';
   Btest=bcol(ones(m,1),:);

   % make sure NaN times 0 equals zero
   % logic matrix is the same size as A matrix

   logic = (A==0 & isnan(Btest)) | (isnan(A) & Btest==0);
   A(logic)=0;
   Btest(logic)=0;

   % perform inner products for this column of C
   C(:,i)=sum(A.*Btest,2);

end



