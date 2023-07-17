function B=cell2blkdiag(A)
%CELLBLX2DIAG    B=cell2blkdiag(A)
%
%Coverts a cell array of matrices to a single, sparse, block
%diagonal matrix.

M=0;
N=0;
for k=1:length(A)
    [m,n]=size(A{k});
    [i{k},j{k},s{k}]=find(A{k});
    i{k}=i{k}+M;
    j{k}=j{k}+N;
    M=M+m;
    N=N+n;
end

B=sparse(cat(1,i{:}),cat(1,j{:}),cat(1,s{:}),M,N);