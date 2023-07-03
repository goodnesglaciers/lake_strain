function [KM]=keepstations(numstats,keepstats)
%keepstations    [KM]=keepstations(numstats,keepstats)
%
%Makes a matrix that keeps the stations listed in keepstats, and
%discards the rest.

i=length(keepstats);
j=numstats*3;

KM=sparse([],[],[],i*3,j,length(keepstats)*3);
I=eye(3,3);

for t=1:i
	q=(keepstats(t)-1)*3+1;
	r=(t-1)*3+1;
	KM(r:r+2,q:q+2)=I;
end
