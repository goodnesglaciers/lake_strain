function	[ts3]=mergegpsts(ts1,ts2);
%
% merges two gps time series,  according to the merge function in the 
% ident toolbox, there is no smart way to do this other than loadding all of the
% two objects data into the arrays and then calling gpsts to create a new object.
% 
% This piece of code is not working. Not sure how to handle siteindex
%

n1=length(ts1.d); n2=length(ts2.d);
d=[ts1.d; ts2.d;];

% The data covariance is why this is tricky.  not enough memory to
% just concatenate the arrays, so we would be basically redoing all the
% loops in readGPSTimeSeries3.m to make a cell array and then do cell2blkdiag
% Maybe there's no faster way to do this because of the memory limits for big arrays?
dcov=[ts1.dcov, zeros(n1,n2); 
      zeros(n2,n1), ts2.dcov];

sites=[ts1.sites, ts2.sites];
epochs=[ts1.epochs;ts2.epochs];
epochindex=[ts1.epochindex;ts2.epochindex];
apcoords=[ts1.apcoords,ts2.apcoords];

masterlist=unique(sites,'rows');

siteindex = zeros(size(d));
I=ones(size(sites,1),1)

for i=size(masterlist,1):-1:1
    site=masterlist(i,:);
    b=(find(~sum(site(I,:)~=sites,2))-1)*3+1;
    b=[b,b+1,b+2]';
    siteindex(b,1)=i;
end

ts3=gpsts(d,dcov,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apcoords);

