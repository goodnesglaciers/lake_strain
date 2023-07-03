function [] = plotpatchslip2(pm,slip,nve)
%
%	 plotpatchslip2(pm,slip,nve)
%
% plots slip on fault patches on 2-D grid
%
% Input:
%	pm = patchmodel (dislocation parameters)
%	pm(:,1) = length;
%	pm(:,2) = width;
%	pm(:,3) = depth;
%	pm(:,4) = dip;
%	pm(:,5) = strike;
%	pm(:,6) = East Offset;
%	pm(:,7) = North Offset;
%	slip = slip values
%	nve  = number of vertical elements
%
% Modified from Susan Owens script by P. Segall to plot 
% distance along fault regardless of strike.


y = -[pm(:, 3) - pm(:,2).*sin(pm(:,4)*pi/180), pm(:,3)];

x = zeros(size(y));
nhe = size(pm,1)/nve;

% for first nve elements
	x(1:nve,2) = pm(1,1);
for k=2:nhe
	index1 = nve*(k-1);
	index2 = nve*(k-2);
	x( index1+1:index1+nve , :) = ...
	[x(index2+1:index2+nve,2), ...
		 x(index2+1:index2+nve,2) + pm(index1+1,1)*ones(nve,1)];
end

xvert = [x(:,1)'; x(:,2)'; x(:,2)'; x(:,1)'];
yvert = [y(:,1)'; y(:,1)'; y(:,2)'; y(:,2)'];

% 
figure
patch(xvert,yvert,slip), axis('equal')
colormap(jet)
colorbar
