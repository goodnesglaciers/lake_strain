	function basis_mplot(pm,Bp,i,j)
%
% 	basis_mplot(pm,Bp,i,j)
%
%
% plots basis functions on multiple 2-D grids
%
% i =  num segments along fault length
% j =  num segments along fault width
%
 
wid = pm(1,2);
len = pm(1,1);
east = pm(:,6);
north = pm(:,7);

mp = [east north];
 
strike = (pi/180)*pm(1,5);
 
Rs = [cos(strike), -sin(strike); sin(strike), cos(strike)];
mp = mp*Rs';
 
x1 = mp(:,2)-(0.5*len);
x2 = mp(:,2)+(0.5*len);

x = [x1 x2];
 
np = size(mp,1);
 
y1 = wid*(0:j-1)'*ones(1,np/j);
y2 = wid*(1:j)'*ones(1,np/j);

y = -1*[y1(:) y2(:)];
 
xvert = [x(:,1)'; x(:,2)'; x(:,2)'; x(:,1)'];
yvert = [y(:,1)'; y(:,1)'; y(:,2)'; y(:,2)'];

%% Determine limits for color scale
        mx = max(Bp(:,1));
        mn = min(Bp(:,1));

for k = 2:size(Bp,2)
        max_now = max(Bp(:,k));
        min_now = min(Bp(:,k));
        mx = max(mx, max_now);
        mn = min(mn, min_now);
end
color_scale = [mn mx];

figure

Nterms = size(Bp,2);
Nplot = floor(Nterms/2);
for i=1:6
        subplot(Nplot,2,i)
        patch(xvert,yvert,Bp(:,i)'), axis('equal')
        caxis(color_scale);
        title(['Basis ',num2str( i )])
end


pos=get(gca,'pos');           % store current axis position
h=colorbar;
set(h,'pos',[.94 .2 .02 .6])  % resize colorbar in normalized figure units
set(gca,'pos',pos)            % restore current axis to original size
