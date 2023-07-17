function slip_mplot(pm,rate,i,j, t, times)
%
%
% plots slip on fault patches on 2-D grid
%
% i =  num segments along fault length
% j =  num segments along fault width
%

% set fault parameters
wid = pm(1,2);
len = pm(1,1);
ai = sqrt(-1);
dipr = pm(1,4)*pi/180;
east = pm(:,6);
north = pm(:,7);
up = -1*pm(:,3);

% corners in local coordinates
ll = -0.5*len; lr =0.5*len;
ul = -0.5*len + ai*wid*cos(dipr); ur = 0.5*len + ai*wid*cos(dipr);
%
% z-component
bottom = up;
top = up + wid*sin(dipr);
%
% transform corners
strkr = (90-pm(1,5))*pi/180;

ll = (east+north*ai) + ll*exp(ai*strkr);
lr = (east+north*ai) + lr*exp(ai*strkr);
ul = (east+north*ai) + ul*exp(ai*strkr);
ur = (east+north*ai) + ur*exp(ai*strkr);

xvert = [real(lr)'; real(ur)'; real(ul)'; real(ll)'];
yvert = [imag(lr)'; imag(ur)'; imag(ul)'; imag(ll)'];
zvert = [bottom'; top';top';bottom'];

%len = pm(1,1);
%wid = pm(1,2);
%east = pm(:,6);
%north = pm(:,7);
%mp = [east north];
%
%strike = (pi/180)*pm(1,5);
%
%Rs = [cos(strike), -sin(strike); sin(strike), cos(strike)];
%mp = mp*Rs';
%
%
%
%x1 = mp(:,2)-(0.5*len);
%x2 = mp(:,2)+(0.5*len); 
%
%offset =  min(  min(x1), min(x2) );
%x1 = x1 - offset;
%x2 = x2 - offset;
%
%
%x = [x1 x2];
%
%[np,foo] = size(mp);
%
%y1 = wid*(0:j-1)'*ones(1,np/j); 
%y2 = wid*(1:j)'*ones(1,np/j); 
%
%y = -1*[y1(:) y2(:)];
%
%xvert = [x(:,1)'; x(:,2)'; x(:,2)'; x(:,1)'];
%yvert = [y(:,1)'; y(:,1)'; y(:,2)'; y(:,2)'];
%
% 


%% Determine limits for color scale
	mx = max(rate(:,1));
	mn = min(rate(:,1));
for k = 2:length(t)
	max_now = max(rate(:,k));
	min_now = min(rate(:,k));
	mx = max(mx, max_now);	
	mn = min(mn, min_now);	
end
color_scale = [mn mx];

figure

for i=1:9
%for i = 1:12
        subplot(3,3,i)
          k = times(i);
%         k = 2*i + 2;
%%        subplot(3,4,i)
%%        k = i + 5;
%patch(xvert,yvert,zvert,rate(:,k)'), axis('equal')
patch(xvert,yvert, zvert, rate(:,k)','EdgeColor',[1 1 1]), axis('equal'); view([5, 40])
%patch(xvert,yvert, zvert, rate(:,k)'), axis('equal'); view([5, 40])
axis([-25   170   -15    10   -50     0])
%xlabel('Distance W to E (km)')
%zlabel('Depth (km)')
        caxis(color_scale);
        day = (t(k)-99 )*365.25;
        title(['Day ',num2str( round(day) )])
end

pos=get(gca,'pos');           % store current axis position
h=colorbar;
set(h,'pos',[.94 .2 .02 .6])  % resize colorbar in normalized figure units
set(gca,'pos',pos)            % restore current axis to original size



