%Number of faults
nf = 1;
%
%Plotting range
xmin = eval( get(getxmin,'string') );
xmax = eval( get(getxmax,'string') );
ymin = eval( get(getymin,'string') );
ymax = eval( get(getymax,'string') );

nu = 0.25;
%
% First fault segment
len = eval( get(getlen,'string') );
wid = eval( get(getwid,'string') );
dep = eval( get(getdep,'string') );
dip = eval( get(getdip,'string') );
strik = eval( get(getstr,'string') );
delE = eval( get(getdle,'string') );
delN = eval( get(getdln,'string') );
ss = eval( get(getsss,'string') );
ds = eval( get(getdss,'string') );
op = eval( get(getops,'string') );
%
dis_geom(1,:) = [len,wid,dep,dip,strik,delE,delN,ss,ds,op];
%
%reciever geometry
y = linspace(ymin, ymax,  20);
x = linspace(xmin, xmax,  20);
%y = linspace(ymin, ymax,  2);
%x = linspace(xmin, xmax,  2);
n = length(x);  m = length(y);
%
%initialize plot arrays to zero
uz = zeros(n,m); ux = zeros(n,m);  uy=zeros(n,m);
%
for i= 1:m,
  for j=1:n,
     for k=1:nf,
     u = disloc(nu, dis_geom(k,:),[x(j),y(i)]);
     uz(i,j)  = uz(i,j) + u(3);
     ux(i,j)  = ux(i,j) + u(1);
     uy(i,j)  = uy(i,j) + u(2);
     end
  end
end
%

%% PLOTTING

figure
C=contour(x,y,uz);
clabel(C);
hold on
quiver(x,y,ux,uy)
hold on
for k=1:nf
  displot(dis_geom(k,:));
end
hold off

figure
mesh(x,y,uz)

figure
quiver(x,y,ux,uy)
hold on
for k=1:nf
  displot(dis_geom(k,:));
end
