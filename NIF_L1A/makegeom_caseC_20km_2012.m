function [G1, G2, G3, G2b, G3b, Rc, Rb, patchesC, patchesB]=makegeom_caseC_20km(tsb,origin)
path('./software/okada85',path)
%
%% LAURA: CRACK GEOMETRY: 2 CRACKS INTERSECT AT APEX OF BEND, mark 2

%       makes the green's function matricies for a set of faults defined in
%       this subroutine given a
%       time series object and an origin.
%
%       G1  nsites*3 (ENU) by n-subfaults matrix for crack strike-slip
%       G2  nsites*3 (ENU) by n-subfaults matrix for crack thrust slip
%       G3  nsites*3 (ENU) by n-subfaults matrix for crack opening
%       G2b nsites*3 (ENU) by n-subfaults matrix for basal thrust slip
%       G3b nsites*3 (ENU) by n-subfaults matrix for basal opening
%       Rc is the laplacian matrix for spatial smoothing of the crack   
%       Rb is the laplacian matrix for spatial smoothing of the base

%hardwired  ASSUMES SAME ORDER AS TS OBJECT 
    load apcoords_lle_2012;
    lats=apcoords_lle_2012(1,:);
    lons=apcoords_lle_2012(2,:);
    hs=apcoords_lle_2012(3,:); % UNITS don't matter, not used
    llh=[lats; lons; hs];
    xy_sta=llh2localxy(llh,origin);        
    Nsites=length(lats);
%%%%%%%%% GOAL IS 3 component Green's functions for each station

%%%%%%%%%%%%%%% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old Segall group mex file setup
% NOTE THE TRICK IS THESE delE and delN are the midpoint of  the lower edge
% I think.  NOT the centroid like the toolbox
%1st get dis_geom from reflection profiles then grid up the fault
             %%  disgeom(1) = len = fault length in strike direction (km)
             %%  disgeom(2) = wid = fault width in dip direction (km)
             %%  disgeom(3) = dep = depth  of lower edge of fault (km)
             %%  disgeom(4) = dip = dip angle (degrees)
             %%  disgeom(5) = strik =  strike, clockwise from N (degrees)
             %%  disgeom(6) = delE  = East offset of midpoint from origin (km)
             %%  disgeom(7) = delN  = North offset of midpoint from origin (km)
             %%  disgeom(8) = ss  =  strike slip motion (m)
             %%  disgeom(9) = ds  =  dip slip motion (m)
             %%  disgeom(10) = op  =  opening  motion (m)
% midpoint,  2 patches in depth
%        clear pm
%        xy_fault = llh2localxy([33.1695, -115.6284]', origin);
%        pm(1,:)=[20,5.2575,10, 108, 36, xy_fault(1), xy_fault(2)];
%        pm(2,:)=[20,5.2575,5, 108, 36, xy_fault(1), xy_fault(2)];	
%        nu=.25;
%        [G1, G2, G3] = MakeKernel(pm, [1,1,0], xy_sta, nu);
% G1=SS, G2=DS, G3=OPENING

% Large Scale Fault Parameters For the Vertical Crack
centroid=[68.723,-49.53,0];
xy_fault = llh2localxy(centroid', origin);
load scarp_ll
crack_h=zeros(464,1);
scarp_ll_h=horzcat(scarp_ll,crack_h);
xy_scarp_llh=llh2localxy(scarp_ll_h', origin);

W=1.0; %km (depth along dip)
dip=90; %vertical (degrees)
buried=0.1; %buried fault (km)

figure
plot(xy_sta(:,1),xy_sta(:,2),'*'); 
hold on; xlabel(' x(km)'); ylabel('y (km)');
plot(xy_fault(1),xy_fault(2),'r*');
title('Local Coordinate System')

% LAURA EDITS CONTINUOUS VERT CRACK
Nvpatches=16;
spacing = 464/Nvpatches;
% location in local xy of plane nodes
scarp_nodes = vertcat(xy_scarp_llh(1:spacing:end,:),xy_scarp_llh(end,:));
scarp_nodes_l=length(scarp_nodes)-1;
for i=1:scarp_nodes_l
    % length of individual segments
    length_segments(i) = sqrt(((scarp_nodes(i,1)-scarp_nodes(i+1,1))^2)+...
        (((scarp_nodes(i,2)-scarp_nodes(i+1,2))^2)));
end
    % strike of individual segments
    delxc = diff(scarp_nodes(:,1)); delyc = diff(scarp_nodes(:,2));
    % piecewise fcn
    del = (90          ) .* (delyc==0) +...
    (atand(delxc./delyc)) .* (delyc>0) + ...
    (90+(atand((abs(delyc))./(delxc)))) .* (delyc<0)

    crack_length=sum(length_segments);
    
% EXTEND AND INTERP       
% start at x = -2
theta = del(1); 
X1=-2; X2= scarp_nodes(1,1); Y2 = scarp_nodes(1,2);
Hstart = abs(((X1-X2))/(sind(theta))); Ydiff = abs(((X1-X2))/(tand(theta)));
Y1=Y2-Ydiff;
% end at x = 2.5
theta2 = 180 - (del(end));
X1_end = 2.5; X2_end = scarp_nodes(end,1); Y2_end = scarp_nodes(end,2);
Hend = abs(((X1_end-X2_end))/(sind(theta)));
Ydiff_end =  abs(((X1_end-X2_end))/(tand(theta)));
Y1_end = Y2_end - Ydiff_end;
% concatenate
XY_start = horzcat(X1,Y1);
XY_end = horzcat(X1_end, Y1_end);
xy_scarp_llh = vertcat(XY_start, xy_scarp_llh, XY_end);
% interp
xq = linspace(-2,2.5,25); %% <--- set how many vert planes wanted HERE!!! n=planes desired+1
vq = interp1(xy_scarp_llh(:,1), xy_scarp_llh(:,2), xq);
xy_scarp_interp = horzcat(xq',vq');

% new del and segment lengths for continuos extended crack
scarp_nodes_l=length(xy_scarp_interp)-1;
for i=1:scarp_nodes_l
    % length of individual segments
    length_segments(i) = sqrt(((xy_scarp_interp(i,1)-xy_scarp_interp(i+1,1))^2)+...
        (((xy_scarp_interp(i,2)-xy_scarp_interp(i+1,2))^2)));
end
    % strike of individual segments
    delxc = diff(xy_scarp_interp(:,1)); delyc = diff(xy_scarp_interp(:,2));
    % piecewise fcn
    del = (90          ) .* (delyc==0) +...
    (atand(delxc./delyc)) .* (delyc>0) + ...
    (90+(atand((abs(delyc))./(delxc)))) .* (delyc<0)

    crack_length=sum(length_segments);

clear disgeom
Nvpatches = 24;

for i=1:1:Nvpatches
    disgeom(i,1) = length_segments(i); % fault length in strike direction (km)
    disgeom(i,2) = W; % fault width in dip direction (km)
    disgeom(i,3) = W+buried; % buried dislocation
    disgeom(i,4) = dip; % dip angle 
    disgeom(i,5) = del(i); % strike, clockwise from N (degrees)
    disgeom(i,6) = (xy_scarp_interp(i,1)+xy_scarp_interp(i+1,1))/2; % East offset of midpoint from origin (km)
    disgeom(i,7) = (xy_scarp_interp(i,2)+xy_scarp_interp(i+1,2))/2; % North offset of midpoint from origin (km)
end

Nvpatches = 24;
Nx = 1; Nz = 6;
Npatches=Nvpatches*Nz;
Nz_tot=Nz;
Nx_tot=Nvpatches;

patches(1:6,:) = patchfault(disgeom(1,1:7),Nx,Nz); %first # is divisions along strike, second is downdip    
patches(7:12,:) = patchfault(disgeom(2,1:7),Nx,Nz);
patches(13:18,:) = patchfault(disgeom(3,1:7),Nx,Nz);
patches(19:24,:) = patchfault(disgeom(4,1:7),Nx,Nz);
patches(25:30,:) = patchfault(disgeom(5,1:7),Nx,Nz);
patches(31:36,:) = patchfault(disgeom(6,1:7),Nx,Nz);
patches(37:42,:) = patchfault(disgeom(7,1:7),Nx,Nz);
patches(43:48,:) = patchfault(disgeom(8,1:7),Nx,Nz);
patches(49:54,:) = patchfault(disgeom(9,1:7),Nx,Nz);
patches(55:60,:) = patchfault(disgeom(10,1:7),Nx,Nz);
patches(61:66,:) = patchfault(disgeom(11,1:7),Nx,Nz);
patches(67:72,:) = patchfault(disgeom(12,1:7),Nx,Nz);
patches(73:78,:) = patchfault(disgeom(13,1:7),Nx,Nz);
patches(79:84,:) = patchfault(disgeom(14,1:7),Nx,Nz);
patches(85:90,:) = patchfault(disgeom(15,1:7),Nx,Nz);
patches(91:96,:) = patchfault(disgeom(16,1:7),Nx,Nz);
patches(97:102,:) = patchfault(disgeom(17,1:7),Nx,Nz);
patches(103:108,:) = patchfault(disgeom(18,1:7),Nx,Nz);
patches(109:114,:) = patchfault(disgeom(19,1:7),Nx,Nz);
patches(115:120,:) = patchfault(disgeom(20,1:7),Nx,Nz);
patches(121:126,:) = patchfault(disgeom(21,1:7),Nx,Nz);
patches(127:132,:) = patchfault(disgeom(22,1:7),Nx,Nz);
patches(133:138,:) = patchfault(disgeom(23,1:7),Nx,Nz);
patches(139:144,:) = patchfault(disgeom(24,1:7),Nx,Nz);


%Material Properties
% Vp=1.8; Vs=0.5;  % in km/s
% nu=(0.5*Vp^2 -Vs^2)/(Vp^2-Vs^2);
nu = 0.3; %PoG p.88 polycrystalline ice

i1=1; i2=i1+Nz;
delx=sqrt((patches(i1,6)-patches(i2,6))^2 + (patches(i1,7)-patches(i2,7))^2); % down strike distance
% delx=abs(patches(1,6)-patches(Nz_tot+1,6));
dely=abs(patches(1,3)-patches(2,3));    % downdip distance change for each patch (crack depth/Nz)
surf=-1; 
[Rc, Rcinv]=modelwt(Nz_tot,Nx_tot,delx,dely,surf);
patchesC=patches;


%%		nve 	= number of vertical elements
%%		nhe 	= number of horizontal elements
%%		delx 	= length of elements in along strike dimension
%%		dely 	= length of elements in dip dimension
%%		surf 	= 1 if fault breaks free surface ~= 1 otherwise
%%  Output:
%%		Lap 	= finite difference Laplacian in two dimensions
%%		Lap_inv = inverse of Lap

figure
plot(xy_sta(:,1),xy_sta(:,2),'*');
hold on
plot(patchesC(:,6),patchesC(:,7),'+')
axis equal
   
%% compute tilts the beaudu toolbox
open=0;
depthc=W/2 + buried;
rake=0;
% FOR RAKE 0=LeftLateral; 180=RL; 90=thrust; -90=normal; A+R convention.
%computes Okada's 1985 solution for displacements, tilts and strains in a
%geographic referential (East,North,Up). Coordinates (0,0) and depth 
%correspond to the fault centroid.
%          E, N   : Coordinates of observation points (relative to fault centroid)
%          DEPTH  : Depth of the fault centroid (DEP > 0)
%          STRIKE : Strike-angle from North (in degrees)
%          DIP    : Dip-angle (in degrees)
%          LENGTH : Fault length in strike direction (LEN > 0)
%          WIDTH  : Fault width in dip direction (WIDTH > 0)
%          RAKE   : Slip-angle direction on fault plane (in degrees)
%          SLIP   : Dislocation in rake direction
%          OPEN   : Dislocation in tensile component
%          NU     : Poisson's ratio
%       returns the following variables (same size as E and N):
%       uE,uN,uZ        : Displacements (Unit of SLIP and OPEN)
%       uZE,uZN         : Tilts (in radian)
%       uNN,uNE,uEN,uEE : Strains (Unit of SLIP and OPEN)/(Unit of N,E,..,WIDTH)

%% BE SUPER CAREfull here;  converting from segall convention where
%% fault patch position is specified as midpoint of lower edge to
%% okada toolbox where position is the centroid of the subfault.   
%% I have hardwired below for a E-W oriented crack and a N-S striking basal fault.
%% have to add in trig functions for other geometries.
%% probably best to stop using patchfault and rewrite that for centroids  

for i=1:Npatches
 % Convert everything to coordinates of fault centroid is at 0,0
 E=xy_sta(:,1)-patches(i,6);
 N=xy_sta(:,2)-patches(i,7);
 L=patches(i,1);
 W=patches(i,2);
 depthc=patches(i,3)-W/2;  % convert from lower edge (segall) to centroid 
 strike=patches(i,5);
 dip=patches(i,4);
 
 %strike-slip
  rake=0; slip=1.0; open=0;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G1(1:Nsites*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           
 
  %thrust
  rake=90; slip=1.0; open=0;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G2(1:Nsites*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           

  %opening
  rake=0; slip=0.0; open=1;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G3(1:Nsites*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           

 %keyboard                  
end
     
plotresults=0;
if(plotresults)
figure
for j=1:8
subplot(4,2,j)
scale=50;

plot(xy_sta(:,1),xy_sta(:,2),'m^','MarkerSize',5); hold on;
  %  plot([patches(j,6)-patches(j,1),patches(j,6)+patches(j,1)],[patches(j,7),patches(j,7)],'g-','LineWidth',3);
     plot([xy_fault(1)+patches(j,6)+((patches(j,1)/Nx)*(sind(patches(j,5)))),...
        xy_fault(1)+patches(j,6)-((patches(j,1)/Nx)*(sind(patches(j,5))))],...
        [xy_fault(1)+patches(j,7)+((patches(j,1)/Nz)*(cosd(patches(j,5)))),...
        xy_fault(1)+patches(j,7)-((patches(j,1)/Nz)*(cosd(patches(j,5))))],...
        'g-','LineWidth',3);
hold on; xlabel(' x(km)'); ylabel('y (km)');
 for i=1:Nsites
    ind1=i; ind2=Nsites+i; ind3=Nsites*2+i;
    quiver(xy_sta(i,1),xy_sta(i,2),scale*G3(ind1,j),scale*G3(ind2,j),'r','filled','LineWidth',2)
    quiver(xy_sta(i,1),xy_sta(i,2),0,scale*G3(ind3,j),'k','filled','LineWidth',2); hold on;
 end
    text(4,4,[num2str(j)])
    axis('equal'); axis([-5 5 -5 5])

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END OF VERTICAL CRACK
%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% START BASE OF ICE SHEET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Large Scale Fault Parameters For the Vertical Crack
centroid=[68.723,-49.53,0];
xy_fault = llh2localxy(centroid', origin);
L=20; %km LAS Nov 2020 -- this used to be 10 km
W=20; %km LAS Nov 2020 -- this used to be 10 km
strike=186;
dip=0.01;  % slight dip down to the west so that is can slip west

% just one patch
buried=.2;
clear disgeom
disgeom(1)=L;
disgeom(2)=W;
disgeom(3)=1.0; % in km?
disgeom(4)=dip;
disgeom(5)=strike;
disgeom(6)=xy_fault(1)*1 -W/2; % east offset of midpoint from origin (km)
disgeom(7)=xy_fault(2)*1;  % north offset of midpoint from origin (km)
                       
%split into subfaults
Nx=12*2; %8 by 8 works good. maybe more is neccesary (Jeff comment)
Ny=12*2; % LAS Nov 2020 -- this used to be 12 for the 10x10km bed fault
patches=patchfault(disgeom,Nx,Ny);  %first # is divisions along strike% second is downdip
Npatchesb=Nx*Ny;
whos patches
%keyboard
%Material Properties
% Vp=1.8; Vs=0.5;  % in km/s
% nu=(0.5*Vp^2 -Vs^2)/(Vp^2-Vs^2);
nu=0.3; %PoG p.88 polycrystalline ice

% Laplacian
%keyboard
i1=1; i2=i1+Nx;
dely=sqrt((patches(i1,6)-patches(i2,6))^2 + (patches(i1,7)-patches(i2,7))^2);
i1=1; i2=2;
delx=sqrt((patches(i1,6)-patches(i2,6))^2 + (patches(i1,7)-patches(i2,7))^2);
surf=-1;
[Rb, Rbinv]=modelwt(Ny,Nx,delx,dely,surf);  % not sure about nx vs ny order?
patchesB=patches;

%keyboard
%% compute tilts the beaudu toolbox
open=0;
depthc=1.0;
rake=0;
% FOR RAKE 0=LeftLateral; 180=RL; 90=thrust; -90=normal; A+R convention.
%computes Okada's 1985 solution for displacements, tilts and strains in a
%geographic referential (East,North,Up). Coordinates (0,0) and depth 
%correspond to the fault centroid.
%          E, N   : Coordinates of observation points (relative to fault centroid)
%          DEPTH  : Depth of the fault centroid (DEP > 0)
%          STRIKE : Strike-angle from North (in degrees)
%          DIP    : Dip-angle (in degrees)
%          LENGTH : Fault length in strike direction (LEN > 0)
%          WIDTH  : Fault width in dip direction (WIDTH > 0)
%          RAKE   : Slip-angle direction on fault plane (in degrees)
%          SLIP   : Dislocation in rake direction
%          OPEN   : Dislocation in tensile component
%          NU     : Poisson's ratio
%       returns the following variables (same size as E and N):
%       uE,uN,uZ        : Displacements (Unit of SLIP and OPEN)
%       uZE,uZN         : Tilts (in radian)
%       uNN,uNE,uEN,uEE : Strains (Unit of SLIP and OPEN)/(Unit of N,E,..,WIDTH)



for i=1:Npatchesb
 % Convert everything to coordinates of fault centroid is at 0,0
 E=xy_sta(:,1)-(patches(i,6)+(patches(i,2)/2)*cosd(strike-180));
 N=xy_sta(:,2)-(patches(i,7) + (patches(i,2)/2)*sind(strike));
 L=patches(i,1);
 W=patches(i,2);
 depthc=patches(i,3)-(W/2)*tan(dip*pi/180);  % convert from lower edge (segall) to centroid 
 depthc=1.0;  % okay to be same for all assume dip is negligble? 

 %figure; plot(E,N,'*'); title(i);
 %keyboard
 %strike-slip
 %rake=0; slip=1.0; open=0;  % strike-slip  in meters
 %[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
 %               depthc,strike,dip,L,W,rake,slip,open,nu);
 % G1(1:Nsites*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           
 
 %thrust
  rake=90; slip=1.0; open=0;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G2b(1:Nsites*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           

  %opening
  rake=0; slip=0.0; open=1;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G3b(1:Nsites*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           
                  
end

%whos G*
%keyboard

plotresults=0;
if(plotresults)
 figure
 for j=1:Npatchesb
  subplot(Nx,Ny,j)
  scale=5;
  plot(xy_sta(:,1),xy_sta(:,2),'m^','MarkerSize',5); hold on;
  rectangle('Position',[patches(j,6),patches(j,7)-patches(j,2)/2,patches(j,1),patches(j,2)]);
  %plot([patches(j,6)-patches(j,1),patches(j,6)+patches(j,1)],[patches(j,7),patches(j,7)],'g-','LineWidth',3);
  hold on; xlabel(' x(km)'); ylabel('y (km)');
  for i=1:Nsites
    ind1=i; ind2=Nsites+i; ind3=Nsites*2+i;
    quiver(xy_sta(i,1),xy_sta(i,2),scale*G3b(ind1,j),scale*G3b(ind2,j),'r','filled','LineWidth',2)
    quiver(xy_sta(i,1),xy_sta(i,2),0,scale*G3b(ind3,j),'k','filled','LineWidth',2); hold on;
  end
    title(['Opening for Base Subfault',num2str(j)])
    axis('equal'); axis([-5 5 -5 5])
 end

 figure
 for j=1:Npatchesb
  subplot(Nx,Ny,j)
  scale=5;
  plot(xy_sta(:,1),xy_sta(:,2),'m^','MarkerSize',5); hold on;
  %keyboard
  rectangle('Position',[patches(j,6),patches(j,7)-patches(j,2)/2,patches(j,1),patches(j,2)]);
  %plot([patches(j,6)-patches(j,1),patches(j,6)+patches(j,1)],[patches(j,7),patches(j,7)],'g-','LineWidth',3);
  hold on; xlabel(' x(km)'); ylabel('y (km)');
  for i=1:Nsites
    ind1=i; ind2=Nsites+i; ind3=Nsites*2+i;
    quiver(xy_sta(i,1),xy_sta(i,2),scale*G2b(ind1,j),scale*G2b(ind2,j),'r','filled','LineWidth',2)
    quiver(xy_sta(i,1),xy_sta(i,2),0,scale*G2b(ind3,j),'k','filled','LineWidth',2); hold on;
  end
    title(['Thrust for Base Subfault',num2str(j)])
    axis('equal'); axis([-5 5 -5 5])
 end


end

%keyboard
% combine  the two; we only want opening of crack,  but both opening and
% dip-slip of the base.

