function [patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(origin,depthc)
path('../software/okada85',path)
%% LAURA: FAULT GEOMOETRY FOR JUST THE BED
% 2021 Feb 25: build Green function matrices for a surface grid
%       makes the Green function matrices for a set of faults defined in
%       this subroutine, without being given a time series object.
%  At locations of the entire 20x20 km surface 
%       G1v  nsurface*3 (ENU) by n-surface-locations matrix for crack strike-slip 
%       G2v  nsurface*3 (ENU) by n-surface-locations matrix for crack thrust slip
%       G3v  nsurface*3 (ENU) by n-surface-locations matrix for crack opening
%       G2bv nsurface*3 (ENU) by n-surface-locations matrix for basal thrust slip
%       G3bv nsurface*3 (ENU) by n-surface-locations matrix for basal opening
%       Rc is the laplacian matrix for spatial smoothing of the crack   
%       Rb is the laplacian matrix for spatial smoothing of the base

% GNSS stations
    load('apcoords_lle.mat'); lats=apcoords_lle(1,:); lons=apcoords_lle(2,:);
    llh=[lats; lons];
    xy_sta=llh2localxy(llh,origin);        
    Nsites=length(lats);

%%%%%%%%% GOAL IS 3-component Greens functions for every surface location

% xy for the 20km x 20 km surface
[xx,yy] = meshgrid(-15:0.25:15, -15:0.25:15); %250-m grid
xy_surf = horzcat(reshape(xx,121*121,1), reshape(yy,121*121,1));
Nsurface=length(xy_surf);

%%%%%%%%%%%%%%% GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% old Segall group mex file setup
% NOTE THE TRICK IS THESE delE and delN are the midpoint of  the lower edge
% I think.  NOT the centroid like the toolbox
% 1st get dis_geom from reflection profiles then grid up the fault
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

%%%%%% START BASE OF ICE SHEET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Large Scale Fault Parameters For the Horizontal Crack
centroid=[68.723,-49.53,0];
xy_fault = llh2localxy(centroid', origin);
L=20; % km 
W=20; % km 
strike=180;
dip=0.01;  % slightly dip down to the west

% just one patch
clear disgeom
disgeom(1)=L;
disgeom(2)=W;
disgeom(3)= depthc; % DEPTH [ km ]
disgeom(4)=dip;
disgeom(5)=strike;
disgeom(6)=xy_fault(1)*1 - (W/2) + 0.5; % east offset of midpoint from origin (km)
disgeom(7)=xy_fault(2)*1 - (0.0846474982469 + 0.25);  % north offset of midpoint from origin (km)
                       
%split into subfaults
Nx=40; % 8 by 8 works good. maybe more is neccesary (Jeff comment)
Ny=40; % LAS Nov 2020 -- this used to be 12 for the 10x10km bed fault
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
% i1=1; i2=i1+Nx;
% dely=sqrt((patches(i1,6)-patches(i2,6))^2 + (patches(i1,7)-patches(i2,7))^2);
% i1=1; i2=2;
% delx=sqrt((patches(i1,6)-patches(i2,6))^2 + (patches(i1,7)-patches(i2,7))^2);
% surf=-1;
% [Rb, Rbinv]=modelwt(Ny,Nx,delx,dely,surf);  % not sure about nx vs ny order?
patchesB=patches;

%keyboard
%% compute tilts the beaudu toolbox
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


% for entire surface
for i=1:Npatchesb
 % Convert everything to coordinates of fault centroid is at 0,0
 E=xy_surf(:,1)-(patches(i,6)+(patches(i,2)/2)*cosd(strike-180));
 N=xy_surf(:,2)-(patches(i,7) + (patches(i,2)/2)*sind(strike));
 L=patches(i,1);
 W=patches(i,2);
 
 %thrust
  rake=90; slip=1.0; open=0;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G2bv(1:Nsurface*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           
  G2bv_strain(1:Nsurface*6,i)=[uZE; uZN; uNN; uNE; uEN; uEE]';   % strains  
 
  %opening
  rake=0; slip=0.0; open=1;  % strike-slip  in meters
  [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = okada85(E,N,...
                depthc,strike,dip,L,W,rake,slip,open,nu);
  G3bv(1:Nsurface*3,i)=[uE; uN; uZ]';   % block order not  ENU for each station           
  G3bv_strain(1:Nsurface*6,i)=[uZE; uZN; uNN; uNE; uEN; uEE]';   % strains  
             
end

%whos G*
%keyboard

% plotresults=0;
% if(plotresults)
%  figure
%  for j=1:Npatchesb
%   subplot(Nx,Ny,j)
%   scale=5;
%   plot(xy_sta(:,1),xy_sta(:,2),'m^','MarkerSize',5); hold on;
%   rectangle('Position',[patches(j,6),patches(j,7)-patches(j,2)/2,patches(j,1),patches(j,2)]);
%   hold on; xlabel(' x(km)'); ylabel('y (km)');
%   for i=1:Nsites
%     ind1=i; ind2=Nsites+i; ind3=Nsites*2+i;
%     quiver(xy_sta(i,1),xy_sta(i,2),scale*G3b(ind1,j),scale*G3b(ind2,j),'r','filled','LineWidth',2)
%     quiver(xy_sta(i,1),xy_sta(i,2),0,scale*G3b(ind3,j),'k','filled','LineWidth',2); hold on;
%   end
%     title(['Opening for Base Subfault',num2str(j)])
%     axis('equal'); axis([-5 5 -5 5])
%  end

 figure
 for j=1:Npatchesb
  subplot(Nx,Ny,j)
  scale=5;
  plot(xy_sta(:,1),xy_sta(:,2),'m^','MarkerSize',5); hold on;
  %keyboard
  rectangle('Position',[patches(j,6),patches(j,7)-patches(j,2)/2,patches(j,1),patches(j,2)]);
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
