%% Nsurface principle stress paper figure
% 16 April 2021 LAS -- single time point, expanded xx_fine 
% 3 May 2021 LAS -- change plot area to 15km by 15 km, add d_max vector
% 21 June 2023 LAS -- updated with Stacy Larochelle infinitesimal strain
% calcs, reviewer comments.
clear all; close all

%% 2011 saved NIF outputs and forward-modeled stresses
load('../NIF_L1A/Gsurface500_strain_surfacecrack.mat') % Green functions for 20 km X 20 km surface using
% 0.5-km spacing across surface and the NL set-up
load time_2011.mat; % time vector

% xy for the 20km x 20 km surface plane (Nsurface solutions)
[xx_fine,yy_fine] = meshgrid(-20:0.5:20, -20:0.5:20); % [ km ]
xy_surf_fine = horzcat(reshape(xx_fine,81*81,1),reshape(yy_fine,81*81,1)); % [ km ]
xy_surf = Gsurface500.xy_surf; % [ km ]
Nsurface = Gsurface500.Nsurface; % length(xy_surf)

% load saved NIF forward-modeled stresses 
load n2011_allfour_surfC.mat % 2011

%% winter velocities load and interpolate to okada85 output
origin = [68.72, -49.53];
load('stresses_nlake_winter_vel_23Sept2020.mat');
    lat_vec = reshape(stresses_nlake_winter_vel.lat_grid_sub, 60*42, 1);
    lon_vec = reshape(stresses_nlake_winter_vel.lon_grid_sub, 60*42, 1);
    Sxx_flow_vec = reshape(stresses_nlake_winter_vel.Sxx_flow, 60*42, 1);
    Sxx_flow_vec(isnan(Sxx_flow_vec))=0;
    Syy_flow_vec = reshape(stresses_nlake_winter_vel.Syy_flow, 60*42, 1);
    Syy_flow_vec(isnan(Syy_flow_vec))=0;
    Sxy_flow_vec = reshape(stresses_nlake_winter_vel.Sxy_flow, 60*42, 1);
    Sxy_flow_vec(isnan(Sxy_flow_vec))=0;
    llh_stresses = [lat_vec'; lon_vec'; zeros(2520,1)'];
    xy_stresses = llh2localxy(llh_stresses,origin);
    
% principals 
princ_sigma1_winter = (0.5.*(Sxx_flow_vec+Syy_flow_vec))+...
    (sqrt(((0.5.*(Sxx_flow_vec-Syy_flow_vec)).^2) + (Sxy_flow_vec.^2)));
princ_sigma2_winter = (0.5.*(Sxx_flow_vec+Syy_flow_vec))-...
    (sqrt(((0.5.*(Sxx_flow_vec-Syy_flow_vec)).^2) + (Sxy_flow_vec.^2)));
theta_winter = (0.5.*(atan2((2.*Sxy_flow_vec),(Sxx_flow_vec-Syy_flow_vec)))); % rotate 7 degrees to get into 270degrees = -y axis
theta_winter_deg = rad2deg(theta_winter)+7; % rotate 7 degrees to get into 270degrees = -y axis
theta_winter_rad = deg2rad(theta_winter_deg);

% interp to X Y of okada85 output
    F = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),Sxx_flow_vec);
    Sxx_flow_vec_vq = F(xx_fine,yy_fine);
    FF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),Syy_flow_vec);
    Syy_flow_vec_vq = FF(xx_fine,yy_fine);
    FF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),Sxy_flow_vec);
    Sxy_flow_vec_vq = FF(xx_fine,yy_fine);
    FFF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),princ_sigma1_winter);
    princ_sigma1_winter_vq = FFF(xx_fine,yy_fine);
    FFFF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),princ_sigma2_winter);
    princ_sigma2_winter_vq = FFFF(xx_fine,yy_fine);
    FFFFF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),theta_winter_deg);
    theta_winter_vq_deg = FFFFF(xx_fine,yy_fine);
    theta_winter_vq_rad = deg2rad(theta_winter_vq_deg);

%% load L1A geographic files
origin = [68.72, -49.53];
load moulin_2011.mat
load out2011_20kmB_third.mat % vertical and bed patches rotated into along-flow
    patches_rotate = out2011_20kmB_third.patches;
    patchesC = out2011_20kmB_third.patches_C;
load out2011_20kmB_surfC.mat % vertical and bed output from 20 km NIF L1A
    patchesB = out2011_20kmB_surfC.patchesB;
load apcoords_lle
    lats=apcoords_lle(1,:); lons=apcoords_lle(2,:); hs=apcoords_lle(3,:);
    llh=[lats; lons; hs];
    xy_sta_11=llh2localxy(llh,origin);
    llh_moulin = [moulin_2011(:,4)'; moulin_2011(:,3)'; zeros(37,1)'];
    xy_moulin = llh2localxy(llh_moulin,origin);
load lake.mat
    llh_lake = [lake(:,5)'; lake(:,4)'; zeros(337,1)'];
    xy_lake = llh2localxy(llh_lake,origin);
load lil_lake_2013_168.mat
    llh_lil_lake = [lil_lake_2013_168(:,4)'; lil_lake_2013_168(:,3)'; zeros(148,1)'];
    xy_lil_lake = llh2localxy(llh_lil_lake,origin);
load('out2011_1day.mat') % GPS displacements at precursor, max open, 1-day
% polar stereographic needed values 
    radius=6378137.0;  eccen=0.08181919; lat_true=70; lon_posy=-45;
NNL=csvread('NNL_20110617.csv',1,0); 
    [phi,lambda]=polarstereo_inv(NNL(:,1),NNL(:,2),radius,eccen,lat_true,lon_posy);
    llh_nnl = [phi,lambda]; 
    xy_nnl_lake = llh2localxy(llh_nnl',origin);
% L1C crack from 2019/August 10th
    NNL2019=csvread('NNL20190810B03_latlon.csv',1,0);
    xy_NNL_crack = llh2localxy(horzcat(NNL2019(:,2), NNL2019(:,1))',origin);
    strike_NNL = atand((xy_NNL_crack(end,2)-xy_NNL_crack(1,2))./...
        (xy_NNL_crack(end,1)-xy_NNL_crack(1,1)) ) +180 % +180 for plotting purposes
% L1B from 2013
load SNL_crack_2013_186.mat
    xy_SNL_crack = llh2localxy(horzcat(SNL_crack_2013_186(:,2),SNL_crack_2013_186(:,1))',origin);   
    strike_SNL = atand((xy_SNL_crack(end,2)-xy_SNL_crack(1,2))./...
        (xy_SNL_crack(end,1)-xy_SNL_crack(1,1)) ) +180 % +180 for plotting purposes
% strike L1A crack in lake    
    x_crack = patchesC(1:6:end,6); y_crack = patchesC(1:6:end,7);
    strike_NL = atand((y_crack(22,1)-y_crack(8,1))./(x_crack(22,1)-x_crack(8,1))) +180  % +180 for plotting purposes
load terravel_winter.mat % winter TSX velocities

%% figure with principal winter background stresses + lake-deformation stresses
fig1=figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[1 1 188 108]./5);
tt=0.22; ss = 188/110; scale=0.65;

axe1 = axes('Position',[0.035 0.748 0.19*scale 0.19*ss*scale],'Box','on');
axe21 = axes('Position',[0.035 0.645 0.19*scale 0.05],'Box','on','xticklabel',[],'yticklabel',[]);
axe2 = axes('Position',[0.035 0.395 0.19*scale 0.19*ss*scale],'Box','on','xticklabel',[],'yticklabel',[]);
axe4 = axes('Position',[0.035 0.145 0.19*scale 0.19*ss*scale],'Box','on','xticklabel',[],'yticklabel',[]);

axe3 = axes('Position',[0.22 0.588 tt tt*ss],'Box','on','xticklabel',[],'yticklabel',[]); 
axe5 = axes('Position',[0.22 0.148 tt tt*ss],'Box','on','yticklabel',[]);

axe32 = axes('Position',[0.445 0.588 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 
axe52 = axes('Position',[0.445 0.148 tt tt*ss],'Box','on','yticklabel',[]);

axe33 = axes('Position',[0.67 0.588 tt tt*ss],'Box','on','xticklabel',[]); 
axe53 = axes('Position',[0.67 0.148 tt tt*ss],'Box','on');

m=10; % fontsize
TriangleSize=5;
scale_factor = 14; scale_factor2 = 10; % quiver scales
v=[-8005:10:8005]; v200=200; % stress contours [ kPa ]
load BWR.mat % colormap
% colorbar ranges
cmin_vert = -0.15; cmax_vert = 0.15; % vertical crack
cmin_bup = -1.0; cmax_bup = 1.0; % bed opening
cmin_bslip = -0.5; cmax_bslip = 0.5; % bed slip

axes(axe1); % winter vel
text(-14.5, 11, 'a','FontWeight','bold','FontSize',12); hold on;
vc=[60:2:120];
[C1,h1]=contourf(terravel_winter.polar.xx_moulin./1e3, terravel_winter.polar.yy_moulin./1e3,...
    terravel_winter.winter.vel,vc); set(h1,'LineColor','none'); 
hold all; caxis([60 120])
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none');
load('V0_enu_good_2011.mat'); % pre-drainage GPS velocities from ENU displacements
v_east = V0_enu(1:3:end); v_north = V0_enu(2:3:end);
quiver(xy_sta_11(1:15,1), xy_sta_11(1:15,2),...
    v_east',v_north','k','AutoScale','on','MaxHeadSize',0.5,'LineWidth',1)
% velocity scale bar
x_qscalebar = [7.5]; y_qscalebar = [-9.3]; u_qscale = [-0.5]; v_qscale = [0];
quiver(x_qscalebar, y_qscalebar, u_qscale*scale_factor, v_qscale*scale_factor,'k','AutoScale','on','MaxHeadSize',0.5,'LineWidth',1)
text(3,-8,'0.5 m d^{-1}','FontSize',9)
% north arrow
x_north = [-8.8]; y_north = [8]; u_north = [0]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(-9.366666666,7.7,'N','FontSize',11)
yl = ylabel(' N-S [ km ]'); xl = xlabel(' E-W [ km ]');
xl.Position(2) = 12;  % Shift x label down
xl.Position(1) = 0;  % Shift x label down
xlim([-10 10]); ylim([-10 10]);
set(gca,'xtick',[-10:5:10],'ytick',[-10:5:10],'tickdir','in','LineWidth',1.1,'FontSize',m-2,...
    'XAxisLocation','top','Layer','Top'); 
colormap(axe1,parula); 
t1=colorbar('EastOutside'); set(t1,'YTick',[40,60,80,100,120]); hold all
set(get(t1,'ylabel'),'String','Surface Velocity [ m yr^{-1} ]','FontSize',m-1);
set(t1, 'Position', [.035+(0.19*scale)+0.005 0.76 0.003 (0.19*ss*scale)-0.03]);
set(gca,'FontName','Avenir');
hold on

axes(axe21); % VERTICAL CRACK
hold on;
z=out2011_20kmB_surfC.vert_crack_open_m;
 Nx=24;
 Nz=6; 
 isf=0;
 for i=1:Nx
  for j=1:Nz
    isf=isf+1;
    x1=patchesC(isf,6)-patchesC(isf,1)/2;  % less patch width
    x2=x1+patchesC(isf,1);
    y1=patchesC(isf,3);
    y2=y1-patchesC(isf,2); % bottom fault goes shallower
    x=[x1, x2, x2, x1];
    y=[y1, y1, y2, y2];
    patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6]); hold on;
 
  end
 end
 colormap(axe21,BWR); 
 caxis([cmin_vert cmax_vert]); 
 set(gca,'ydir','rev');  
 xlim([-2.01 2.51]); ylim([0.099 1.11]);
 set(gca,'FontSize',m-1)
text(-2.8, 0.2, 'b','FontWeight','bold','FontSize',12);
% 1 km scale bar
x_qscalebar = [-1.0]; y_qscalebar = [0.99]; u_qscalebar = [-1]; v_qscalebar = [0];
quiver(x_qscalebar, y_qscalebar, u_qscalebar, v_qscalebar,'k','AutoScale','on','ShowArrowHead','off','LineWidth',2)
text(-1.78,0.8,'1 km','FontSize',9,'FontWeight','bold')
% north arrow
x_north = [-1.9]; y_north = [0.55]; u_north = [0]; v_north = [-0.4]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(-1.75,0.35,'z','FontSize',9)
axis off
title('Hydro-fracture Opening');
tt21=colorbar('EastOutside'); set(tt21,'YTick',[-0.15, 0, 0.15]); hold all
set(get(tt21,'ylabel'),'String','[ m ]','FontSize',m-1);
set(tt21, 'Position', [.035+(0.19*scale)+0.005 0.645 0.003 .05]);
set(gca,'FontName','Avenir');

axes(axe2); % Basal Thrust at Max H-F Opening
hold on;
z=-1.*out2011_20kmB_surfC.bed_crack_slip_m;
 Nx=24;
 Ny=24; 
 isf=0;
 
 for i=1:Nx
  for j=1:Ny
    isf=isf+1;
    
    x1 = patches_rotate(isf,1);
    x2 = x1 - patchesB(isf,1); % less patch width
    y1 = patches_rotate(isf,2);
    y2 = y1 - patchesB(isf,2); % less patch width
    x=[x1, x2, x2, x1];
    y=[y1, y1, y2, y2];
    patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6]); hold on;
  
  end
 end
colormap(axe2,BWR); 
caxis([-0.5 0.5]);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',2)

quiver(xy_sta_11(1:15,1), xy_sta_11(1:15,2),...
    out2011.u_opening_GPS'+out2011.u_precursor_GPS',...
    out2011.v_opening_GPS'+out2011.v_precursor_GPS',...
    'k','AutoScale','on','MaxHeadSize',0.5,'LineWidth',1)
% 5 km scale bar
x_qscalebar = [-4.5]; y_qscalebar = [-8.5]; u_qscalebar = [-5]; v_qscalebar = [0];
quiver(x_qscalebar, y_qscalebar, u_qscalebar, v_qscalebar,'k','AutoScale','on','ShowArrowHead','off','LineWidth',2)
text(-8.3,-7.5,'5 km','FontSize',9,'FontWeight','bold')
% north arrow
x_north = [-8.8]; y_north = [8]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(-9.366666666,7.7,'N','FontSize',11)
x_qscale = [8.5]; y_qscale = [-8.5]; u_qscale = [-0.5]; v_qscale = [0];u_qscale2 = [-0.1];
quiver(x_qscale, y_qscale, u_qscale*scale_factor2, v_qscale*scale_factor2,'k',...
    'AutoScale','on','MaxHeadSize',0.5,'LineWidth',1.0)
text(5,-7.5,'0.5 m','FontSize',9)
xlim([-10 10]); ylim([-10 10]);
text(-13.5, 9.5, 'c','FontWeight','bold','FontSize',12);
axis off
title('Extra Basal Slip');
ttt=colorbar('EastOutside'); set(ttt,'YTick',[-0.5, -0.25, 0, 0.25, 0.5]); hold all
set(get(ttt,'ylabel'),'String','[ m ]','FontSize',m-1);
set(ttt, 'Position', [.035+(0.19*scale)+0.005 0.405 0.003 (0.19*ss*scale)-0.02]);
set(gca,'FontName','Avenir'); set(gca,'FontSize',m-1)


axes(axe4); % Basal Opening at Max H-F Opening
hold on;
z=out2011_20kmB_surfC.bed_crack_open_m;
 Nx=24;
 Ny=24; 
 isf=0;
 
 for i=1:Nx
  for j=1:Ny
    isf=isf+1;
    
    x1 = patches_rotate(isf,1);
    x2 = x1 - patchesB(isf,1); % less patch width
    y1 = patches_rotate(isf,2);
    y2 = y1 - patchesB(isf,2); % less patch width
    x=[x1, x2, x2, x1];
    y=[y1, y1, y2, y2];
    patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6]); hold on;
  
  end
 end
colormap(axe4,BWR); 
caxis([-1 1]);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',2)
% 5 km scale bar
x_qscalebar = [-4.5]; y_qscalebar = [-8.5]; u_qscalebar = [-5]; v_qscalebar = [0];
quiver(x_qscalebar, y_qscalebar, u_qscalebar, v_qscalebar,'k','AutoScale','on','ShowArrowHead','off','LineWidth',2)
text(-8.3,-7.5,'5 km','FontSize',9,'FontWeight','bold')
% north arrow
x_north = [-8.8]; y_north = [8]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(-9.366666666,7.7,'N','FontSize',11)
xlim([-10 10]); ylim([-10 10]);
text(-13.5, 9.5, 'd','FontWeight','bold','FontSize',12);
axis off
title('Basal Cavity Opening');
tt1=colorbar('EastOutside'); set(tt1,'YTick',[-1, -0.5, 0, 0.5, 1]); hold all
set(get(tt1,'ylabel'),'String','[ m ]','FontSize',m-1);
set(tt1, 'Position', [.035+(0.19*scale)+0.005 .155 0.003 (0.19*ss*scale)-0.02]);
set(gca,'FontName','Avenir');  set(gca,'FontSize',m-1)

% winter stresses 
axes(axe3); % winter sigma 1
hold on; 
[~,h2]=contourf(xx_fine,yy_fine,(princ_sigma1_winter_vq./1e3),v); set(h2,'LineColor','none'); 
contour(xx_fine,yy_fine,(princ_sigma1_winter_vq./1e3),[v200 v200],'k','LineWidth',1.1);
% north arrow
x_north = [14]; y_north = [-13.25]; u_north = 1.2*[0.1577]; v_north = 1.2*[1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(13.5,-14.05,'N','FontSize',10)
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
xlim([-15 15]); ylim([-15 15]); caxis([-1000 1000])
text(-13.5, 13.5, 'e','FontWeight','bold','FontSize',12);
rectangle('Position',[-10 -9.7 20 20],'EdgeColor',[0.7 0.7 0.7],...
    'LineWidth',0.9,'LineStyle','-');
plot([18.5,8],[-20,20],'--k','LineWidth',0.9); %edge
plot([-13,-20],[-20,5],'--k','LineWidth',0.9); %edge
text(11,13.5,'Edge of Winter Velocities','Rotation',-76,'FontSize',m-1); %edge
% lake names
text(1.25, 5.0, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.95, -3.75, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')
colormap(axe3,BWR); 
cc2 = colorbar('SouthOutside'); 
title('\sigma_{1,winter}','FontSize',m+2,'FontName','Avenir');
set(cc2,'Position',[0.22+0.135 0.075 (tt*2)-0.04 .006],'xtick',[-1000,-800,-600,-400,-200,0,200,400,600,800,1000]);
text(28,-59.5,'\sigma  [ kPa ]','FontSize',m+1,'FontName','Avenir');
set(gca,'xtick',[-10:5:15],'ytick',[-15:5:15],'tickdir','in','LineWidth',1.1,'FontSize',m,'FontName','Avenir','Layer','Top');  box on

axes(axe5); % winter sigma 2
hold on; 
[c2,h2]=contourf(xx_fine,yy_fine,(princ_sigma2_winter_vq./1e3),v); set(h2,'LineColor','none'); 
contour(xx_fine,yy_fine,(princ_sigma2_winter_vq./1e3),[v200 v200],'k','LineWidth',1.1);
% north arrow
x_north = [14]; y_north = [-13.25]; u_north = 1.2*[0.1577]; v_north = 1.2*[1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(13.5,-14.05,'N','FontSize',10)
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
xlim([-15 15]); ylim([-15 15]); caxis([-1000 1000])
% lake names
text(1.25, 5.0, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.95, -3.75, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-13.5, 13.5, 'f','FontWeight','bold','FontSize',12);
rectangle('Position',[-10 -9.7 20 20],'EdgeColor',[0.7 0.7 0.7],...
    'LineWidth',0.9,'LineStyle','-');
plot([18.5,8],[-20,20],'--k','LineWidth',0.9); %edge
plot([-13,-20],[-20,5],'--k','LineWidth',0.9); %edge
colormap(axe5,BWR); 
title('\sigma_{2,winter}','FontSize',m+2,'FontName','Avenir');
xlabel('x [ km ]');
set(gca,'xtick',[-10:5:20],'ytick',[-20:5:20],'tickdir','in','LineWidth',1.1,'FontSize',m,'FontName','Avenir','Layer','Top');  box on

% lake deformation stresses
axes(axe32); % lake sigma 1
hold on; 
i=1280; % 2011 t_{3} index 
[~,h6]=contourf(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1(i,:,:))-princ_sigma1_winter_vq)./1e3,v); 
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200],'k','LineWidth',1.1);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_39mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200],'k--','LineWidth',0.8);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_032mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200],'w','LineWidth',1.5);

% %calculate r_m of 200 kPa contour across \lambda and \mu values
% %dmax 2011 all quadrants
% [Cdmax,~]=contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200]); 
% Cdmax(Cdmax > 15) = NaN;
% NLdmax2011_mu15 = nanmax(sqrt((Cdmax(1,2:end).^2)+(Cdmax(2,2:end).^2)))  
% 
% % % %dmax 2011 inland
% [Cdmax,~]=contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200]); 
% Cdmax(Cdmax > 15) = NaN;
% test = Cdmax(1,:); test(test < 0.5) = NaN;
% NLdmax2011_mu15_inland = nanmax(sqrt((test(1,2:end).^2)+(Cdmax(2,2:end).^2)))  
 
% contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_39mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200],'k--','LineWidth',0.8);
% % %dmax 2011
% [Cdmax,~]=contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_39mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200]); 
% Cdmax(Cdmax > 15) = NaN;
% NLdmax2011_mu39 = nanmax(sqrt((Cdmax(1,2:end).^2)+(Cdmax(2,2:end).^2)))  
% 
% % %dmax 2011 inland 39
% [Cdmax,~]=contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_39mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200]); 
% Cdmax(Cdmax > 15) = NaN;
% test = Cdmax(1,:); test(test < 0.5) = NaN;
% NLdmax2011_mu39_inland = nanmax(sqrt((test(1,2:end).^2)+(Cdmax(2,2:end).^2))) 
 
% contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_032mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200],'w','LineWidth',1.5);
% [Cdmax,~]=contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_032mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200]); 
% Cdmax(Cdmax > 15) = NaN;
% NLdmax2011_mu032 = nanmax(sqrt((Cdmax(1,2:end).^2)+(Cdmax(2,2:end).^2)))
% 
% % %dmax 2011 inland
% [Cdmax,~]=contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_032mu(i,:,:))-princ_sigma1_winter_vq)./1e3,[v200 v200]); 
% Cdmax(Cdmax > 15) = NaN;
% test = Cdmax(1,:); test(test < 0.5) = NaN;
% NLdmax2011_mu032_inland = nanmax(sqrt((test(1,2:end).^2)+(Cdmax(2,2:end).^2))) 
% 
% % save for plotting values on fig for idealized simulations
% NLdmax2011_surfC.mu15 = NLdmax2011_mu15;
% NLdmax2011_surfC.mu39 = NLdmax2011_mu39;
% NLdmax2011_surfC.mu032 = NLdmax2011_mu032;
% NLdmax2011_surfC.mu15_inland = NLdmax2011_mu15_inland;
% NLdmax2011_surfC.mu032_inland = NLdmax2011_mu032_inland;
% NLdmax2011_surfC.mu39_inland = NLdmax2011_mu39_inland;
% %save NLdmax2011_surfC_200kPa_20230721.mat NLdmax2011_surfC

% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
text(-13.5, 13.5, 'g','FontWeight','bold','FontSize',12);
rectangle('Position',[-10 -9.7 20 20],'EdgeColor',[0.7 0.7 0.7],...
    'LineWidth',0.9,'LineStyle','-');
% north arrow
x_north = [14]; y_north = [-13.25]; u_north = 1.2*[0.1577]; v_north = 1.2*[1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(13.5,-14.05,'N','FontSize',10)
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k'); hold on;
xlim([-15 15]); ylim([-15 15]); caxis([-1000 1000])
colormap(axe32,BWR); 
title('\sigma_{1,lake}(2011:{\itt}_{3})','FontSize',m+2,'FontName','Avenir');
set(gca,'xtick',[-10:5:20],'ytick',[-20:5:20],'tickdir','in','LineWidth',1.1,'FontSize',m,...
    'yticklabel',[],'FontName','Avenir','Layer','Top');  box on

axes(axe52); % lake sigma 2
hold on; 
[C5,h6]=contourf(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2(i,:,:))-princ_sigma2_winter_vq)./1e3,v);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2(i,:,:))-princ_sigma2_winter_vq)./1e3,[v200 v200],'k','LineWidth',1.1);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2_39mu(i,:,:))-princ_sigma2_winter_vq)./1e3,[v200 v200],'k--','LineWidth',0.8);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2_032mu(i,:,:))-princ_sigma2_winter_vq)./1e3,[v200 v200],'w','LineWidth',1.5);
% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
text(-13.5, 13.5, 'h','FontWeight','bold','FontSize',12);
rectangle('Position',[-10 -9.7 20 20],'EdgeColor',[0.7 0.7 0.7],...
    'LineWidth',0.9,'LineStyle','-');
% north arrow
x_north = [14]; y_north = [-13.25]; u_north = 1.2*[0.1577]; v_north = 1.2*[1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(13.5,-14.05,'N','FontSize',10)
xlim([-15 15]); ylim([-15 15]); caxis([-1000 1000]);
title('\sigma_{2,lake}(2011:{\itt}_{3})','FontSize',m+2,'FontName','Avenir');
colormap(axe52,BWR); 
set(gca,'xtick',[-10:5:20],'ytick',[-20:5:20],'tickdir','in','LineWidth',1.1,'FontSize',m,'FontName','Avenir','Layer','Top'); 
xlabel('x [ km ]');

axes(axe33); % all four sigma 1
hold on; 
[~,h6]=contourf(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1(i,:,:)))./1e3,v);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1(i,:,:)))./1e3,[v200 v200],'k','LineWidth',1.1);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_39mu(i,:,:)))./1e3,[v200 v200],'k--','LineWidth',0.8);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma1_032mu(i,:,:)))./1e3,[v200 v200],'w','LineWidth',1.5);
% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
text(-13.5, 13.5, 'i','FontWeight','bold','FontSize',12);
plot([18.5,8],[-20,20],'--k','LineWidth',0.9); %edge
plot([-13,-20],[-20,5],'--k','LineWidth',0.9); %edge
rectangle('Position',[-10 -9.7 20 20],'EdgeColor',[0.7 0.7 0.7],...
    'LineWidth',0.9,'LineStyle','-');
% north arrow
x_north = [14]; y_north = [-13.25]; u_north = 1.2*[0.1577]; v_north = 1.2*[1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(13.5,-14.05,'N','FontSize',10)
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k'); hold on;
xlim([-15 15]); ylim([-15 15]); caxis([-1000 1000])
colormap(axe33,BWR); 
title('\sigma_{1}(2011:{\itt}_{3})','FontSize',m+2,'FontName','Avenir'); 
set(gca,'xtick',[-10:5:20],'ytick',[-20:5:20],'tickdir','in','LineWidth',1.1,'FontSize',m,...
    'FontName','Avenir','Layer','Top','yaxislocation','right');  box on
ylabel('y [ km ]');

axes(axe53); % all four sigma 2
hold on;
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2_032mu(i,:,:)))./1e3,[v200 v200],'w','LineWidth',1.8);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2(i,:,:)))./1e3,[v200 v200],'k','LineWidth',1.2);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2_39mu(i,:,:)))./1e3,[v200 v200],'k--','LineWidth',1.2);
[~,h6]=contourf(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2(i,:,:)))./1e3,v);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2(i,:,:)))./1e3,[v200 v200],'k','LineWidth',1.1);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2_39mu(i,:,:)))./1e3,[v200 v200],'k--','LineWidth',0.8);
contour(xx_fine,yy_fine,(squeeze(n2011_allfour.sigma2_032mu(i,:,:)))./1e3,[v200 v200],'w','LineWidth',1.5);
text(-13.5, 13.5, 'j','FontWeight','bold','FontSize',12);
% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
plot([18.5,8],[-20,20],'--k','LineWidth',0.9); %edge
plot([-13,-20],[-20,5],'--k','LineWidth',0.9); %edge
rectangle('Position',[-10 -9.7 20 20],'EdgeColor',[0.7 0.7 0.7],...
    'LineWidth',0.9,'LineStyle','-');
% north arrow
x_north = [14]; y_north = [-13.25]; u_north = 1.2*[0.1577]; v_north = 1.2*[1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
text(13.5,-14.05,'N','FontSize',10)

lgd_i = legend('0.32 GPa','1.5 GPa','3.9 GPa');
legend boxoff
set(lgd_i,'Position',[0.755 0.178 0.06 0.01]);
    
xlim([-15 15]); ylim([-15 15]); caxis([-1000 1000]);
title('\sigma_{2}(2011:{\itt}_{3})','FontSize',m+2,'FontName','Avenir');
colormap(axe53,BWR); 
set(gca,'xtick',[-10:5:20],'ytick',[-20:5:20],'tickdir','in','LineWidth',1.1,...
    'FontSize',m,'FontName','Avenir','Layer','Top','yaxislocation','right'); 
xlabel('x [ km ]'); ylabel('y [ km ]');

% save figure
% figurename=sprintf('paperfigS1_nif2011_max20kmB_200kPa_principals_20240208.png');
% print(gcf,'-dpng','-r500',figurename);
