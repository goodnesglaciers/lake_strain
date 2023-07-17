%% Movie frames for 6 km by 6 km region from pre-saved NIF outputs
% Uses saved NIF L1A outputs for 20 km x 20 km bed plane and vertical 
%  crack along L1A hydro-fracture scarp.
% Laura A. Stevens Oct 26, 2015 (Oxford visit)
% 2021 Mar 04 LAS -- now with strains!
% 2023 May 31 LAS -- updated with Stacy Larochelle infinitesimal strain 
% calcs; finalized for revisions
clear all; close all;

%% 2011 NIF L1A outputs
load time_2011.mat % time vector
% 0.5-km spacing across surface and the NL set-up
[xx_fine,yy_fine] = meshgrid(-20:0.5:20, -20:0.5:20); % [ km ]
xy_surf_fine = horzcat(reshape(xx_fine,81*81,1),reshape(yy_fine,81*81,1)); % [ km ]

% indices for timeseries averages for NNL, NL, and SNL (matrix)
[indNL_vec,indNL2_vec] = find(xx_fine(:,:)<1.75 & xx_fine(:,:)>0 & ...
    yy_fine(:,:)<1.2 & yy_fine(:,:)>-0.2);

[indNNL_vec,indNNL2_vec] = find(xx_fine(:,:)<0.85 & xx_fine(:,:)>-0.55 & ...
    yy_fine(:,:)<4.25 & yy_fine(:,:)>3);

[indSNL_vec,indSNL2_vec] = find(xx_fine(:,:)<0.25 & xx_fine(:,:)>-1.10 & ...
    yy_fine(:,:)<-1.5 & yy_fine(:,:)>-2.4);

% indices for timeseries averages for NNL, NL, and SNL (vector)
[indNL_surf] = find(xy_surf_fine(:,1)<1.75 & xy_surf_fine(:,1)>0 & ...
    xy_surf_fine(:,2)<1.2 & xy_surf_fine(:,2)>-0.2);

[indNNL_surf] = find(xy_surf_fine(:,1)<0.85 & xy_surf_fine(:,1)>-0.55 & ...
    xy_surf_fine(:,2)<4.25 & xy_surf_fine(:,2)>3);

[indSNL_surf] = find(xy_surf_fine(:,1)<0.25 & xy_surf_fine(:,1)>-1.10 & ...
    xy_surf_fine(:,2)<-1.5 & xy_surf_fine(:,2)>-2.4);

% indices for timeseries averages for NNL, NL, and SNL (single point)
[indNL_surf1] = find(xy_surf_fine(:,1)<1.25 & xy_surf_fine(:,1)>0.75 & ...
    xy_surf_fine(:,2)<0.75 & xy_surf_fine(:,2)>0.25);

[indNNL_surf1] = find(xy_surf_fine(:,1)<0.25 & xy_surf_fine(:,1)>-0.25 & ...
    xy_surf_fine(:,2)<3.75 & xy_surf_fine(:,2)>3.25);

[indSNL_surf1] = find(xy_surf_fine(:,1)<-0.35 & xy_surf_fine(:,1)>-0.65 & ...
    xy_surf_fine(:,2)<-1.80 & xy_surf_fine(:,2)>-2.05);

%% load saved surface stresses
load n2011_allfour.mat % all 3 sources + winter = "ALLFOUR"

% calculate \sigma_{2} within lake basins
for i=1:1995 %length(surface_disp_2011.Disp_E(:,1))
    % NL n2011
    n2011_allfour.princ_sigma2_NL(i,1) = nanmean(nanmean(n2011_allfour.sigma2(i,indNL_surf)));
    
    % NNL n2011 
    n2011_allfour.princ_sigma2_NNL(i,1) = nanmean(nanmean(n2011_allfour.sigma2(i,indNNL_surf)));
    
    % SNL n2011
    n2011_allfour.princ_sigma2_SNL(i,1) = nanmean(nanmean(n2011_allfour.sigma2(i,indSNL_surf)));
end

%% winter velocities: load and interpolate to Nsurface locations
load('stresses_nlake_winter_vel_23Sept2020.mat');
origin = [68.72, -49.53]; % M1 moulin
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
    
% principals winter
    princ_sigma1_winter = (0.5.*(Sxx_flow_vec+Syy_flow_vec))+...
        (sqrt(((0.5.*(Sxx_flow_vec-Syy_flow_vec)).^2) + (Sxy_flow_vec.^2)));
    princ_sigma2_winter = (0.5.*(Sxx_flow_vec+Syy_flow_vec))-...
        (sqrt(((0.5.*(Sxx_flow_vec-Syy_flow_vec)).^2) + (Sxy_flow_vec.^2)));
    theta_winter = (0.5.*(atan2((2.*Sxy_flow_vec),(Sxx_flow_vec-Syy_flow_vec))));
    theta_winter_deg = rad2deg(theta_winter)+7; % [deg] rotate 7 degrees to get into 270degrees = -y axis of NIF set-up
    theta_winter_rad = deg2rad(theta_winter_deg); % [rad]

% interp to X Y of Nsurface locations
    F = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),Sxx_flow_vec);
    Sxx_flow_vec_vq = F(xx_fine,yy_fine);
    Sxx_flow_vec_vq_vec = reshape(Sxx_flow_vec_vq,81*81,1); 
    FF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),Syy_flow_vec);
    Syy_flow_vec_vq = FF(xx_fine,yy_fine);
    Syy_flow_vec_vq_vec = reshape(Syy_flow_vec_vq,81*81,1);
    FF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),Sxy_flow_vec);
    Sxy_flow_vec_vq = FF(xx_fine,yy_fine);
    Sxy_flow_vec_vq_vec = reshape(Sxy_flow_vec_vq,81*81,1);
    FFF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),princ_sigma1_winter);
    princ_sigma1_winter_vq = FFF(xx_fine,yy_fine);
    princ_sigma1_winter_vq_vec = reshape(princ_sigma1_winter_vq,81*81,1);
    FFFF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),princ_sigma2_winter);
    princ_sigma2_winter_vq = FFFF(xx_fine,yy_fine);
    FFFFF = scatteredInterpolant(xy_stresses(:,1),xy_stresses(:,2),theta_winter_deg);
    theta_winter_vq_deg = FFFFF(xx_fine,yy_fine); 
    theta_winter_vq_rad = deg2rad(theta_winter_vq_deg);
    
% theta (winter) at lake locations [ degrees ]
theta_NL_winter =  theta_winter_vq_deg(indNL_vec,indNL2_vec); 
theta_NNL_winter =  theta_winter_vq_deg(indNNL_vec,indNNL2_vec); 
theta_SNL_winter =  theta_winter_vq_deg(indSNL_vec,indSNL2_vec);    

%% load L1A geographic files
origin = [68.72, -49.53]; % M1 moulin
load moulin_2011.mat
load out2011_20kmB_third.mat % vertical and bed patches from 20 km NIF L1A
    patchesB = out2011_20kmB_third.patches;
    patchesC = out2011_20kmB_third.patches_C;
load apcoords_lle
    lats=apcoords_lle(1,:); lons=apcoords_lle(2,:); hs=apcoords_lle(3,:);
    llh=[lats; lons; hs];
    xy_sta_11=llh2localxy(llh,origin);
load apcoords_lle_2012
    lats=apcoords_lle_2012(1,:); lons=apcoords_lle_2012(2,:); hs=apcoords_lle_2012(3,:);
    llh=[lats; lons; hs];
    xy_sta_12=llh2localxy(llh,origin);
    llh_moulin = [moulin_2011(:,4)'; moulin_2011(:,3)'; zeros(37,1)'];
    xy_moulin = llh2localxy(llh_moulin,origin);
load lake.mat % L1A
    llh_lake = [lake(:,5)'; lake(:,4)'; zeros(337,1)'];
    xy_lake = llh2localxy(llh_lake,origin);
load lil_lake_2013_168.mat % L1B
    llh_lil_lake = [lil_lake_2013_168(:,4)'; lil_lake_2013_168(:,3)'; zeros(148,1)'];
    xy_lil_lake = llh2localxy(llh_lil_lake,origin);
% polarstereo conversion needed values
    radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
NNL=csvread('NNL_20110617.csv',1,0);  % L1C
    [phi,lambda]=polarstereo_inv(NNL(:,1),NNL(:,2),radius,eccen,lat_true,lon_posy);
    llh_nnl = [phi,lambda]; 
    xy_nnl_lake = llh2localxy(llh_nnl',origin);
% L1C crack from 2016/222
    NNL2016=csvread('nnl2016222_latlon.csv',1,0);
    xy_NNL_crack = llh2localxy(horzcat(NNL2016(:,2), NNL2016(:,1))',origin);
    strike_NNL = 90 - atand((xy_NNL_crack(end,2)-xy_NNL_crack(1,2))./...
        (xy_NNL_crack(end,1)-xy_NNL_crack(1,1)) )  %  for plotting purposes
% L1B crack from 2011 June 17
    SNL=csvread('snl_crack2_2011June17.csv',1,0);
    [phi,lambda]=polarstereo_inv(SNL(:,1),SNL(:,2),radius,eccen,lat_true,lon_posy);
    xy_SNL_crack = llh2localxy([phi,lambda]',origin); 
    strike_SNL = atand((xy_SNL_crack(end,2)-xy_SNL_crack(1,2))./...
        (xy_SNL_crack(end,1)-xy_SNL_crack(1,1)) ) 
% strike L1A crack in lake    
    x_crack = patchesC(1:6:end,6); y_crack = patchesC(1:6:end,7);
    strike_NL = atand((y_crack(22,1)-y_crack(8,1))./(x_crack(22,1)-x_crack(8,1))) +180  % +180 for plotting purposes
    
%%
%%%%%%% plot forward problem surface stresses movie %%%%%%%%%%%%%%%%%%%%%%
for i=130:3:1854 
fig1 = figure('Units','centimeters','Position',1.1.*[30 1 20.5 16]); clf;
    % mapview
    axe1 = axes('Position',[0.099 0.50 0.36 0.46],'Box','on','NextPlot','add');
    axe2 = axes('Position',[0.529 0.50 0.36 0.46],'Box','on','NextPlot','add');
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    % timeseries
    axe3 = axes('Position',[0.099 0.28 0.36 0.14],'Box','on','NextPlot','add');
    axe4 = axes('Position',[0.529 0.28 0.36 0.14],'Box','on','NextPlot','add');
    % rose plots
    axe51 = axes('Position',[0.01 0.0 0.144 0.19],'Box','on','NextPlot','add'); 
    axe52 = axes('Position',[0.17 0.0 0.144 0.19],'Box','on','NextPlot','add');
    axe53 = axes('Position',[0.33 0.0 0.144 0.19],'Box','on','NextPlot','add');
    axe61 = axes('Position',[0.52 0.0 0.145 0.19],'Box','on','NextPlot','add'); 
    axe62 = axes('Position',[0.68 0.0 0.145 0.19],'Box','on','NextPlot','add'); 
    axe63 = axes('Position',[0.84 0.0 0.145 0.19],'Box','on','NextPlot','add'); 

    load BWR.mat % colormap
    vvKPA=[-8005:10:8005]; % kPa contours
    v200=200; % 200 kPa theshold
    TriangleSize = 7; % marker size
    % dock at eel pond 
    sky_blue = [111, 169, 228]./255; % SNL
    metal = [87, 115, 131]./255; % NNL
    oar = [251, 219, 154]./255; % NENL
    handle = [161, 37, 49]./255; % NL
    dark_oar = [164, 114, 63]./255;

axes(axe1) % principle stress \sigma_{1} -- all 3 sources + winter = ALLFOUR
hold on;
[C5,h6]=contourf(xx_fine,yy_fine,squeeze(n2011_allfour.sigma1(i,:,:)./1e3),vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,squeeze(n2011_allfour.sigma1(i,:,:)./1e3),[v200 v200],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% boxes single point
rectangle('position',[-0.25 3.75 0.5 0.5],'EdgeColor','w','FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[0.75 0.25 0.5 0.5],'EdgeColor','w','FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-0.75 -2.25 0.5 0.5],'EdgeColor','w','FaceColor','none','LineWidth',1.5); % L1B
% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)
% tick marks 
        hypot = 0.2;
% lake + winter P2
        theta_15GPa_P2_plotX1(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2011_allfour.theta(i,:))+90),81*81,1));  
        theta_15GPa_P2_plotX2(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2011_allfour.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotY1(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2011_allfour.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotY2(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2011_allfour.theta(i,:))+90),81*81,1)); 
%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2(1,j) theta_15GPa_P2_plotX1(1,j)],...
        [theta_15GPa_P2_plotY2(1,j) theta_15GPa_P2_plotY1(1,j)],...
        '-','LineWidth',1.1,'Color',[0.6 0.6 0.6]);   
    end 
caxis([-1000 1000]);
colormap(axe1,BWR); 

% lake names
text(1.25, 4.5, 'L1C','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',13,'FontName','Avenir','FontAngle','italic')

text(-8.75, 5.5, 'a','FontSize',10,'FontWeight','bold');
text(-9.1, 6.65, '{\itStevens et al.} (2023) Movie M1','FontSize',9,'FontName','Avenir');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10); 
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('[\sigma,\theta]_{1}','FontSize',10,'FontName','Avenir')
% time string top in middle
text(5.5,6.6,sprintf('2011/%6.3f',time_2011(i)),'FontWeight','bold','FontSize',10);

axes(axe2) % principle stress \sigma_{2} -- all 3 sources + winter = ALLFOUR
hold on;
[C2,h2]=contourf(xx_fine,yy_fine,squeeze(n2011_allfour.sigma2(i,:,:)./1e3),vvKPA);
set(h2,'LineColor','none'); 
contour(xx_fine,yy_fine,squeeze(n2011_allfour.sigma2(i,:,:)./1e3),[v200 v200],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% boxes single point
rectangle('position',[-0.25 3.75 0.5 0.5],'EdgeColor','w','FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[0.75 0.25 0.5 0.5],'EdgeColor','w','FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-0.75 -2.25 0.5 0.5],'EdgeColor','w','FaceColor','none','LineWidth',1.5); % L1B
% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)
% tick marks 
        hypot = 0.2;
% lake + winter P2
        theta_15GPa_P2_plotX1(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2011_allfour.theta(i,:))),81*81,1));  
        theta_15GPa_P2_plotX2(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2011_allfour.theta(i,:))),81*81,1)); 
        theta_15GPa_P2_plotY1(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2011_allfour.theta(i,:))),81*81,1)); 
        theta_15GPa_P2_plotY2(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2011_allfour.theta(i,:))),81*81,1)); 
%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2(1,j) theta_15GPa_P2_plotX1(1,j)],...
        [theta_15GPa_P2_plotY2(1,j) theta_15GPa_P2_plotY1(1,j)],...
        '-','LineWidth',1.1,'Color',[0.6 0.6 0.6]);   
    end 
caxis([-1000 1000]);
colormap(axe2,BWR); 

% lake names
text(1.25, 4.5, 'L1C','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',13,'FontName','Avenir','FontAngle','italic')

text(-7.25, 5.5, 'b','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('[\sigma,\theta]_{2}','FontSize',10,'FontName','Avenir')

t2=colorbar('EastOutside');
set(t2,'YTick',[-1500,-1000,-500,0,500,1000,1500],'TickDirection','out'); 
hold all
set(get(t2,'ylabel'),'String','\sigma  [ kPa ]','FontSize',10);
set(t2, 'Position', [.91 0.51 .007 (0.46)-0.02]);

%% timeseries

axes(axe3) % \sigma_{1} timeseries
plot(time_2011(1:3:1855),(n2011_allfour.princ_sigma1_NNL(1:3:1855))./1e3,'-','LineWidth',1.2,'Color',metal); 
plot(time_2011(1:3:1855),(n2011_allfour.princ_sigma1_NL(1:3:1855))./1e3,'-','LineWidth',1.2,'Color',handle); 
plot(time_2011(1:3:1855),(n2011_allfour.princ_sigma1_SNL(1:3:1855))./1e3,'-','LineWidth',1.2,'Color',sky_blue); 

plot([168.85 168.85],[-499 999],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-499 999],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-499 999],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32+(6/24) 169.32+(6/24)],[-499 999],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); % maxwell

text(168.85-0.05,1200,'{\itt}_{1}','FontSize',7,'FontName','Avenir'); 
text(169.21-0.05,1200,'{\itt}_{2}','FontSize',7,'FontName','Avenir'); 
text(169.32-0.05,1200,'{\itt}_{3}','FontSize',7,'FontName','Avenir'); 
text(170.32-0.05,1200,'{\itt}_{4}','FontSize',8,'FontName','Avenir'); 
text(169.32-0.05+(4/24),1200,'{\it\tau} = {6 hr}','FontSize',8,'FontName','Avenir','Color',[0.2 0.55 0.2]); % maxwell

plot(time_2011(1:1995),(n2011_allfour.princ_sigma1_NNL(1:1995))./1e3,'-','LineWidth',1.2,'Color',metal); 
plot(time_2011(1:1995),(n2011_allfour.princ_sigma1_SNL(1:1995))./1e3,'-','LineWidth',1.2,'Color',sky_blue); 
plot(time_2011(1:1995),(n2011_allfour.princ_sigma1_NL(1:1995))./1e3,'-','LineWidth',1.2,'Color',handle); 

% moving time marker
plot([time_2011(i) time_2011(i)],[-499 999],'-','LineWidth',1.05,'Color',[0.4940 0.1840 0.5560]);
plot(time_2011(i),1000,'v','MarkerSize',TriangleSize-3,'MarkerFaceColor',[0.4940 0.1840 0.5560],'Color',[0.4940 0.1840 0.5560]);

xlabel('Day of Year, 2011'); ylabel('\sigma_{1} [ kPa ]');
lgd_s1=legend('L1C','L1A','L1B');
legend boxoff
set(lgd_s1,'Position',[0.1 0.36 0.1 0.03]);
ylim([-500 1000]); xlim([167.32 170.32]);
text(166.65, 1200, 'c','FontSize',10,'FontWeight','bold');
text(166.7, -1400, 'e','FontSize',10,'FontWeight','bold');
set(gca,'tickdir','in','LineWidth',1.05,'FontSize',9); 
set(gca,'FontName','Avenir','ytick',[-500,0,500,1000]);
grid on


axes(axe4) % \sigma_{2} timeseries
plot(time_2011(1:3:1855),(n2011_allfour.princ_sigma2_NNL(1:3:1855))./1e3,'-','LineWidth',1.2,'Color',metal); 
plot(time_2011(1:3:1855),(n2011_allfour.princ_sigma2_NL(1:3:1855))./1e3,'-','LineWidth',1.2,'Color',handle); 
plot(time_2011(1:3:1855),(n2011_allfour.princ_sigma2_SNL(1:3:1855))./1e3,'-','LineWidth',1.2,'Color',sky_blue); 

plot([168.85 168.85],[-499 999],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-499 999],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-499 999],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32+(6/24) 169.32+(6/24)],[-499 999],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); % maxwell

text(168.85-0.05,1200,'{\itt}_{1}','FontSize',7,'FontName','Avenir'); 
text(169.21-0.05,1200,'{\itt}_{2}','FontSize',7,'FontName','Avenir'); 
text(169.32-0.05,1200,'{\itt}_{3}','FontSize',7,'FontName','Avenir'); 
text(170.32-0.05,1200,'{\itt}_{4}','FontSize',8,'FontName','Avenir'); 
text(169.32-0.05+(4/24),1200,'{\it\tau} = {6 hr}','FontSize',8,'FontName','Avenir','Color',[0.2 0.55 0.2]); % maxwell

plot(time_2011(1:1995),(n2011_allfour.princ_sigma2_NNL(1:1995))./1e3,'-','LineWidth',1.2,'Color',metal); 
plot(time_2011(1:1995),(n2011_allfour.princ_sigma2_SNL(1:1995))./1e3,'-','LineWidth',1.2,'Color',sky_blue); 
plot(time_2011(1:1995),(n2011_allfour.princ_sigma2_NL(1:1995))./1e3,'-','LineWidth',1.2,'Color',handle); 

% moving time marker
plot([time_2011(i) time_2011(i)],[-499 999],'-','LineWidth',1.05,'Color',[0.4940 0.1840 0.5560]);
plot(time_2011(i),1000,'v','MarkerSize',TriangleSize-3,'MarkerFaceColor',[0.4940 0.1840 0.5560],'Color',[0.4940 0.1840 0.5560]);

xlabel('Day of Year, 2011'); ylabel('\sigma_{2} [ kPa ]');
lgd_s2=legend('L1C','L1A','L1B');
legend boxoff
set(lgd_s2,'Position',[0.53 0.36 0.1 0.03]);
set(gca,'tickdir','in','LineWidth',1.05,'FontSize',9); 
set(gca, 'YAxisLocation', 'right')
set(gca,'FontName','Avenir','ytick',[-500,0,500,1000]);
ylim([-500 1000]); xlim([167.32 170.32]);
text(167.05, 1200, 'd','FontSize',10,'FontWeight','bold');
text(167.1, -1400, 'f','FontSize',10,'FontWeight','bold');
grid on 


%% polar histograms
polar_cmap = BWR;
cmin = -1000; cmax = 1000;
m = length(polar_cmap); 
 
polaraxes(axe51) % polarhistogram is in radians -- L1C
hold off
polarhistogram(deg2rad(strike_NNL),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NNL+180),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NNL_winter(2,3)+90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NNL_winter(2,3)-90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% % call values; select location on HF; assign color; 
tt_NNL = rad2deg(n2011_allfour.theta(i,indNNL_surf1))+90;
color = round(n2011_allfour.princ_sigma1_NNL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_NNL),'BinWidth',pi./23,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_NNL+180),'BinWidth',pi./23,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 180]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('\theta_{1}(L1C)','FontWeight','bold','FontSize',10)


polaraxes(axe52) % polarhistogram is in radians -- L1A
hold off
polarhistogram(deg2rad(strike_NL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NL_winter(1,4)+90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NL_winter(1,4)-90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% call values; assign color
tt = rad2deg(n2011_allfour.theta(i,indNL_surf1))+90;
color = round(n2011_allfour.princ_sigma1_NL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt+180),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 180]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('\theta_{1}(L1A)','FontWeight','bold','FontSize',10)

    
polaraxes(axe53) % polarhistogram is in radians - L1B
hold off
polarhistogram(deg2rad(strike_SNL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_SNL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_SNL_winter(1,2)+90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_SNL_winter(1,2)-90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% call values; assign color
tt_SNL = rad2deg(n2011_allfour.theta(i,indSNL_surf1))+90;
color = round(n2011_allfour.princ_sigma1_SNL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_SNL),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_SNL+180),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 180]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('\theta_{1}(L1B)','FontWeight','bold','FontSize',10)


polaraxes(axe61) % polarhistogram is in radians -- L1C
hold off
polarhistogram(deg2rad(strike_NNL),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NNL+180),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NNL_winter(2,3)),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NNL_winter(2,3)+180),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% % call values; select location on HF; assign color; 
tt_NNL = rad2deg(n2011_allfour.theta(i,indNNL_surf1));
color = round(n2011_allfour.princ_sigma2_NNL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_NNL),'BinWidth',pi./23,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_NNL+180),'BinWidth',pi./23,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 180]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('\theta_{2}(L1C)','FontWeight','bold','FontSize',10)


polaraxes(axe62) % polarhistogram is in radians -- L1A
hold off
polarhistogram(deg2rad(strike_NL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NL_winter(1,4)),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NL_winter(1,4)+180),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% call values; assign color
tt = rad2deg(n2011_allfour.theta(i,indNL_surf1));
color = round(n2011_allfour.princ_sigma2_NL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt+180),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 180]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('\theta_{2}(L1A)','FontWeight','bold','FontSize',10)

    
polaraxes(axe63) % polarhistogram is in radians - L1B
hold off
polarhistogram(deg2rad(strike_SNL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_SNL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','y','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_SNL_winter(1,2)),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_SNL_winter(1,2)+180),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% call values; assign color
tt_SNL = rad2deg(n2011_allfour.theta(i,indSNL_surf1));
color = round(n2011_allfour.princ_sigma2_SNL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_SNL),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_SNL+180),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 180]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('\theta_{2}(L1B)','FontWeight','bold','FontSize',10)


% print movie frame
rez=500; %resolution (dpi) of final graphic
figpos=getpixelposition(fig1); %dont need to change anything here
resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
set(fig1,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
path='Animation2011_dials_2023July'; % folder path
name=sprintf('img%04d.png',i-128); % file name
print(gcf,'-dpng','-r500',fullfile(path,name))
close(gcf)

end

%% concatenate .png in ffmpeg to make movie
%ffmpeg -r 5 -f image2 -s 4439x3465 -pattern_type glob -i "*.png" -vcodec libx264 -crf 25 -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p output.mp4

