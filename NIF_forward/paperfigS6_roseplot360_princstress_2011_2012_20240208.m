%% Roseplot for surface stress change direction
% 26 Feb 2021 LAS
% 16 March 2023 LAS -- only plot theta_1 for one location on the HFs
% 19 June 2023 LAS -- Include S. Larochelle stress derivation, read out the
% directions of the stress change 
clear all; close all;

%% 2011
load('../NIF_L1A/Gsurface500_strain_surfacecrack.mat') % Green functions for 20 km X 20 km surface using
% 0.5-km spacing across surface and the NL set-up

%% load north lake geographic files
origin = [68.72, -49.53]; % M1 moulin
load moulin_2011.mat
load out2011_20kmB_third.mat
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
    NNL=csvread('NNL_20110617.csv',1,0); 
    [phi,lambda]=polarstereo_inv(NNL(:,1),NNL(:,2),radius,eccen,lat_true,lon_posy);
    llh_nnl = [phi,lambda]; 
    xy_nnl_lake = llh2localxy(llh_nnl',origin); 
% L1C crack from 2016/222
    NNL2016=csvread('nnl2016222_latlon.csv',1,0);
    xy_NNL_crack = llh2localxy(horzcat(NNL2016(:,2), NNL2016(:,1))',origin);
    strike_NNL = 90 - atand((xy_NNL_crack(end,2)-xy_NNL_crack(1,2))./...
        (xy_NNL_crack(end,1)-xy_NNL_crack(1,1)) ) ; %  for plotting purposes
% L1B crack from 2011 June 17
    SNL=csvread('snl_crack2_2011June17.csv',1,0);
    [phi,lambda]=polarstereo_inv(SNL(:,1),SNL(:,2),radius,eccen,lat_true,lon_posy);
    xy_SNL_crack = llh2localxy([phi,lambda]',origin); 
    strike_SNL = atand((xy_SNL_crack(end,2)-xy_SNL_crack(1,2))./...
        (xy_SNL_crack(end,1)-xy_SNL_crack(1,1)) ) ;
% strike L1A crack in lake    
    x_crack = patchesC(1:6:end,6); y_crack = patchesC(1:6:end,7);
    strike_NL = atand((y_crack(22,1)-y_crack(8,1))./(x_crack(22,1)-x_crack(8,1))) +180 ; % +180 for plotting purposes

 TriangleSize = 7;

%% time and spatial info
load time_2011.mat % time vector
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

% indices for timeseries averages for NNL, NL, and SNL (vector, single point)
[indNL_surf1] = find(xy_surf_fine(:,1)<1.25 & xy_surf_fine(:,1)>0.75 & ...
    xy_surf_fine(:,2)<0.75 & xy_surf_fine(:,2)>0.25);

[indNNL_surf1] = find(xy_surf_fine(:,1)<0.25 & xy_surf_fine(:,1)>-0.25 & ...
    xy_surf_fine(:,2)<3.75 & xy_surf_fine(:,2)>3.25);

[indSNL_surf1] = find(xy_surf_fine(:,1)<-0.35 & xy_surf_fine(:,1)>-0.65 & ...
    xy_surf_fine(:,2)<-1.80 & xy_surf_fine(:,2)>-2.05);

%% winter velocities: load and interpolate to Nsurface locations
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

%% load pre-calculated all-four component stresses
load n2011_allfour_surfC.mat
load n2012_allfour_surfC.mat

% theta at t_3 at lake locations (radians)
theta_2011_mat = squeeze(n2011_allfour.theta_mat(1280,:,:));
theta_2011_NNL = theta_2011_mat(indNNL_vec,indNNL2_vec);
theta_2011_NL = theta_2011_mat(indNL_vec,indNL2_vec);
theta_2011_SNL = theta_2011_mat(indSNL_vec,indSNL2_vec);

theta_2012_mat = squeeze(n2012_allfour.theta_mat(544,:,:));
theta_2012_NNL = theta_2012_mat(indNNL_vec,indNNL2_vec);
theta_2012_NL = theta_2012_mat(indNL_vec,indNL2_vec);
theta_2012_SNL = theta_2012_mat(indSNL_vec,indSNL2_vec);

%% plot surface stresses and rose plots
load BWR.mat 
    
%%%%%%% 2011 %%%%%%% 
for i=1280 % 2011 t_3 index

    fig1 = figure('Units','centimeters','Position',1.0.*[2 1 25 15]);
    clf
    axe1 = axes('Position',[0.045 0.50 0.264 0.44],'Box','on','NextPlot','add'); 
    
    axe2 = axes('Position',[0.345 0.50 0.264 0.44],'Box','on','NextPlot','add');
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe3 = axes('Position',[0.645 0.50 0.264 0.44],'Box','on','NextPlot','add');
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe51 = axes('Position',[0.1 0.12 0.145 0.19],'Box','on','NextPlot','add'); 
    axe52 = axes('Position',[0.25 0.12 0.145 0.19],'Box','on','NextPlot','add');
    axe53 = axes('Position',[0.4 0.12 0.145 0.19],'Box','on','NextPlot','add');
    
    axe61 = axes('Position',[0.55 0.12 0.145 0.19],'Box','on','NextPlot','add'); 
    axe62 = axes('Position',[0.7 0.12 0.145 0.19],'Box','on','NextPlot','add'); 
    axe63 = axes('Position',[0.85 0.12 0.145 0.19],'Box','on','NextPlot','add'); 
       
% dock at eel pond 
sky_blue = [111, 169, 228]./255; % SNL
metal = [87, 115, 131]./255; % NNL
oar = [251, 219, 154]./255; % NENL
handle = [161, 37, 49]./255; % NL
dark_oar = [164, 114, 63]./255;

% principle stress 1 -- winter
axes(axe1)
vvKPA=[-2005:10:2005]; v150=200;
[C5,h6]=contourf(xx_fine,yy_fine,princ_sigma1_winter_vq./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,princ_sigma1_winter_vq./1e3,[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.25 3.75 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %NNL
rectangle('position',[0.75 0.25 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %NL
rectangle('position',[-0.75 -2.25 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %SNL
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k','LineWidth',1.5)
% tick marks winter
        hypot = 0.2;
% just winter P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(theta_winter_vq_deg(:,:)+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(theta_winter_vq_deg(:,:)+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(theta_winter_vq_deg(:,:)+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(theta_winter_vq_deg(:,:)+90),81*81,1));

%theta_P tick marks winter
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.1,'Color',[0.5 0.5 0.5]);   
    end    
caxis([-1000 1000]);
colormap(axe1,BWR); 
% lake names
text(1.25, 4.5, 'L1C','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',13,'FontName','Avenir','FontAngle','italic')
text(-5.5, 5.5, 'a','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',8,'FontName','Avenir'); 
xlabel(' x [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',1.01,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;
t2=colorbar('EastOutside');
set(t2,'YTick',[-1000,-500,0,500,1000],'TickDirection','out','FontSize',9); 
hold all
set(get(t2,'ylabel'),'String','\sigma_{1}  [ kPa ]','FontSize',10);
set(t2, 'Position', [0.924 0.50 .006 0.44]);
title('[\sigma,\theta]_{1,winter}','FontSize',9,'FontName','Avenir')

% principle stress 1 -- 2011
axes(axe2)
vvKPA=[-7005:10:7005]; v150=200;
[C5,h6]=contourf(xx_fine,yy_fine,squeeze(n2011_allfour.sigma1(i,:,:)./1e3),vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,squeeze(n2011_allfour.sigma1(i,:,:)./1e3),[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.25 3.75 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %NNL
rectangle('position',[0.75 0.25 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %NL
rectangle('position',[-0.75 -2.25 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %SNL
% lakes
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k','LineWidth',1.5)
% tick marks 
        hypot = 0.2;
% lake + winter P2
        theta_15GPa_P2_plotX1(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2011_allfour.theta(i,:,:))+90),81*81,1));  
        theta_15GPa_P2_plotX2(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2011_allfour.theta(i,:,:))+90),81*81,1)); 
        theta_15GPa_P2_plotY1(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2011_allfour.theta(i,:,:))+90),81*81,1)); 
        theta_15GPa_P2_plotY2(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2011_allfour.theta(i,:,:))+90),81*81,1)); 
%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2(1,j) theta_15GPa_P2_plotX1(1,j)],...
        [theta_15GPa_P2_plotY2(1,j) theta_15GPa_P2_plotY1(1,j)],...
        '-','LineWidth',1.1,'Color',[0.5 0.5 0.5]);   
    end    
caxis([-1000 1000]);
colormap(axe2,BWR); 
text(-5.5, 5.5, 'b','FontSize',10,'FontWeight','bold');
xlabel(' x [ km ]','FontSize',10,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir','Layer','Top'); 
grid on;
title('[\sigma,\theta]_{1}(2011:{\itt}_{3})','FontSize',11,'FontName','Avenir')

% polar histograms
polar_cmap = BWR;
cmin = -1000; cmax = 1000;
m = length(polar_cmap); 
 
% NNL P1 (L1C)
polaraxes(axe53) % polarhistogram is in radians
hold off
polarhistogram(deg2rad(strike_NNL),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NNL+180),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
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
thetalim([0 360]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
%title({'L1C';'\theta_{1}(2011:{\itt}_{3})'},'FontWeight','bold','FontSize',10)
title('L1C 2011','FontWeight','bold','FontSize',10)

% NL P1 (L1A)
polaraxes(axe51) % polarhistogram is in radians
hold off
polarhistogram(deg2rad(strike_NL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
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
thetalim([0 360]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('L1A 2011','FontWeight','bold','FontSize',10)
    
% SNL P1 2011 (L1B)
polaraxes(axe52) % polarhistogram is in radians
hold off
polarhistogram(deg2rad(strike_SNL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_SNL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
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
thetalim([0 360]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('L1B 2011','FontWeight','bold','FontSize',10)

% read out the directions of winter and t_3 theta_2 (most compressive)
winter_theta2 = [theta_NL_winter(1,4)+90;theta_SNL_winter(1,2)+90;theta_NNL_winter(2,3)+90]
n2011_t3_theta2 = [tt;tt_SNL;tt_NNL]
end

%% 2012
%%%%%%% 2012 %%%%%%%
for i = 544 % 2012 t_3 index
    
% principle stress 1 -- 2012
axes(axe3)
vvKPA=[-7005:10:7005]; v150=200;
[C5,h6]=contourf(xx_fine,yy_fine,squeeze(n2012_allfour.sigma1(i,:,:)./1e3),vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,squeeze(n2012_allfour.sigma1(i,:,:)./1e3),[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.25 3.75 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %NNL
rectangle('position',[0.75 0.25 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %NL
rectangle('position',[-0.75 -2.25 0.5 0.5],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); %SNL
% lakes
plot(xy_sta_12(:,1),xy_sta_12(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k','LineWidth',1.5)
% tick marks 
        hypot = 0.2;
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2012_allfour.theta(i,:,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2012_allfour.theta(i,:,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2012_allfour.theta(i,:,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2012_allfour.theta(i,:,:))+90),81*81,1));
%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.1,'Color',[0.5 0.5 0.5]);    
    end   
caxis([-1000 1000]);
colormap(axe3,BWR); 
text(-5.5, 5.5, 'c','FontSize',10,'FontWeight','bold'); 
xlabel(' x [ km ]','FontSize',10,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir','Layer','Top'); 
grid on;

title('[\sigma,\theta]_{1}(2012:{\itt}_{3})','FontSize',11,'FontName','Avenir')

% polar histograms
polar_cmap = BWR;
cmin = -1000; cmax = 1000;
m = length(polar_cmap); 

% NNL P1 2012
polaraxes(axe63) % polarhistogram is in radians
hold off
polarhistogram(deg2rad(strike_NNL),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NNL+180),'BinWidth',pi./25,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NNL_winter(2,3)+90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NNL_winter(2,3)-90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);

% % call values; select location on HF; assign color; 
tt_NNL_2012 = rad2deg(n2012_allfour.theta(i,indNNL_surf1))+90;
color = round(n2012_allfour.princ_sigma1_NNL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_NNL_2012),'BinWidth',pi./23,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_NNL_2012+180),'BinWidth',pi./23,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 360]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('L1C 2012','FontWeight','bold','FontSize',10)

% NL P1 2012
polaraxes(axe61) % polarhistogram is in radians
hold off
polarhistogram(deg2rad(strike_NL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_NL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NL_winter(1,4)+90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_NL_winter(1,4)-90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% call values; assign color
tt_2012 = rad2deg(n2012_allfour.theta(i,indNL_surf1))+90;
color = round(n2012_allfour.princ_sigma1_NL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_2012),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_2012+180),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 360]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('L1A 2012','FontWeight','bold','FontSize',10)

% SNL P1 2012
polaraxes(axe62) % polarhistogram is in radians
hold off
polarhistogram(deg2rad(strike_SNL),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
hold on;
polarhistogram(deg2rad(strike_SNL+180),'BinWidth',pi./30,'Normalization','probability',...
    'FaceColor','k','FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_SNL_winter(1,2)+90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
polarhistogram(deg2rad(theta_SNL_winter(1,2)-90),'BinWidth',pi./33,'Normalization','probability',...
    'FaceColor',[0.6 0.6 0.6],'FaceAlpha',.7,'EdgeColor',[0.1 0.1 0.1]);
% call values; assign color
tt_SNL_2012 = rad2deg(n2012_allfour.theta(i,indSNL_surf1))+90;
color = round(n2012_allfour.princ_sigma1_SNL(i)./1e3);
index = fix((color-cmin)/(cmax-cmin)*m)+1;
RGB = ind2rgb(index,polar_cmap);
% plot
polarhistogram(deg2rad(tt_SNL_2012),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
polarhistogram(deg2rad(tt_SNL_2012+180),'BinWidth',pi./30,'Normalization','probability',...
        'FaceColor',RGB,'FaceAlpha',0.9,'EdgeColor',[0.4 0.4 0.4]);
thetalim([0 360]); rlim([0 0.3]);
set(gca,'FontName','Avenir','FontSize',8,'LineWidth',1.05,'ThetaTickLabel',[],...
    'Rtick',[0, 0.3],'Rticklabel',[]);
title('L1B 2012','FontWeight','bold','FontSize',10)

% read out the directions of t_3 theta_2 (most compressive)
n2012_t3_theta2 = [tt_2012;tt_SNL_2012;tt_NNL_2012]

end

%% print figure
%print(gcf,'-dpng','-r500','paperfigS6_roseplot360_20240208.png')