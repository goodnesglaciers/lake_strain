%% Nsurface principle stress source decomposition 2011, 2012
% calculates and saves n2011_allfour.mat and n2012_allfour.mat
% 14 June 2021 LAS -- decomposition into three NIF sources
% 24 March 2022 LAS -- decomposition only the NIF sources (no winter viscous)
% 19 June 2023 LAS -- Include S. Larochelle stress derivation; plotting
% updates from revisions; surface crack
clear all; close all;

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
        (xy_NNL_crack(end,1)-xy_NNL_crack(1,1)) ) ; %  for plotting purposes
% L1B crack from 2011 June 17
    SNL=csvread('snl_crack2_2011June17.csv',1,0);
    [phi,lambda]=polarstereo_inv(SNL(:,1),SNL(:,2),radius,eccen,lat_true,lon_posy);
    xy_SNL_crack = llh2localxy([phi,lambda]',origin); 
    strike_SNL = atand((xy_SNL_crack(end,2)-xy_SNL_crack(1,2))./...
        (xy_SNL_crack(end,1)-xy_SNL_crack(1,1)) ) ;
% strike L1A crack in lake    
    x_crack = patchesC(1:6:end,6); y_crack = patchesC(1:6:end,7);
    strike_NL = atand((y_crack(22,1)-y_crack(8,1))./(x_crack(22,1)-x_crack(8,1))) +180;  % +180 for plotting purposes

%% time and spatial info from 2011 NIF output
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

%% hydro-fracture: load stresses, calculate mean within lake basins 
load sigma_ij_2011_SLder_HF_surfC_032mu.mat
n2011_HF_032mu = sigma_ij_2011_SLder_HF;
load sigma_ij_2011_SLder_HF_surfC_39mu.mat
n2011_HF_39mu = sigma_ij_2011_SLder_HF;
load sigma_ij_2011_SLder_HF_surfC.mat
n2011_HF = sigma_ij_2011_SLder_HF;

for i=600:1995 %length(surface_disp_2011.Disp_E(:,1))
% NL n2011_HF. 
n2011_HF.princ_sigma1_NL(i,1) = nanmean(nanmean(n2011_HF.sigma1(i,indNL_surf)));
n2011_HF.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2011_HF_032mu.sigma1(i,indNL_surf)));
n2011_HF.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2011_HF_39mu.sigma1(i,indNL_surf)));

% NNL n2011_HF. 
n2011_HF.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2011_HF.sigma1(i,indNNL_surf)));
n2011_HF.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2011_HF_032mu.sigma1(i,indNNL_surf)));
n2011_HF.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2011_HF_39mu.sigma1(i,indNNL_surf)));

% SNL n2011_HF.  
n2011_HF.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2011_HF.sigma1(i,indSNL_surf)));
n2011_HF.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2011_HF_032mu.sigma1(i,indSNL_surf)));
n2011_HF.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2011_HF_39mu.sigma1(i,indSNL_surf)));
end
 
%% basal opening: load stresses, calculate mean within lake basins
load sigma_ij_2011_SLder_bedopen_surfC_032mu.mat
n2011_bedopen_032mu = sigma_ij_2011_SLder_bedopen;
load sigma_ij_2011_SLder_bedopen_surfC_39mu.mat
n2011_bedopen_39mu = sigma_ij_2011_SLder_bedopen;
load sigma_ij_2011_SLder_bedopen_surfC.mat
n2011_bedopen = sigma_ij_2011_SLder_bedopen;

for i=600:1995 %length(surface_disp_2011.Disp_E(:,1))
% NL n2011_bedopen. 
n2011_bedopen.princ_sigma1_NL(i,1) = nanmean(nanmean(n2011_bedopen.sigma1(i,indNL_surf)));
n2011_bedopen.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2011_bedopen_032mu.sigma1(i,indNL_surf)));
n2011_bedopen.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2011_bedopen_39mu.sigma1(i,indNL_surf)));

% NNL n2011_bedopen. 
n2011_bedopen.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2011_bedopen.sigma1(i,indNNL_surf)));
n2011_bedopen.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2011_bedopen_032mu.sigma1(i,indNNL_surf)));
n2011_bedopen.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2011_bedopen_39mu.sigma1(i,indNNL_surf)));

% SNL n2011_bedopen.  
n2011_bedopen.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2011_bedopen.sigma1(i,indSNL_surf)));
n2011_bedopen.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2011_bedopen_032mu.sigma1(i,indSNL_surf)));
n2011_bedopen.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2011_bedopen_39mu.sigma1(i,indSNL_surf)));
end

%% basal slip: load stresses, calculate mean within lake basins
load sigma_ij_2011_SLder_bedslip_surfC_032mu.mat
n2011_bedslip_032mu = sigma_ij_2011_SLder_bedslip;
load sigma_ij_2011_SLder_bedslip_surfC_39mu.mat
n2011_bedslip_39mu = sigma_ij_2011_SLder_bedslip;
load sigma_ij_2011_SLder_bedslip_surfC.mat
n2011_bedslip = sigma_ij_2011_SLder_bedslip;

for i=600:1995 %length(surface_disp_2011.Disp_E(:,1))
% NL n2011_bedslip. 
n2011_bedslip.princ_sigma1_NL(i,1) = nanmean(nanmean(n2011_bedslip.sigma1(i,indNL_surf)));
n2011_bedslip.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2011_bedslip_032mu.sigma1(i,indNL_surf)));
n2011_bedslip.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2011_bedslip_39mu.sigma1(i,indNL_surf)));

% NNL n2011_bedslip. 
n2011_bedslip.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2011_bedslip.sigma1(i,indNNL_surf)));
n2011_bedslip.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2011_bedslip_032mu.sigma1(i,indNNL_surf)));
n2011_bedslip.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2011_bedslip_39mu.sigma1(i,indNNL_surf)));

% SNL n2011_bedslip.  
n2011_bedslip.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2011_bedslip.sigma1(i,indSNL_surf)));
n2011_bedslip.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2011_bedslip_032mu.sigma1(i,indSNL_surf)));
n2011_bedslip.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2011_bedslip_39mu.sigma1(i,indSNL_surf)));
end

%% winter+HF+bed_open+basal_slip: calculate sigma1 + theta for all four contributions

for i=1:length(time_2011)
    
n2011_allfour.sigma1(i,:,:) = reshape(n2011_HF.sigma1(i,:)+n2011_bedopen.sigma1(i,:)+n2011_bedslip.sigma1(i,:),81,81)+princ_sigma1_winter_vq;
n2011_allfour.sigma2(i,:,:) = reshape(n2011_HF.sigma2(i,:)+n2011_bedopen.sigma2(i,:)+n2011_bedslip.sigma2(i,:),81,81)+princ_sigma2_winter_vq;
n2011_allfour.sigma1_032mu(i,:,:) = reshape(n2011_HF_032mu.sigma1(i,:)+n2011_bedopen_032mu.sigma1(i,:)+n2011_bedslip_032mu.sigma1(i,:),81,81)+princ_sigma1_winter_vq;
n2011_allfour.sigma2_032mu(i,:,:) = reshape(n2011_HF_032mu.sigma2(i,:)+n2011_bedopen_032mu.sigma2(i,:)+n2011_bedslip_032mu.sigma2(i,:),81,81)+princ_sigma2_winter_vq;
n2011_allfour.sigma1_39mu(i,:,:) = reshape(n2011_HF_39mu.sigma1(i,:)+n2011_bedopen_39mu.sigma1(i,:)+n2011_bedslip_39mu.sigma1(i,:),81,81)+princ_sigma1_winter_vq;
n2011_allfour.sigma2_39mu(i,:,:) = reshape(n2011_HF_39mu.sigma2(i,:)+n2011_bedopen_39mu.sigma2(i,:)+n2011_bedslip_39mu.sigma2(i,:),81,81)+princ_sigma2_winter_vq;

% NL n2011_allfour. 
n2011_allfour.princ_sigma1_NL(i,1) = nanmean(nanmean(n2011_HF.sigma1(i,indNL_surf)+...
    n2011_bedopen.sigma1(i,indNL_surf)+n2011_bedslip.sigma1(i,indNL_surf)+princ_sigma1_winter_vq_vec(indNL_surf)));
n2011_allfour.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2011_HF_032mu.sigma1(i,indNL_surf)+...
    n2011_bedopen_032mu.sigma1(i,indNL_surf)+n2011_bedslip_032mu.sigma1(i,indNL_surf)+princ_sigma1_winter_vq_vec(indNL_surf)));
n2011_allfour.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2011_HF_39mu.sigma1(i,indNL_surf)+...
    n2011_bedopen_39mu.sigma1(i,indNL_surf)+n2011_bedslip_39mu.sigma1(i,indNL_surf)+princ_sigma1_winter_vq_vec(indNL_surf)));

% NNL n2011_allfour. 
n2011_allfour.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2011_HF.sigma1(i,indNNL_surf)+...
    n2011_bedopen.sigma1(i,indNNL_surf)+n2011_bedslip.sigma1(i,indNNL_surf)+princ_sigma1_winter_vq_vec(indNNL_surf)));
n2011_allfour.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2011_HF_032mu.sigma1(i,indNNL_surf)+...
    n2011_bedopen_032mu.sigma1(i,indNNL_surf)+n2011_bedslip_032mu.sigma1(i,indNNL_surf)+princ_sigma1_winter_vq_vec(indNNL_surf)));
n2011_allfour.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2011_HF_39mu.sigma1(i,indNNL_surf)+...
    n2011_bedopen_39mu.sigma1(i,indNNL_surf)+n2011_bedslip_39mu.sigma1(i,indNNL_surf)+princ_sigma1_winter_vq_vec(indNNL_surf)));

% SNL n2011_allfour. 
n2011_allfour.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2011_HF.sigma1(i,indSNL_surf)+...
    n2011_bedopen.sigma1(i,indSNL_surf)+n2011_bedslip.sigma1(i,indSNL_surf)+princ_sigma1_winter_vq_vec(indSNL_surf)));
n2011_allfour.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2011_HF_032mu.sigma1(i,indSNL_surf)+...
    n2011_bedopen_032mu.sigma1(i,indSNL_surf)+n2011_bedslip_032mu.sigma1(i,indSNL_surf)+princ_sigma1_winter_vq_vec(indSNL_surf)));
n2011_allfour.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2011_HF_39mu.sigma1(i,indSNL_surf)+...
    n2011_bedopen_39mu.sigma1(i,indSNL_surf)+n2011_bedslip_39mu.sigma1(i,indSNL_surf)+princ_sigma1_winter_vq_vec(indSNL_surf)));

% theta all four
n2011_allfour.theta(i,:) = 0.5.*atan2(2.*(Sxy_flow_vec_vq_vec'+n2011_HF.sigma_xy_vert(i,:)+n2011_bedopen.sigma_xy_open(i,:)+n2011_bedslip.sigma_xy_slip(i,:)),...
    ((Sxx_flow_vec_vq_vec'+n2011_HF.sigma_xx_vert(i,:)+n2011_bedopen.sigma_xx_open(i,:)+n2011_bedslip.sigma_xx_slip(i,:))-...
    (Syy_flow_vec_vq_vec'+n2011_HF.sigma_yy_vert(i,:)+n2011_bedopen.sigma_yy_open(i,:)+n2011_bedslip.sigma_yy_slip(i,:))));
n2011_allfour.theta_mat(i,:,:) = reshape(n2011_allfour.theta(i,:),81,81);
n2011_allfour.theta_NL(i,1) = nanmean(nanmean(n2011_allfour.theta(i,indNL_surf)));
n2011_allfour.theta_NNL(i,1) = nanmean(nanmean(n2011_allfour.theta(i,indNNL_surf)));
n2011_allfour.theta_SNL(i,1) = nanmean(nanmean(n2011_allfour.theta(i,indSNL_surf)));

n2011_allfour.theta_032mu(i,:) = 0.5.*atan2(2.*(Sxy_flow_vec_vq_vec'+n2011_HF_032mu.sigma_xy_vert(i,:)+n2011_bedopen_032mu.sigma_xy_open(i,:)+n2011_bedslip_032mu.sigma_xy_slip(i,:)),...
    ((Sxx_flow_vec_vq_vec'+n2011_HF_032mu.sigma_xx_vert(i,:)+n2011_bedopen_032mu.sigma_xx_open(i,:)+n2011_bedslip_032mu.sigma_xx_slip(i,:))-...
    (Syy_flow_vec_vq_vec'+n2011_HF_032mu.sigma_yy_vert(i,:)+n2011_bedopen_032mu.sigma_yy_open(i,:)+n2011_bedslip_032mu.sigma_yy_slip(i,:))));
n2011_allfour.theta_NL_032mu(i,1) = nanmean(nanmean(n2011_allfour.theta_032mu(i,indNL_surf)));
n2011_allfour.theta_NNL_032mu(i,1) = nanmean(nanmean(n2011_allfour.theta_032mu(i,indNNL_surf)));
n2011_allfour.theta_SNL_032mu(i,1) = nanmean(nanmean(n2011_allfour.theta_032mu(i,indSNL_surf)));

n2011_allfour.theta_39mu(i,:) = 0.5.*atan2(2.*(Sxy_flow_vec_vq_vec'+n2011_HF_39mu.sigma_xy_vert(i,:)+n2011_bedopen_39mu.sigma_xy_open(i,:)+n2011_bedslip_39mu.sigma_xy_slip(i,:)),...
    ((Sxx_flow_vec_vq_vec'+n2011_HF_39mu.sigma_xx_vert(i,:)+n2011_bedopen_39mu.sigma_xx_open(i,:)+n2011_bedslip_39mu.sigma_xx_slip(i,:))-...
    (Syy_flow_vec_vq_vec'+n2011_HF_39mu.sigma_yy_vert(i,:)+n2011_bedopen_39mu.sigma_yy_open(i,:)+n2011_bedslip_39mu.sigma_yy_slip(i,:))));
n2011_allfour.theta_NL_39mu(i,1) = nanmean(nanmean(n2011_allfour.theta_39mu(i,indNL_surf)));
n2011_allfour.theta_NNL_39mu(i,1) = nanmean(nanmean(n2011_allfour.theta_39mu(i,indNNL_surf)));
n2011_allfour.theta_SNL_39mu(i,1) = nanmean(nanmean(n2011_allfour.theta_39mu(i,indSNL_surf)));
end
%%
% save n2011_allfour_surfC.mat n2011_allfour

%% Paper Fig8(?) decomposition into hf, bed open, bed slip, winter, 2011 
for i=1280 % 2011 t_{3}

load BWR.mat    
TriangleSize = 5;

% dock at eel pond 
sky_blue = [111, 169, 228]./255; % SNL
metal = [87, 115, 131]./255; % NNL
oar = [251, 219, 154]./255; % NENL
handle = [161, 37, 49]./255; % NL
dark_oar = [164, 114, 63]./255;
 
fig1 = figure('Units','centimeters','Position',0.655.*[2 0 40 50.5]); clf
    % 2011    
    axe1 = axes('Position',[0.03 0.84 0.1768 0.14],'Box','on','NextPlot','add');   
    axe2 = axes('Position',[0.03 0.62 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    axe4 = axes('Position',[0.26 0.62 0.235 0.14],'Box','on','NextPlot','add');
    
    %  column 3
    axe9 = axes('Position',[0.03 0.45 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'TickLabel', 'none', 'XTick', []);    
    axe10 = axes('Position',[0.26 0.45 0.235 0.14],'Box','on','NextPlot','add');
    
    %  column 4
    axe13 = axes('Position',[0.03 0.28 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'TickLabel', 'none', 'XTick', []);    
    axe14 = axes('Position',[0.26 0.28 0.235 0.14],'Box','on','NextPlot','add');
    
    %  column 5
    axe17 = axes('Position',[0.03 0.11 0.1768 0.14],'Box','on','NextPlot','add');    
    axe18 = axes('Position',[0.26 0.11 0.235 0.14],'Box','on','NextPlot','add');

    % 2012        
    axe22 = axes('Position',[0.53 0.62 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
  
    axe24 = axes('Position',[0.76 0.62 0.235 0.14],'Box','on','NextPlot','add');
    
    %  column 3
    axe29 = axes('Position',[0.53 0.45 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'TickLabel', 'none', 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', [])   
    axe210 = axes('Position',[0.76 0.45 0.235 0.14],'Box','on','NextPlot','add');
    
    %  column 4
    axe213 = axes('Position',[0.53 0.28 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'TickLabel', 'none', 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', [])    
    axe214 = axes('Position',[0.76 0.28 0.235 0.14],'Box','on','NextPlot','add');
 
    %  column 5
    axe217 = axes('Position',[0.53 0.11 0.1768 0.14],'Box','on','NextPlot','add');
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', [])    
    axe218 = axes('Position',[0.76 0.11 0.235 0.14],'Box','on','NextPlot','add');

% principle stress 1 -- winter
axes(axe1)
vvKPA=[-2005:10:2005]; v150=200;
[C5,h6]=contourf(xx_fine,yy_fine,princ_sigma1_winter_vq./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,princ_sigma1_winter_vq./1e3,[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)
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
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end  
caxis([-1000 1000]);
colormap(axe1,BWR); 
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-5.5, 5.5, 'a','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',8,'FontName','Avenir'); 
xlabel(' x [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;
title('[\sigma,\theta]_{1,winter}','FontSize',9,'FontName','Avenir')

% principle stress 1 -- HF
axes(axe2)
[C5,h6]=contourf(xx_fine,yy_fine,reshape(n2011_HF.sigma1(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,reshape(n2011_HF.sigma1(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
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
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2011_HF.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2011_HF.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2011_HF.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2011_HF.theta(i,:))+90),81*81,1));
%theta_P tick marks lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end    
caxis([-1000 1000]);
colormap(axe2,BWR); 
text(-5.5, 5.5, 'b','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
%xlabel(' x [ km ]','FontSize',8,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;
title('[\sigma,\theta]_{1,HF}(2011:{\itt}_{3})','FontSize',9,'FontName','Avenir')
text(4,9,'2011 L1A Drainage','FontSize',15,'FontName','Avenir','FontAngle','italic')
text(38,9,'2012 L1A Drainage','FontSize',15,'FontName','Avenir','FontAngle','italic')

% principle stress 1, HF open each lake basin
axes(axe4)
text(168.05+0.32, 2800, 'f','FontSize',9,'FontWeight','bold');
plot(time_2011(600:1995),n2011_HF.princ_sigma1_NNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',metal); hold on;
plot(time_2011(600:1995),n2011_HF.princ_sigma1_NL(600:1995)./1e3,'-','LineWidth',1.5,'Color',handle); hold on;
plot(time_2011(600:1995),n2011_HF.princ_sigma1_SNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',sky_blue); hold on;
plot([168.85 168.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
text(168.85-0.05,3250,'{\itt}_{1}','Rotation',0,'FontSize',8,'FontName','Avenir'); 
text(169.21-0.05,3250,'{\itt}_{2}','Rotation',0,'FontSize',8,'FontName','Avenir'); 
text(169.32-0.05,3250,'{\itt}_{3}','FontSize',8,'FontName','Avenir'); 
text(170.32-0.05,3250,'{\itt}_{4}','Rotation',0,'FontSize',8,'FontName','Avenir'); 
plot([169.32+(6/24) 169.32+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
text(169.32-0.07+(4/24),3250,'{\it\tau} = {6 hr}','Rotation',0,'FontSize',8,'FontName','Avenir','Color',[0.2 0.55 0.2]); 
% plus minus 032GPs and 39GPa
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_HF.princ_sigma1_NNL_39mu(600:1995)./1e3;flipud(n2011_HF.princ_sigma1_NNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_HF.princ_sigma1_SNL_39mu(600:1995)./1e3;flipud(n2011_HF.princ_sigma1_SNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_HF.princ_sigma1_NL_39mu(600:1995)./1e3;flipud(n2011_HF.princ_sigma1_NL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);
plot(time_2011(600:1995),n2011_HF.princ_sigma1_NNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',metal); 
plot(time_2011(600:1995),n2011_HF.princ_sigma1_SNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',sky_blue); 
plot(time_2011(600:1995),n2011_HF.princ_sigma1_NL(600:1995)./1e3,'-','LineWidth',1.5,'Color',handle); 
plot(time_2011(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);
legend('L1C','L1A','L1B','Location','NorthEast'); legend boxoff
ylim([-500 3000]);  xlim([170.32-2 170.32]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[168.5:0.5:170],'xticklabel',[]); 
grid on
ylabel('\sigma_{1,HF} [ kPa ]','FontSize',8,'FontName','Avenir');

% principle stress 1 -- bed open
axes(axe9)
[C5,h6]=contourf(xx_fine,yy_fine,reshape(n2011_bedopen.sigma1(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,reshape(n2011_bedopen.sigma1(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
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
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2011_bedopen.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2011_bedopen.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2011_bedopen.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2011_bedopen.theta(i,:))+90),81*81,1));
%theta_P tick marks lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end   
caxis([-1000 1000]);
colormap(axe9,BWR); 
text(-5.5, 5.5, 'c','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;
title('[\sigma,\theta]_{1,bedopen}(2011:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed open, each lake basin
axes(axe10)
text(168.05+0.32, 2800, 'g','FontSize',9,'FontWeight','bold');
plot(time_2011(600:1995),n2011_bedopen.princ_sigma1_NNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',metal); hold on;
plot(time_2011(600:1995),n2011_bedopen.princ_sigma1_NL(600:1995)./1e3,'-','LineWidth',1.5,'Color',handle); hold on;
plot(time_2011(600:1995),n2011_bedopen.princ_sigma1_SNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',sky_blue); hold on;
plot([168.85 168.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
text(169.7, 800, '1.5 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-5);
text(169.7, 250, '0.32 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',0);
text(169.7, 1650, '3.9 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-10);
plot([169.32+(6/24) 169.32+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
% plus minus 032GPs and 39GPa
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_bedopen.princ_sigma1_NNL_39mu(600:1995)./1e3;flipud(n2011_bedopen.princ_sigma1_NNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_bedopen.princ_sigma1_SNL_39mu(600:1995)./1e3;flipud(n2011_bedopen.princ_sigma1_SNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_bedopen.princ_sigma1_NL_39mu(600:1995)./1e3;flipud(n2011_bedopen.princ_sigma1_NL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);
plot(time_2011(600:1995),n2011_bedopen.princ_sigma1_NNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',metal); 
plot(time_2011(600:1995),n2011_bedopen.princ_sigma1_SNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',sky_blue); 
plot(time_2011(600:1995),n2011_bedopen.princ_sigma1_NL(600:1995)./1e3,'-','LineWidth',1.5,'Color',handle); 
plot(time_2011(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);
ylim([-500 3000]);  xlim([170.32-2 170.32]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[168.5:0.5:170],'xticklabel',[]); 
grid on
ylabel('\sigma_{1,bedopen} [ kPa ]','FontSize',8,'FontName','Avenir');


% principle stress 1 -- bed slip 
axes(axe13)
[C5,h6]=contourf(xx_fine,yy_fine,reshape(n2011_bedslip.sigma1(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,reshape(n2011_bedslip.sigma1(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
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
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2011_bedslip.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2011_bedslip.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2011_bedslip.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2011_bedslip.theta(i,:))+90),81*81,1));
%theta_P tick lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end  
caxis([-1000 1000]);
colormap(axe13,BWR); 
text(-5.5, 5.5, 'd','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;
title('[\sigma,\theta]_{1,bedslip}(2011:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed slip, each lake basin
axes(axe14)
text(168.05+0.32, 2800, 'h','FontSize',9,'FontWeight','bold');
plot(time_2011(600:1995),n2011_bedslip.princ_sigma1_NNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',metal); hold on;
plot(time_2011(600:1995),n2011_bedslip.princ_sigma1_NL(600:1995)./1e3,'-','LineWidth',1.5,'Color',handle); hold on;
plot(time_2011(600:1995),n2011_bedslip.princ_sigma1_SNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',sky_blue); hold on;
plot([168.85 168.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32+(6/24) 169.32+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
% plus minus 032GPs and 39GPa
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_bedslip.princ_sigma1_NNL_39mu(600:1995)./1e3;flipud(n2011_bedslip.princ_sigma1_NNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_bedslip.princ_sigma1_SNL_39mu(600:1995)./1e3;flipud(n2011_bedslip.princ_sigma1_SNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_bedslip.princ_sigma1_NL_39mu(600:1995)./1e3;flipud(n2011_bedslip.princ_sigma1_NL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);
plot(time_2011(600:1995),n2011_bedslip.princ_sigma1_NNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',metal); 
plot(time_2011(600:1995),n2011_bedslip.princ_sigma1_SNL(600:1995)./1e3,'-','LineWidth',1.5,'Color',sky_blue); 
plot(time_2011(600:1995),n2011_bedslip.princ_sigma1_NL(600:1995)./1e3,'-','LineWidth',1.5,'Color',handle); 
plot(time_2011(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);
legend('L1C','L1A','L1B','Location','NorthEast'); legend boxoff
ylim([-500 3000]);  xlim([170.32-2 170.32]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[168.5:0.5:170],'xticklabel',[]); 
grid on
ylabel('\sigma_{1,bedslip} [ kPa ]','FontSize',8,'FontName','Avenir');

% principle stress 1 -- all 3 sources + winter = ALLFOUR
axes(axe17)
[C5,h6]=contourf(xx_fine,yy_fine,squeeze(n2011_allfour.sigma1(i,:,:)./1e3),vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,squeeze(n2011_allfour.sigma1(i,:,:)./1e3),[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
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
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end 
caxis([-1000 1000]);
colormap(axe17,BWR); 
text(-5.5, 5.5, 'e','FontSize',9,'FontWeight','bold');
xlabel(' x [ km ]','FontSize',8,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',8,'FontName','Avenir'); 
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;
t3=colorbar('EastOutside');
set(t3,'YTick',[-1000,-500,0,500,1000],'TickDirection','out','FontSize',9); 
hold all
set(get(t3,'xlabel'),'String','\sigma_{1}  [ kPa ]','FontSize',10);
set(t3, 'Position', [0.22 0.84 0.005 0.14]);
title('[\sigma,\theta]_{1}(2011:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed slip, each lake basin
axes(axe18)
text(168.05+0.32, 2800, 'i','FontSize',9,'FontWeight','bold');
plot(time_2011(600:1995),(n2011_allfour.princ_sigma1_NNL(600:1995))./1e3,'-','LineWidth',1.5,'Color',metal); 
plot(time_2011(600:1995),(n2011_allfour.princ_sigma1_SNL(600:1995))./1e3,'-','LineWidth',1.5,'Color',sky_blue); 
plot(time_2011(600:1995),(n2011_allfour.princ_sigma1_NL(600:1995))./1e3,'-','LineWidth',1.5,'Color',handle); 
plot([168.85 168.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
text(169.7, 950, '1.5 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-5);
text(169.7, 350, '0.32 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',0);
text(169.7, 2100, '3.9 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-10);
plot([169.32+(6/24) 169.32+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
% plus minus 032GPs and 39GPa
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_allfour.princ_sigma1_NNL_39mu(600:1995)./1e3;flipud(n2011_allfour.princ_sigma1_NNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_allfour.princ_sigma1_SNL_39mu(600:1995)./1e3;flipud(n2011_allfour.princ_sigma1_SNL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2011(600:1995);flipud(time_2011(600:1995))],...
    [n2011_allfour.princ_sigma1_NL_39mu(600:1995)./1e3;flipud(n2011_allfour.princ_sigma1_NL_032mu(600:1995)./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);
plot(time_2011(600:1995),(n2011_allfour.princ_sigma1_NNL(600:1995))./1e3,'-','LineWidth',1.5,'Color',metal); 
plot(time_2011(600:1995),(n2011_allfour.princ_sigma1_SNL(600:1995))./1e3,'-','LineWidth',1.5,'Color',sky_blue); 
plot(time_2011(600:1995),(n2011_allfour.princ_sigma1_NL(600:1995))./1e3,'-','LineWidth',1.5,'Color',handle); 
plot(time_2011(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);
ylim([-500 3000]);  xlim([170.32-2 170.32]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[168.5:0.5:170]); 
grid on
ylabel('\sigma_{1} [ kPa ]','FontSize',8,'FontName','Avenir');
xlabel('Day of Year, 2011','FontSize',8,'FontName','Avenir');
end

%% START OF 2012 PLOTTING
% time and spatial info
load time_2012.mat

%% hydro-fracture: load stresses, calculate mean within lake basins (2012)
load sigma_ij_2012_SLder_HF_surfC_032mu.mat
n2012_HF_032mu = sigma_ij_2012_SLder_HF;
load sigma_ij_2012_SLder_HF_surfC_39mu.mat
n2012_HF_39mu = sigma_ij_2012_SLder_HF;
load sigma_ij_2012_SLder_HF_surfC.mat
n2012_HF = sigma_ij_2012_SLder_HF;

for i=1:length(time_2012)
% NL n2012_HF. 
n2012_HF.princ_sigma1_NL(i,1) = nanmean(nanmean(n2012_HF.sigma1(i,indNL_surf)));
n2012_HF.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2012_HF_032mu.sigma1(i,indNL_surf)));
n2012_HF.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2012_HF_39mu.sigma1(i,indNL_surf)));

% NNL n2012_HF. 
n2012_HF.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2012_HF.sigma1(i,indNNL_surf)));
n2012_HF.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2012_HF_032mu.sigma1(i,indNNL_surf)));
n2012_HF.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2012_HF_39mu.sigma1(i,indNNL_surf)));

% SNL n2012_HF.  
n2012_HF.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2012_HF.sigma1(i,indSNL_surf)));
n2012_HF.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2012_HF_032mu.sigma1(i,indSNL_surf)));
n2012_HF.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2012_HF_39mu.sigma1(i,indSNL_surf)));
end

%% basal opening: load stresses, calculate mean within lake basins (2012)
load sigma_ij_2012_SLder_bedopen_surfC_032mu.mat
n2012_bedopen_032mu = sigma_ij_2012_SLder_bedopen;
load sigma_ij_2012_SLder_bedopen_surfC_39mu.mat
n2012_bedopen_39mu = sigma_ij_2012_SLder_bedopen;
load sigma_ij_2012_SLder_bedopen_surfC.mat
n2012_bedopen = sigma_ij_2012_SLder_bedopen;

for i=1:length(time_2012)
% NL n2012_bedopen. 
n2012_bedopen.princ_sigma1_NL(i,1) = nanmean(nanmean(n2012_bedopen.sigma1(i,indNL_surf)));
n2012_bedopen.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2012_bedopen_032mu.sigma1(i,indNL_surf)));
n2012_bedopen.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2012_bedopen_39mu.sigma1(i,indNL_surf)));

% NNL n2012_bedopen. 
n2012_bedopen.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2012_bedopen.sigma1(i,indNNL_surf)));
n2012_bedopen.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2012_bedopen_032mu.sigma1(i,indNNL_surf)));
n2012_bedopen.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2012_bedopen_39mu.sigma1(i,indNNL_surf)));

% SNL n2012_bedopen.  
n2012_bedopen.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2012_bedopen.sigma1(i,indSNL_surf)));
n2012_bedopen.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2012_bedopen_032mu.sigma1(i,indSNL_surf)));
n2012_bedopen.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2012_bedopen_39mu.sigma1(i,indSNL_surf)));
end

%% basal slip: load stresses, calculate mean within lake basins
load sigma_ij_2012_SLder_bedslip_surfC_032mu.mat
n2012_bedslip_032mu = sigma_ij_2012_SLder_bedslip;
load sigma_ij_2012_SLder_bedslip_surfC_39mu.mat
n2012_bedslip_39mu = sigma_ij_2012_SLder_bedslip;
load sigma_ij_2012_SLder_bedslip_surfC.mat
n2012_bedslip = sigma_ij_2012_SLder_bedslip;

for i=1:length(time_2012)
% NL n2012_bedslip. 
n2012_bedslip.princ_sigma1_NL(i,1) = nanmean(nanmean(n2012_bedslip.sigma1(i,indNL_surf)));
n2012_bedslip.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2012_bedslip_032mu.sigma1(i,indNL_surf)));
n2012_bedslip.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2012_bedslip_39mu.sigma1(i,indNL_surf)));

% NNL n2012_bedslip. 
n2012_bedslip.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2012_bedslip.sigma1(i,indNNL_surf)));
n2012_bedslip.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2012_bedslip_032mu.sigma1(i,indNNL_surf)));
n2012_bedslip.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2012_bedslip_39mu.sigma1(i,indNNL_surf)));

% SNL n2012_bedslip.  
n2012_bedslip.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2012_bedslip.sigma1(i,indSNL_surf)));
n2012_bedslip.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2012_bedslip_032mu.sigma1(i,indSNL_surf)));
n2012_bedslip.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2012_bedslip_39mu.sigma1(i,indSNL_surf)));
end

%% winter+HF+bed_open+basal_slip: calculate sigma1 + theta for all four contributions

for i=1:length(time_2012)
    
n2012_allfour.sigma1(i,:,:) = reshape(n2012_HF.sigma1(i,:)+n2012_bedopen.sigma1(i,:)+n2012_bedslip.sigma1(i,:),81,81)+princ_sigma1_winter_vq;
n2012_allfour.sigma2(i,:,:) = reshape(n2012_HF.sigma2(i,:)+n2012_bedopen.sigma2(i,:)+n2012_bedslip.sigma2(i,:),81,81)+princ_sigma2_winter_vq;
n2012_allfour.sigma1_032mu(i,:,:) = reshape(n2012_HF_032mu.sigma1(i,:)+n2012_bedopen_032mu.sigma1(i,:)+n2012_bedslip_032mu.sigma1(i,:),81,81)+princ_sigma1_winter_vq;
n2012_allfour.sigma2_032mu(i,:,:) = reshape(n2012_HF_032mu.sigma2(i,:)+n2012_bedopen_032mu.sigma2(i,:)+n2012_bedslip_032mu.sigma2(i,:),81,81)+princ_sigma2_winter_vq;
n2012_allfour.sigma1_39mu(i,:,:) = reshape(n2012_HF_39mu.sigma1(i,:)+n2012_bedopen_39mu.sigma1(i,:)+n2012_bedslip_39mu.sigma1(i,:),81,81)+princ_sigma1_winter_vq;
n2012_allfour.sigma2_39mu(i,:,:) = reshape(n2012_HF_39mu.sigma2(i,:)+n2012_bedopen_39mu.sigma2(i,:)+n2012_bedslip_39mu.sigma2(i,:),81,81)+princ_sigma2_winter_vq;

% NL n2012_allfour. 
n2012_allfour.princ_sigma1_NL(i,1) = nanmean(nanmean(n2012_HF.sigma1(i,indNL_surf)+...
    n2012_bedopen.sigma1(i,indNL_surf)+n2012_bedslip.sigma1(i,indNL_surf)+princ_sigma1_winter_vq_vec(indNL_surf)));
n2012_allfour.princ_sigma1_NL_032mu(i,1) = nanmean(nanmean(n2012_HF_032mu.sigma1(i,indNL_surf)+...
    n2012_bedopen_032mu.sigma1(i,indNL_surf)+n2012_bedslip_032mu.sigma1(i,indNL_surf)+princ_sigma1_winter_vq_vec(indNL_surf)));
n2012_allfour.princ_sigma1_NL_39mu(i,1) = nanmean(nanmean(n2012_HF_39mu.sigma1(i,indNL_surf)+...
    n2012_bedopen_39mu.sigma1(i,indNL_surf)+n2012_bedslip_39mu.sigma1(i,indNL_surf)+princ_sigma1_winter_vq_vec(indNL_surf)));

% NNL n2012_allfour. 
n2012_allfour.princ_sigma1_NNL(i,1) = nanmean(nanmean(n2012_HF.sigma1(i,indNNL_surf)+...
    n2012_bedopen.sigma1(i,indNNL_surf)+n2012_bedslip.sigma1(i,indNNL_surf)+princ_sigma1_winter_vq_vec(indNNL_surf)));
n2012_allfour.princ_sigma1_NNL_032mu(i,1) = nanmean(nanmean(n2012_HF_032mu.sigma1(i,indNNL_surf)+...
    n2012_bedopen_032mu.sigma1(i,indNNL_surf)+n2012_bedslip_032mu.sigma1(i,indNNL_surf)+princ_sigma1_winter_vq_vec(indNNL_surf)));
n2012_allfour.princ_sigma1_NNL_39mu(i,1) = nanmean(nanmean(n2012_HF_39mu.sigma1(i,indNNL_surf)+...
    n2012_bedopen_39mu.sigma1(i,indNNL_surf)+n2012_bedslip_39mu.sigma1(i,indNNL_surf)+princ_sigma1_winter_vq_vec(indNNL_surf)));

% SNL n2012_allfour. 
n2012_allfour.princ_sigma1_SNL(i,1) = nanmean(nanmean(n2012_HF.sigma1(i,indSNL_surf)+...
    n2012_bedopen.sigma1(i,indSNL_surf)+n2012_bedslip.sigma1(i,indSNL_surf)+princ_sigma1_winter_vq_vec(indSNL_surf)));
n2012_allfour.princ_sigma1_SNL_032mu(i,1) = nanmean(nanmean(n2012_HF_032mu.sigma1(i,indSNL_surf)+...
    n2012_bedopen_032mu.sigma1(i,indSNL_surf)+n2012_bedslip_032mu.sigma1(i,indSNL_surf)+princ_sigma1_winter_vq_vec(indSNL_surf)));
n2012_allfour.princ_sigma1_SNL_39mu(i,1) = nanmean(nanmean(n2012_HF_39mu.sigma1(i,indSNL_surf)+...
    n2012_bedopen_39mu.sigma1(i,indSNL_surf)+n2012_bedslip_39mu.sigma1(i,indSNL_surf)+princ_sigma1_winter_vq_vec(indSNL_surf)));

% theta all four
n2012_allfour.theta(i,:) = 0.5.*atan2(2.*(Sxy_flow_vec_vq_vec'+n2012_HF.sigma_xy_vert(i,:)+n2012_bedopen.sigma_xy_open(i,:)+n2012_bedslip.sigma_xy_slip(i,:)),...
    ((Sxx_flow_vec_vq_vec'+n2012_HF.sigma_xx_vert(i,:)+n2012_bedopen.sigma_xx_open(i,:)+n2012_bedslip.sigma_xx_slip(i,:))-...
    (Syy_flow_vec_vq_vec'+n2012_HF.sigma_yy_vert(i,:)+n2012_bedopen.sigma_yy_open(i,:)+n2012_bedslip.sigma_yy_slip(i,:))));
n2012_allfour.theta_mat(i,:,:) = reshape(n2012_allfour.theta(i,:),81,81);
n2012_allfour.theta_NL(i,1) = nanmean(nanmean(n2012_allfour.theta(i,indNL_surf)));
n2012_allfour.theta_NNL(i,1) = nanmean(nanmean(n2012_allfour.theta(i,indNNL_surf)));
n2012_allfour.theta_SNL(i,1) = nanmean(nanmean(n2012_allfour.theta(i,indSNL_surf)));

n2012_allfour.theta_032mu(i,:) = 0.5.*atan2(2.*(Sxy_flow_vec_vq_vec'+n2012_HF_032mu.sigma_xy_vert(i,:)+n2012_bedopen_032mu.sigma_xy_open(i,:)+n2012_bedslip_032mu.sigma_xy_slip(i,:)),...
    ((Sxx_flow_vec_vq_vec'+n2012_HF_032mu.sigma_xx_vert(i,:)+n2012_bedopen_032mu.sigma_xx_open(i,:)+n2012_bedslip_032mu.sigma_xx_slip(i,:))-...
    (Syy_flow_vec_vq_vec'+n2012_HF_032mu.sigma_yy_vert(i,:)+n2012_bedopen_032mu.sigma_yy_open(i,:)+n2012_bedslip_032mu.sigma_yy_slip(i,:))));
n2012_allfour.theta_NL_032mu(i,1) = nanmean(nanmean(n2012_allfour.theta_032mu(i,indNL_surf)));
n2012_allfour.theta_NNL_032mu(i,1) = nanmean(nanmean(n2012_allfour.theta_032mu(i,indNNL_surf)));
n2012_allfour.theta_SNL_032mu(i,1) = nanmean(nanmean(n2012_allfour.theta_032mu(i,indSNL_surf)));

n2012_allfour.theta_39mu(i,:) = 0.5.*atan2(2.*(Sxy_flow_vec_vq_vec'+n2012_HF_39mu.sigma_xy_vert(i,:)+n2012_bedopen_39mu.sigma_xy_open(i,:)+n2012_bedslip_39mu.sigma_xy_slip(i,:)),...
    ((Sxx_flow_vec_vq_vec'+n2012_HF_39mu.sigma_xx_vert(i,:)+n2012_bedopen_39mu.sigma_xx_open(i,:)+n2012_bedslip_39mu.sigma_xx_slip(i,:))-...
    (Syy_flow_vec_vq_vec'+n2012_HF_39mu.sigma_yy_vert(i,:)+n2012_bedopen_39mu.sigma_yy_open(i,:)+n2012_bedslip_39mu.sigma_yy_slip(i,:))));
n2012_allfour.theta_NL_39mu(i,1) = nanmean(nanmean(n2012_allfour.theta_39mu(i,indNL_surf)));
n2012_allfour.theta_NNL_39mu(i,1) = nanmean(nanmean(n2012_allfour.theta_39mu(i,indNNL_surf)));
n2012_allfour.theta_SNL_39mu(i,1) = nanmean(nanmean(n2012_allfour.theta_39mu(i,indSNL_surf)));
end

%%
save n2012_allfour_surfC.mat n2012_allfour

%% PLOT
for i = 544 % 2012 t_{3} index
%% principle stress 1 -- hydro-fracture (2012)
axes(axe22)
[~,h6]=contourf(xx_fine,yy_fine,reshape(n2012_HF.sigma1(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,reshape(n2012_HF.sigma1(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);

% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% lakes
plot(xy_sta_12(:,1),xy_sta_12(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)

% tick marks combine okada85 and winter
        hypot = 0.2;
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2012_HF.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2012_HF.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2012_HF.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2012_HF.theta(i,:))+90),81*81,1));
%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end
    
caxis([-1000 1000]);
colormap(axe22,BWR); 
text(-5.5, 5.5, 'j','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;

title('[\sigma,\theta]_{1,HF}(2012:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed slip, each lake basin
axes(axe24)
text(160.55+0.35, 2800, 'n','FontSize',9,'FontWeight','bold');
plot(time_2012,n2012_HF.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); hold on;
plot(time_2012,n2012_HF.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); hold on;
plot(time_2012,n2012_HF.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); hold on;

plot([161.20 161.20],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([162.85 162.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85+(6/24) 161.85+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell

text(161.20-0.05,3250,'{\itt}_{1}','Rotation',0,'FontSize',8,'FontName','Avenir'); 
text(161.72-0.05,3250,'{\itt}_{2}','Rotation',0,'FontSize',8,'FontName','Avenir'); 
text(161.85-0.05,3250,'{\itt}_{3}','FontSize',8,'FontName','Avenir'); 
text(162.85-0.05,3250,'{\itt}_{4}','Rotation',0,'FontSize',8,'FontName','Avenir'); 
text(161.85-0.07+(4/24),3250,'{\it\tau} = {6 hr}','Rotation',0,'FontSize',8,'FontName','Avenir','Color',[0.2 0.55 0.2]); 

% plus minus 032GPa and 39GPa
fill([time_2012;flipud(time_2012)],...
    [n2012_HF.princ_sigma1_NNL_39mu./1e3;flipud(n2012_HF.princ_sigma1_NNL_032mu./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_HF.princ_sigma1_SNL_39mu./1e3;flipud(n2012_HF.princ_sigma1_SNL_032mu./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_HF.princ_sigma1_NL_39mu./1e3;flipud(n2012_HF.princ_sigma1_NL_032mu./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);

plot(time_2012,n2012_HF.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); 
plot(time_2012,n2012_HF.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); 
plot(time_2012,n2012_HF.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); 

plot(time_2012(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);

legend('L1C','L1A','L1B','Location','NorthEast'); legend boxoff
ylim([-500 3000]);  xlim([162.85-2 162.85]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[161:0.5:162.5]); 
grid on
ylabel('\sigma_{1,HF} [ kPa ]','FontSize',8,'FontName','Avenir');


%% principle stress 1 -- bed open
axes(axe29)
[~,h6]=contourf(xx_fine,yy_fine,reshape(n2012_bedopen.sigma1(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,reshape(n2012_bedopen.sigma1(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);

% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% lakes
plot(xy_sta_12(:,1),xy_sta_12(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)

% tick marks combine okada85 and winter
        hypot = 0.2;
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2012_bedopen.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2012_bedopen.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2012_bedopen.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2012_bedopen.theta(i,:))+90),81*81,1));

%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end

caxis([-1000 1000]);
colormap(axe29,BWR); 

text(-5.5, 5.5, 'k','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
%xlabel(' x [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;

title('[\sigma,\theta]_{1,bedopen}(2012:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed open, each lake basin
axes(axe210)
text(160.55+0.35, 2800, 'o','FontSize',9,'FontWeight','bold'); hold on
plot(time_2012,n2012_bedopen.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); hold on;
plot(time_2012,n2012_bedopen.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); hold on;
plot(time_2012,n2012_bedopen.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); hold on;

plot([161.20 161.20],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([162.85 162.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85+(6/24) 161.85+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell

% plus minus 032GPs and 39GPa
fill([time_2012;flipud(time_2012)],...
    [n2012_bedopen.princ_sigma1_NNL_39mu./1e3;flipud(n2012_bedopen.princ_sigma1_NNL_032mu./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_bedopen.princ_sigma1_SNL_39mu./1e3;flipud(n2012_bedopen.princ_sigma1_SNL_032mu./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_bedopen.princ_sigma1_NL_39mu./1e3;flipud(n2012_bedopen.princ_sigma1_NL_032mu./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);

plot(time_2012,n2012_bedopen.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); 
plot(time_2012,n2012_bedopen.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); 
plot(time_2012,n2012_bedopen.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); 

plot(time_2012(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);

ylim([-500 3000]);  xlim([162.85-2 162.85]);
legend('L1C','L1A','L1B','Location','NorthEast'); legend boxoff
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[161:0.5:162.5],'xticklabel',[]); 
grid on
ylabel('\sigma_{1,bedopen} [ kPa ]','FontSize',8,'FontName','Avenir');


%% principle stress 1 -- bed slip 
axes(axe213)
[~,h6]=contourf(xx_fine,yy_fine,reshape(n2012_bedslip.sigma1(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,reshape(n2012_bedslip.sigma1(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);

% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% lakes
plot(xy_sta_12(:,1),xy_sta_12(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)

% tick marks combine okada85 and winter
        hypot = 0.2;
% just lake P2
        theta_15GPa_P2_plotX1_w(1,:) = xy_surf_fine(:,1) + hypot.*(reshape(cosd(rad2deg(n2012_bedslip.theta(i,:))+90),81*81,1)); 
        theta_15GPa_P2_plotX2_w(1,:) = xy_surf_fine(:,1) - hypot.*(reshape(cosd(rad2deg(n2012_bedslip.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY1_w(1,:) = xy_surf_fine(:,2) + hypot.*(reshape(sind(rad2deg(n2012_bedslip.theta(i,:))+90),81*81,1));
        theta_15GPa_P2_plotY2_w(1,:) = xy_surf_fine(:,2) - hypot.*(reshape(sind(rad2deg(n2012_bedslip.theta(i,:))+90),81*81,1));

%theta_P tick marks winter + lake
    for j=2200:1:4400
    plot([theta_15GPa_P2_plotX2_w(1,j) theta_15GPa_P2_plotX1_w(1,j)],...
        [theta_15GPa_P2_plotY2_w(1,j) theta_15GPa_P2_plotY1_w(1,j)],...
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end
    
caxis([-1000 1000]);
colormap(axe213,BWR); 

text(-5.5, 5.5, 'l','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
grid on;

title('[\sigma,\theta]_{1,bedslip}(2012:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed slip, each lake basin
axes(axe214)
text(160.55+0.35, 2800, 'p','FontSize',9,'FontWeight','bold');
plot(time_2012,n2012_bedslip.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); hold on;
plot(time_2012,n2012_bedslip.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); hold on;
plot(time_2012,n2012_bedslip.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); hold on;

plot([161.20 161.20],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([162.85 162.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85+(6/24) 161.85+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell

% plus minus 032GPs and 39GPa
fill([time_2012;flipud(time_2012)],...
    [n2012_bedslip.princ_sigma1_NNL_39mu./1e3;flipud(n2012_bedslip.princ_sigma1_NNL_032mu./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_bedslip.princ_sigma1_SNL_39mu./1e3;flipud(n2012_bedslip.princ_sigma1_SNL_032mu./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_bedslip.princ_sigma1_NL_39mu./1e3;flipud(n2012_bedslip.princ_sigma1_NL_032mu./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);

plot(time_2012,n2012_bedslip.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); 
plot(time_2012,n2012_bedslip.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); 
plot(time_2012,n2012_bedslip.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); 

plot(time_2012(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);

legend('L1C','L1A','L1B','Location','NorthEast'); legend boxoff
ylim([-500 3000]);  xlim([162.85-2 162.85]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[161:0.5:162.5],'xticklabel',[]); 
grid on
ylabel('\sigma_{1,bedslip} [ kPa ]','FontSize',8,'FontName','Avenir');


%% principle stress 1 -- all 3 sources + winter = ALL FOUR
axes(axe217)
[C5,h6]=contourf(xx_fine,yy_fine,squeeze(n2012_allfour.sigma1(i,:,:)./1e3),vvKPA);
set(h6,'LineColor','none'); 
contour(xx_fine,yy_fine,squeeze(n2012_allfour.sigma1(i,:,:)./1e3),[v150 v150],'k','LineWidth',1.1);
% boxes
rectangle('position',[-0.55 3 1.4 1.25],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1C
rectangle('position',[-0.25 -0.2 2 1.0],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1A
rectangle('position',[-1.10 -2.4 1.35 0.9],'EdgeColor',[0.2 0.75 0.2],'FaceColor','none','LineWidth',1.5); % L1B
% lakes
plot(xy_sta_12(:,1),xy_sta_12(:,2),'k^','MarkerSize',TriangleSize-1,'MarkerFaceColor','none'); hold on;
scatter(xy_lake(:,1),xy_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bo','filled','MarkerEdgeColor','none')
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',2.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'y-','LineWidth',2.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'y-','LineWidth',2.5)

% tick marks combine okada85 and winter
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
        '-','LineWidth',1.0,'Color',[0.7 0.7 0.7]);   
    end    
    
caxis([-1000 1000]);
colormap(axe217,BWR); 

text(-5.5, 5.5, 'm','FontSize',9,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlabel(' x [ km ]','FontSize',8,'FontName','Avenir');
xlim([-6 6]); ylim([-6 6])
set(gca,'xtick',[-6:2:6],'ytick',[-6:2:6],'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','Layer','Top'); 
set(gca, 'YAxisLocation', 'right')
grid on;

title('[\sigma,\theta]_{1}(2012:{\itt}_{3})','FontSize',9,'FontName','Avenir')

% principle stress 1, bed slip, each lake basin
axes(axe218)
text(160.55+0.35, 2800, 'q','FontSize',9,'FontWeight','bold');
plot(time_2012,n2012_allfour.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); hold on;
plot(time_2012,n2012_allfour.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); hold on;
plot(time_2012,n2012_allfour.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); hold on;
plot([161.20 161.20],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([162.85 162.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85+(6/24) 161.85+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
text(162.2, 800, '1.5 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-10);
text(162.2, 300, '0.32 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-1);
text(162.2, 1750, '3.9 GPa','FontSize',9,'FontName','Avenir','Color',handle,'Rotation',-20);
% plus minus 032GPs and 39GPa
fill([time_2012;flipud(time_2012)],...
    [n2012_allfour.princ_sigma1_NNL_39mu./1e3;flipud(n2012_allfour.princ_sigma1_NNL_032mu./1e3)],'w',...
    'EdgeColor',metal,'FaceColor',metal,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_allfour.princ_sigma1_SNL_39mu./1e3;flipud(n2012_allfour.princ_sigma1_SNL_032mu./1e3)],'w',...
    'EdgeColor',sky_blue,'FaceColor',sky_blue,'FaceAlpha',0.05,'LineWidth',0.9);
fill([time_2012;flipud(time_2012)],...
    [n2012_allfour.princ_sigma1_NL_39mu./1e3;flipud(n2012_allfour.princ_sigma1_NL_032mu./1e3)],'w',...
    'EdgeColor',handle,'FaceColor',handle,'FaceAlpha',0.05,'LineWidth',0.9);
plot(time_2012,n2012_allfour.princ_sigma1_NNL./1e3,'-','LineWidth',1.4,'Color',metal); hold on;
plot(time_2012,n2012_allfour.princ_sigma1_NL./1e3,'-','LineWidth',1.4,'Color',handle); hold on;
plot(time_2012,n2012_allfour.princ_sigma1_SNL./1e3,'-','LineWidth',1.4,'Color',sky_blue); hold on;
plot(time_2012(i),3000,'v','MarkerSize',4,'MarkerFaceColor','k','Color','k','LineWidth',0.1);
ylim([-500 3000]);  xlim([162.85-2 162.85]);
set(gca,'tickdir','in','LineWidth',0.50,'FontSize',8,'FontName','Avenir','xtick',[161:0.5:162.5]);
grid on
ylabel('\sigma_{1} [ kPa ]','FontSize',8,'FontName','Avenir');
xlabel('Day of Year, 2012','FontSize',8,'FontName','Avenir'); 

end

%% stress values within L1A box at t_{2}
%n2011_allfour.princ_sigma1_NL_032mu(1215)./1e3
n2011_allfour.princ_sigma1_NL(1215)./1e3
%n2011_allfour.princ_sigma1_NL_39mu(1215)./1e3
 
%n2012_allfour.princ_sigma1_SNL_032mu(519)./1e3
n2012_allfour.princ_sigma1_SNL(519)./1e3
%n2012_allfour.princ_sigma1_SNL_39mu(519)./1e3    

%% print figure
%figurename=sprintf('paperfig7_nif20112012_decomp_atmaxHF_surfC_20230724.png');
%print(gcf,'-dpng','-r600',figurename);
