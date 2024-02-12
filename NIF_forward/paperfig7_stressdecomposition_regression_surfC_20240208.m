%% regress stress decomposition against all four stress contributors, 2011 and 2012
% 28 June 2021 LAS -- initial file
% 19 June 2023 LAS -- Include S. Larochelle stress derivation
% 12 Feb 2924 LAS -- let's have a look at the regressions
clear all; close all;

%% 2011/2012 saved NIF outputs and forward model geometry
load('../NIF_L1A/Gsurface500_strain_surfacecrack.mat') % Green functions for 20 km X 20 km surface using
% 0.5-km spacing across surface and the NL set-up
load time_2011.mat; % time vector 
time_2011 = time_2011(600:1995);
load time_2012.mat; % time vector

%% dock at eel pond 
sky_blue = [111, 169, 228]./255; % SNL, winter
metal = [87, 115, 131]./255; % NNL, bedslip
oar = [251, 219, 154]./255; % NENL, HF
handle = [161, 37, 49]./255; % NL, bed open
dark_oar = [164, 114, 63]./255;

%% calculate xy grid for comparison
% Gsurface500 (green's function matrix) --> displacement
% xy for the 20 km x 20 km surface
[xx_fine,yy_fine] = meshgrid(-20:0.5:20, -20:0.5:20); % [ km ]
xy_surf_fine = horzcat(reshape(xx_fine,81*81,1),reshape(yy_fine,81*81,1));
xy_surf = Gsurface500.xy_surf; %horzcat(reshape(xx,81*81,1), reshape(yy,81*81,1));
Nsurface = Gsurface500.Nsurface; %length(xy_surf);
% xy for the 6-km grid [these points are included in regression]
[xx_fine6,yy_fine6] = meshgrid(-6:0.5:6, -6:0.5:6); %[ km ]
xy_surf_fine6 = horzcat(reshape(xx_fine6,25*25,1),reshape(yy_fine6,25*25,1));
Nsurface6 = length(xy_surf_fine6);
X1 = 29; Y1 = X1; % index boundaries for matrix output
X2 = 53; Y2 = X2; % index boundaries for matrix output
% indices for timeseries in vector format
[ind_vec] = find(xy_surf_fine(:,1)<=6.0 & xy_surf_fine(:,1)>=-6.0 & ...
    xy_surf_fine(:,2)<=6.0 & xy_surf_fine(:,2)>=-6.0);

%% winter velocities: load and interpolate to Nsurface locations
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

%% 2011 stress decompositions (20 km x 20 km surface) to (6 km x 6 km)
% winter
W = reshape(princ_sigma1_winter_vq(X1:X2,Y1:Y2),1,Nsurface6);
% HF
load sigma_ij_2011_SLder_HF_surfC.mat
HF = sigma_ij_2011_SLder_HF.sigma1(600:1995,ind_vec); % size = 1396 X 625
clear sigma_ij_2011_SLder_HF
% bedslip
load sigma_ij_2011_SLder_bedslip_surfC.mat
BS = sigma_ij_2011_SLder_bedslip.sigma1(600:1995,ind_vec);
clear sigma_ij_2011_SLder_bedslip
% bedopen
load sigma_ij_2011_SLder_bedopen_surfC.mat
BO = sigma_ij_2011_SLder_bedopen.sigma1(600:1995,ind_vec);
clear sigma_ij_2011_SLder_bedopen
% total 
load n2011_allfour_surfC.mat
ALL = n2011_allfour.sigma1(600:1995,ind_vec);
clear n2011_allfour

%% Pearson's correlation coefficient r (from -1 to +1) and p-values
for i=1:length(time_2011)
    % W vs. ALL
    [R_w,p_w] = corrcoef(W,ALL(i,:));
    W_ALL_stats(i,1) = R_w(1,2); W_ALL_stats(i,2) = p_w(1,2);
    
    % HF vs. ALL
    [R_HF,p_HF] = corrcoef(HF(i,:),ALL(i,:));
    HF_ALL_stats(i,1) = R_HF(1,2); HF_ALL_stats(i,2) = p_HF(1,2);
    
    % BS vs. ALL
    [R_BS,p_BS] = corrcoef(BS(i,:),ALL(i,:));
    BS_ALL_stats(i,1) = R_BS(1,2); BS_ALL_stats(i,2) = p_BS(1,2);

    % BO vs. ALL
    [R_BO,p_BO] = corrcoef(BO(i,:),ALL(i,:));
    BO_ALL_stats(i,1) = R_BO(1,2); BO_ALL_stats(i,2) = p_BO(1,2);
end

% r values for (p-values < 0.01) turn into NaNs
threshold = 0.01; % 99% significance
I = find(W_ALL_stats(:,2) >= threshold); W_ALL_stats(:,3) = W_ALL_stats(:,1); W_ALL_stats(I,3) = NaN;
I = find(HF_ALL_stats(:,2) >= threshold); HF_ALL_stats(:,3) = HF_ALL_stats(:,1); HF_ALL_stats(I,3) = NaN;
I = find(BS_ALL_stats(:,2) >= threshold); BS_ALL_stats(:,3) = BS_ALL_stats(:,1); BS_ALL_stats(I,3) = NaN;
I = find(BO_ALL_stats(:,2) >= threshold); BO_ALL_stats(:,3) = BO_ALL_stats(:,1); BO_ALL_stats(I,3) = NaN;


%% 2011 spot check regressions
for i=208 % 168.5
    % W vs. ALL
    [R_w,p_w] = corrcoef(W,ALL(i,:));
    W_ALL_stats(i,1) = R_w(1,2); W_ALL_stats(i,2) = p_w(1,2);
    
    % HF vs. ALL
    [R_HF,p_HF] = corrcoef(HF(i,:),ALL(i,:));
    HF_ALL_stats(i,1) = R_HF(1,2); HF_ALL_stats(i,2) = p_HF(1,2);
    
    % BS vs. ALL
    [R_BS,p_BS] = corrcoef(BS(i,:),ALL(i,:));
    BS_ALL_stats(i,1) = R_BS(1,2); BS_ALL_stats(i,2) = p_BS(1,2);

    % BO vs. ALL
    [R_BO,p_BO] = corrcoef(BO(i,:),ALL(i,:));
    BO_ALL_stats(i,1) = R_BO(1,2); BO_ALL_stats(i,2) = p_BO(1,2);

    figure(1) % = figure('Units','centimeters','Position',[2 1 20 20]);
    subplot(221)
    plot(W./1e3,ALL(i,:)./1e3,'.','Color',metal); hold on
    figurename1 = sprintf('winter (time = 2011 168.5) (R = %5.4f)',R_w(1,2));
    title(figurename1)
    xlabel('\sigma_{1,winter} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
    
    subplot(222)
    plot(HF(i,:)./1e3,ALL(i,:)./1e3,'.','Color',[0.2 0.55 0.2]); hold on
    figurename1 = sprintf('HF (time = 2011 168.5) (R = %5.4f)',R_HF(1,2));
    title(figurename1)
    xlabel('\sigma_{1,HF} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')


    subplot(223)
    plot(BO(i,:)./1e3,ALL(i,:)./1e3,'.','Color',handle); hold on
    figurename1 = sprintf('bedopen (time = 2011 168.5) (R = %5.4f)',R_BO(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedopen} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')

    subplot(224)
    plot(BS(i,:)./1e3,ALL(i,:)./1e3,'.','Color',sky_blue); hold on
    figurename1 = sprintf('bedslip (time = 2011 168.5) (R = %5.4f)',R_BS(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedslip} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
end
% figurename1=sprintf('paperfig8_decomp_regression_surfC_check_1685_20240208.png');
% print(gcf,'-dpng','-r600',figurename1);
%% 2011 169.0
for i=496 % 169.0
    % W vs. ALL
    [R_w,p_w] = corrcoef(W,ALL(i,:));
    W_ALL_stats(i,1) = R_w(1,2); W_ALL_stats(i,2) = p_w(1,2);
    
    % HF vs. ALL
    [R_HF,p_HF] = corrcoef(HF(i,:),ALL(i,:));
    HF_ALL_stats(i,1) = R_HF(1,2); HF_ALL_stats(i,2) = p_HF(1,2);
    
    % BS vs. ALL
    [R_BS,p_BS] = corrcoef(BS(i,:),ALL(i,:));
    BS_ALL_stats(i,1) = R_BS(1,2); BS_ALL_stats(i,2) = p_BS(1,2);

    % BO vs. ALL
    [R_BO,p_BO] = corrcoef(BO(i,:),ALL(i,:));
    BO_ALL_stats(i,1) = R_BO(1,2); BO_ALL_stats(i,2) = p_BO(1,2);

    figure(2) % = figure('Units','centimeters','Position',[2 1 20 20]);
    subplot(221)
    plot(W./1e3,ALL(i,:)./1e3,'.','Color',metal); hold on
    figurename1 = sprintf('winter (time = 2011 169.0) (R = %5.4f)',R_w(1,2));
    title(figurename1)
    xlabel('\sigma_{1,winter} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
    
    subplot(222)
    plot(HF(i,:)./1e3,ALL(i,:)./1e3,'.','Color',[0.2 0.55 0.2]); hold on
    figurename1 = sprintf('HF (time = 2011 169.0) (R = %5.4f)',R_HF(1,2));
    title(figurename1)
    xlabel('\sigma_{1,HF} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')


    subplot(223)
    plot(BO(i,:)./1e3,ALL(i,:)./1e3,'.','Color',handle); hold on
    figurename1 = sprintf('bedopen (time = 2011 169.0) (R = %5.4f)',R_BO(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedopen} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')

    subplot(224)
    plot(BS(i,:)./1e3,ALL(i,:)./1e3,'.','Color',sky_blue); hold on
    figurename1 = sprintf('bedslip (time = 2011 169.0) (R = %5.4f)',R_BS(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedslip} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
end
% figurename2=sprintf('paperfig8_decomp_regression_surfC_check_169_20240208.png');
% print(gcf,'-dpng','-r600',figurename2);
%% 2011 169.5
for i=783 % 169.5
    % W vs. ALL
    [R_w,p_w] = corrcoef(W,ALL(i,:));
    W_ALL_stats(i,1) = R_w(1,2); W_ALL_stats(i,2) = p_w(1,2);
    
    % HF vs. ALL
    [R_HF,p_HF] = corrcoef(HF(i,:),ALL(i,:));
    HF_ALL_stats(i,1) = R_HF(1,2); HF_ALL_stats(i,2) = p_HF(1,2);
    
    % BS vs. ALL
    [R_BS,p_BS] = corrcoef(BS(i,:),ALL(i,:));
    BS_ALL_stats(i,1) = R_BS(1,2); BS_ALL_stats(i,2) = p_BS(1,2);

    % BO vs. ALL
    [R_BO,p_BO] = corrcoef(BO(i,:),ALL(i,:));
    BO_ALL_stats(i,1) = R_BO(1,2); BO_ALL_stats(i,2) = p_BO(1,2);

    figure(3) % = figure('Units','centimeters','Position',[2 1 20 20]);
    subplot(221)
    plot(W./1e3,ALL(i,:)./1e3,'.','Color',metal); hold on
    figurename1 = sprintf('winter (time = 2011 169.5) (R = %5.4f)',R_w(1,2));
    title(figurename1)
    xlabel('\sigma_{1,winter} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
    
    subplot(222)
    plot(HF(i,:)./1e3,ALL(i,:)./1e3,'.','Color',[0.2 0.55 0.2]); hold on
    figurename1 = sprintf('HF (time = 2011 169.5) (R = %5.4f)',R_HF(1,2));
    title(figurename1)
    xlabel('\sigma_{1,HF} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')


    subplot(223)
    plot(BO(i,:)./1e3,ALL(i,:)./1e3,'.','Color',handle); hold on
    figurename1 = sprintf('bedopen (time = 2011 169.5) (R = %5.4f)',R_BO(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedopen} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')

    subplot(224)
    plot(BS(i,:)./1e3,ALL(i,:)./1e3,'.','Color',sky_blue); hold on
    figurename1 = sprintf('bedslip (time = 2011 169.5) (R = %5.4f)',R_BS(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedslip} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
end
% figurename3=sprintf('paperfig8_decomp_regression_surfC_check_1695_20240208.png');
% print(gcf,'-dpng','-r600',figurename3);
%% 2011 170.0
for i=1071 % 170.0
    % W vs. ALL
    [R_w,p_w] = corrcoef(W,ALL(i,:));
    W_ALL_stats(i,1) = R_w(1,2); W_ALL_stats(i,2) = p_w(1,2);
    
    % HF vs. ALL
    [R_HF,p_HF] = corrcoef(HF(i,:),ALL(i,:));
    HF_ALL_stats(i,1) = R_HF(1,2); HF_ALL_stats(i,2) = p_HF(1,2);
    
    % BS vs. ALL
    [R_BS,p_BS] = corrcoef(BS(i,:),ALL(i,:));
    BS_ALL_stats(i,1) = R_BS(1,2); BS_ALL_stats(i,2) = p_BS(1,2);

    % BO vs. ALL
    [R_BO,p_BO] = corrcoef(BO(i,:),ALL(i,:));
    BO_ALL_stats(i,1) = R_BO(1,2); BO_ALL_stats(i,2) = p_BO(1,2);

    figure(4) % = figure('Units','centimeters','Position',[2 1 20 20]);
    subplot(221)
    plot(W./1e3,ALL(i,:)./1e3,'.','Color',metal); hold on
    figurename1 = sprintf('winter (time = 2011 170.0) (R = %5.4f)',R_w(1,2));
    title(figurename1)
    xlabel('\sigma_{1,winter} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
    
    subplot(222)
    plot(HF(i,:)./1e3,ALL(i,:)./1e3,'.','Color',[0.2 0.55 0.2]); hold on
    figurename1 = sprintf('HF (time = 2011 170.0) (R = %5.4f)',R_HF(1,2));
    title(figurename1)
    xlabel('\sigma_{1,HF} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')


    subplot(223)
    plot(BO(i,:)./1e3,ALL(i,:)./1e3,'.','Color',handle); hold on
    figurename1 = sprintf('bedopen (time = 2011 170.0) (R = %5.4f)',R_BO(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedopen} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')

    subplot(224)
    plot(BS(i,:)./1e3,ALL(i,:)./1e3,'.','Color',sky_blue); hold on
    figurename1 = sprintf('bedslip (time = 2011 170.0) (R = %5.4f)',R_BS(1,2));
    title(figurename1)
    xlabel('\sigma_{1,bedslip} value [ kPa ]')
    ylabel('\sigma_{1} value [ kPa ]')
end
% figurename4=sprintf('paperfig8_decomp_regression_surfC_check_170_20240208.png');
% print(gcf,'-dpng','-r600',figurename4);

%% 2012 stress decompositions (20 km x 20 km surface) to (6 km x 6 km)
% HF
load sigma_ij_2012_SLder_HF_surfC.mat
HF_2012 = sigma_ij_2012_SLder_HF.sigma1(:,ind_vec); % size = 764 X 625
clear sigma_ij_2012_SLder_HF
% bedslip
load sigma_ij_2012_SLder_bedslip_surfC.mat
BS_2012 = sigma_ij_2012_SLder_bedslip.sigma1(:,ind_vec);
clear sigma_ij_2012_SLder_bedslip
% bedopen
load sigma_ij_2012_SLder_bedopen_surfC.mat
BO_2012 = sigma_ij_2012_SLder_bedopen.sigma1(:,ind_vec);
clear sigma_ij_2012_SLder_bedopen
% total 
load n2012_allfour_surfC.mat
ALL_2012 = n2012_allfour.sigma1(:,ind_vec);
clear n2012_allfour

%% Pearson's correlation coefficient r (from -1 to +1) and p-values
for i=1:length(time_2012)
    % W vs. ALL
    [R_w,p_w] = corrcoef(W,ALL_2012(i,:));
    W_ALL_stats_2012(i,1) = R_w(1,2); W_ALL_stats_2012(i,2) = p_w(1,2);
    
    % HF vs. ALL
    [R_HF,p_HF] = corrcoef(HF_2012(i,:),ALL_2012(i,:));
    HF_ALL_stats_2012(i,1) = R_HF(1,2); HF_ALL_stats_2012(i,2) = p_HF(1,2);
    
    % BS vs. ALL
    [R_BS,p_BS] = corrcoef(BS_2012(i,:),ALL_2012(i,:));
    BS_ALL_stats_2012(i,1) = R_BS(1,2); BS_ALL_stats_2012(i,2) = p_BS(1,2);

    % BO vs. ALL
    [R_BO,p_BO] = corrcoef(BO_2012(i,:),ALL_2012(i,:));
    BO_ALL_stats_2012(i,1) = R_BO(1,2); BO_ALL_stats_2012(i,2) = p_BO(1,2);
end

% r values for (p-values < 0.01) turn into NaNs
threshold = 0.01; % 99% significance
I = find(W_ALL_stats_2012(:,2) >= threshold); W_ALL_stats_2012(:,3) = W_ALL_stats_2012(:,1); W_ALL_stats_2012(I,3) = NaN;
I = find(HF_ALL_stats_2012(:,2) >= threshold); HF_ALL_stats_2012(:,3) = HF_ALL_stats_2012(:,1); HF_ALL_stats_2012(I,3) = NaN;
I = find(BS_ALL_stats_2012(:,2) >= threshold); BS_ALL_stats_2012(:,3) = BS_ALL_stats_2012(:,1); BS_ALL_stats_2012(I,3) = NaN;
I = find(BO_ALL_stats_2012(:,2) >= threshold); BO_ALL_stats_2012(:,3) = BO_ALL_stats_2012(:,1); BO_ALL_stats_2012(I,3) = NaN;

%% Paper Figure correlation coefficient through time   

m=13; % fontsize

fig1 = figure('Units','centimeters','Position',0.18.*[2 1 85 120]);
    clf
    axe1 = axes('Position',[0.1 0.57 0.87 0.40],'Box','on','NextPlot','add');
    axe2 = axes('Position',[0.1 0.06 0.87 0.40],'Box','on','NextPlot','add'); 

% 2011
axes(axe1)
text(168.35-0.22, 1.05, 'a','FontSize',10,'FontWeight','bold'); hold on;

plot(time_2011,W_ALL_stats(:,3),'-','LineWidth',1.8,'Color',metal); hold on;
plot(time_2011,HF_ALL_stats(:,3),'-','LineWidth',1.8,'Color',[0.2 0.55 0.2]); hold on;
plot(time_2011,BO_ALL_stats(:,3),'-','LineWidth',1.8,'Color',handle); hold on;
plot(time_2011,BS_ALL_stats(:,3),'-','LineWidth',1.8,'Color',sky_blue); hold on;

plot([168.85 168.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);

text(168.85-0.03,1.05,'{\itt}_{1}','Rotation',0,'FontSize',10,'FontName','Avenir'); 
text(169.21-0.03,1.05,'{\itt}_{2}','Rotation',0,'FontSize',10,'FontName','Avenir'); 
text(169.32-0.03,1.05,'{\itt}_{3}','Rotation',0,'FontSize',10,'FontName','Avenir'); 
text(170.32-0.03,1.05,'{\itt}_{4}','Rotation',0,'FontSize',10,'FontName','Avenir'); 

plot([169.32+(6/24) 169.32+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
text(169.32-0.03+(5/24),1.05,'{\it\tau} = {6 hr}','Rotation',0,'FontSize',10,'FontName','Avenir','Color',[0.2 0.55 0.2]); 

plot(time_2011,W_ALL_stats(:,3),'-','LineWidth',1.8,'Color',metal); hold on;
plot(time_2011,HF_ALL_stats(:,3),'-','LineWidth',1.8,'Color',[0.2 0.55 0.2]); hold on;
plot(time_2011,BO_ALL_stats(:,3),'-','LineWidth',1.8,'Color',handle); hold on;
plot(time_2011,BS_ALL_stats(:,3),'-','LineWidth',1.8,'Color',sky_blue); hold on;

lgd = legend('\sigma_{1,winter}','\sigma_{1,HF}','\sigma_{1,bedopen}','\sigma_{1,bedslip}');
legend boxoff
set(lgd,'Position',[0.82 0.618 0.1 0.03],'FontSize',m-3);

ylim([-0.5 1]);   xlim([170.32-2 170.32]);
ylabel(['Correlation coefficient, {\it r},  when {\it p} ' char(8804) ' 0.01'],...
    'FontName','Avenir','FontSize',m);
xlabel('Day of Year, 2011  [ UTC ]','FontName','Avenir','FontSize',m);
set(gca,'tickdir','in','LineWidth',1.05,'FontSize',10,'FontName','Avenir'); 
set(gca,'FontName','Avenir','LineWidth',1.05,'xtick',[168.5:0.5:170]);
grid on 
hold all

% 2012
axes(axe2)
text(160.88-0.22, 1.05, 'b','FontSize',10,'FontWeight','bold'); hold on;

plot(time_2012,W_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',metal); hold on;
plot(time_2012,HF_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',[0.2 0.55 0.2]); hold on;
plot(time_2012,BO_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',handle); hold on;
plot(time_2012,BS_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',sky_blue); hold on;

plot([161.20 161.20],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([162.85 162.85],[-2000 4000],'-','LineWidth',1,'Color',[0.4 0.4 0.4]);

text(161.20-0.03,1.05,'{\itt}_{1}','Rotation',0,'FontSize',10,'FontName','Avenir'); 
text(161.72-0.03,1.05,'{\itt}_{2}','Rotation',0,'FontSize',10,'FontName','Avenir'); 
text(161.85-0.03,1.05,'{\itt}_{3}','Rotation',0,'FontSize',10,'FontName','Avenir'); 
text(162.85-0.03,1.05,'{\itt}_{4}','Rotation',0,'FontSize',10,'FontName','Avenir'); 

plot([161.85+(6/24) 161.85+(6/24)],[-2000 4000],'-','LineWidth',1.2,'Color',[0.2 0.55 0.2]); %maxwell
text(161.85-0.03+(5/24),1.05,'{\it\tau} = {6 hr}','Rotation',0,'FontSize',10,'FontName','Avenir','Color',[0.2 0.55 0.2]); 

plot(time_2012,W_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',metal); hold on;
plot(time_2012,HF_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',[0.2 0.55 0.2]); hold on;
plot(time_2012,BO_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',handle); hold on;
plot(time_2012,BS_ALL_stats_2012(:,3),'-','LineWidth',1.8,'Color',sky_blue); hold on;

lgd2 = legend('\sigma_{1,winter}','\sigma_{1,HF}','\sigma_{1,bedopen}','\sigma_{1,bedslip}');
legend boxoff
set(lgd2,'Position',[0.82 0.118 0.1 0.03],'FontSize',m-3);
 
ylim([-0.5 1]);  xlim([162.85-2 162.85]);
ylabel(['Correlation coefficient, {\it r},  when {\it p} ' char(8804) ' 0.01'],...
    'FontName','Avenir','FontSize',m);
xlabel('Day of Year, 2012  [ UTC ]','FontName','Avenir','FontSize',m);
set(gca,'tickdir','in','LineWidth',1.05,'FontSize',10,'FontName','Avenir'); 
set(gca,'FontName','Avenir','LineWidth',1.05,'xtick',[161:0.5:162.5]);
grid on 

% PRINT FIGURE
%figurename=sprintf('paperfig7_decomp_regression_surfC_20240208.png');
%print(gcf,'-dpng','-r600',figurename);
