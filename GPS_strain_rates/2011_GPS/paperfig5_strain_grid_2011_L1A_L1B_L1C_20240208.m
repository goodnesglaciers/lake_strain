%% strain rate plots for 2011 drainage with Track position outputs archived in neu_nif_2011.mat
% June 2018: reimagined for proposal, added strain rates for L1A, L1C
% May 2019: Added a little map
% Oct 2020: Added duration of elevated positive strain rates
% Oct 2021: Corrected strain-rate error estimates based on Jonny's derivation
% June 2023: Figure updates from reviews, new supplementary figure on
% strain rates across two basins
% February 2024: replace yellow, 1-sigma errors from TRACK

%% interGPS station distance, station locations in 2011 
clear all; close all;
load('S_neu_err_nif_2011.mat'); % load TRACK NEU position archive

%% interpolate to get consistent time vector, smooth over time window
interptime = 0.0138;      % 20 minutes [decimal day]
xi12 = (166.85:interptime:171.0)'; yi12 = xi12; zi12 = yi12; 
length_xi12 = length(xi12);
length_xi12_minus1 = length_xi12-1;

Q = SS; % include interpolated values within SS archive
nfQ = 14; % number of stations

timestep = nanmean(diff(Q(1).time_new2(4000:25000))); % decimal days
timestep_minutes = timestep*24*60; % minutes

span = 120;  % window width [points]
timestep_smooth = span.*timestep_minutes; % window width [minutes]

%return

% work through stations; avoid later portions of the timeseries
AA = horzcat(Q(1).time_new2(4000:25000), smooth(Q(1).e12(4000:25000,1),span),...
   smooth(Q(1).n12(4000:25000,1),span),...
   smooth(Q(1).u12(4000:25000,1),span) , ...
   smooth(Q(1).e12(4000:25000,2),span), smooth(Q(1).n12(4000:25000,2),span), ...
   smooth(Q(1).u12(4000:25000,2),span));
AA1 = sortrows(AA,1);
F=griddedInterpolant(AA1(:,1),AA1(:,2));
Q(1).e12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,3));
Q(1).n12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,4));
Q(1).u12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,5));
Q(1).de12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,6));
Q(1).dn12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,7));
Q(1).du12i = F(xi12);

for i=2:1:4
    Q(i,1).e12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,1),span),xi12);
    Q(i,1).n12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,1),span),yi12);
    Q(i,1).u12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,1),span),zi12);
    Q(i,1).de12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,2),span),xi12);
    Q(i,1).dn12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,2),span),yi12);
    Q(i,1).du12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,2),span),zi12);
end

AA = horzcat(Q(5).time_new2(4000:25000), smooth(Q(5).e12(4000:25000,1),span),...
   smooth(Q(5).n12(4000:25000,1),span),...
   smooth(Q(5).u12(4000:25000,1),span) , ...
   smooth(Q(5).e12(4000:25000,2),span), smooth(Q(5).n12(4000:25000,2),span), ...
   smooth(Q(5).u12(4000:25000,2),span));
AA1 = sortrows(AA,1);
F=griddedInterpolant(AA1(:,1),AA1(:,2));
Q(5).e12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,3));
Q(5).n12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,4));
Q(5).u12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,5));
Q(5).de12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,6));
Q(5).dn12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,7));
Q(5).du12i = F(xi12);

AA = horzcat(Q(6).time_new2(4000:25000), smooth(Q(6).e12(4000:25000,1),span),...
   smooth(Q(6).n12(4000:25000,1),span),...
   smooth(Q(6).u12(4000:25000,1),span) , ...
   smooth(Q(6).e12(4000:25000,2),span), smooth(Q(6).n12(4000:25000,2),span), ...
   smooth(Q(6).u12(4000:25000,2),span));
AA1 = sortrows(AA,1);
F=griddedInterpolant(AA1(:,1),AA1(:,2));
Q(6).e12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,3));
Q(6).n12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,4));
Q(6).u12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,5));
Q(6).de12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,6));
Q(6).dn12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,7));
Q(6).du12i = F(xi12);

for i=7
    Q(i,1).e12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,1),span),xi12);
    Q(i,1).n12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,1),span),yi12);
    Q(i,1).u12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,1),span),zi12);
    Q(i,1).de12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,2),span),xi12);
    Q(i,1).dn12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,2),span),yi12);
    Q(i,1).du12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,2),span),zi12);
end

AA = horzcat(Q(8).time_new2(4000:25000), smooth(Q(8).e12(4000:25000,1),span),...
   smooth(Q(8).n12(4000:25000,1),span),...
   smooth(Q(8).u12(4000:25000,1),span) , ...
   smooth(Q(8).e12(4000:25000,2),span), smooth(Q(8).n12(4000:25000,2),span), ...
   smooth(Q(8).u12(4000:25000,2),span));
AA1 = sortrows(AA,1);
F=griddedInterpolant(AA1(:,1),AA1(:,2));
Q(8).e12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,3));
Q(8).n12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,4));
Q(8).u12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,5));
Q(8).de12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,6));
Q(8).dn12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,7));
Q(8).du12i = F(xi12);

AA = horzcat(Q(9).time_new2(4000:25000), smooth(Q(9).e12(4000:25000,1),span),...
   smooth(Q(9).n12(4000:25000,1),span),...
   smooth(Q(9).u12(4000:25000,1),span) , ...
   smooth(Q(9).e12(4000:25000,2),span), smooth(Q(9).n12(4000:25000,2),span), ...
   smooth(Q(9).u12(4000:25000,2),span));
AA1 = sortrows(AA,1);
F=griddedInterpolant(AA1(:,1),AA1(:,2));
Q(9).e12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,3));
Q(9).n12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,4));
Q(9).u12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,5));
Q(9).de12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,6));
Q(9).dn12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,7));
Q(9).du12i = F(xi12);

for i=10:11
    Q(i,1).e12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,1),span),xi12);
    Q(i,1).n12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,1),span),yi12);
    Q(i,1).u12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,1),span),zi12);
    Q(i,1).de12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,2),span),xi12);
    Q(i,1).dn12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,2),span),yi12);
    Q(i,1).du12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,2),span),zi12);
end

AA = horzcat(Q(12).time_new2(4000:25000), smooth(Q(12).e12(4000:25000,1),span),...
   smooth(Q(12).n12(4000:25000,1),span),...
   smooth(Q(12).u12(4000:25000,1),span) , ...
   smooth(Q(12).e12(4000:25000,2),span), smooth(Q(12).n12(4000:25000,2),span), ...
   smooth(Q(12).u12(4000:25000,2),span));
AA1 = sortrows(AA,1);
F=griddedInterpolant(AA1(:,1),AA1(:,2));
Q(12).e12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,3));
Q(12).n12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,4));
Q(12).u12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,5));
Q(12).de12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,6));
Q(12).dn12i = F(xi12);
F=griddedInterpolant(AA1(:,1),AA1(:,7));
Q(12).du12i = F(xi12);

for i=13:14
    Q(i,1).e12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,1),span),xi12);
    Q(i,1).n12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,1),span),yi12);
    Q(i,1).u12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,1),span),zi12);
    Q(i,1).de12i=interp1(Q(i).time_new2,smooth(Q(i).e12(:,2),span),xi12);
    Q(i,1).dn12i=interp1(Q(i).time_new2,smooth(Q(i).n12(:,2),span),yi12);
    Q(i,1).du12i=interp1(Q(i).time_new2,smooth(Q(i).u12(:,2),span),zi12);
end

%% find positions in between stations
Loc=[];
failedFiles={};
for i=1:1:nfQ
    for j=1:1:nfQ
        try
        Loc.x(i,j)=((Q(i,1).x+Q(j,1).x)/2);
        Loc.y(i,j)=((Q(i,1).y+Q(j,1).y)/2);
            catch 
    failedFiles = [failedFiles; Q(i,1).site];
    continue
    end
    end
end
Loc.x=triu(Loc.x); Loc.y=triu(Loc.y);

% find velocity of stations
for i=1:1:nfQ
    Q(i,1).xvel=(gradient(Q(i,1).e12i))./interptime; % m/day
    Q(i,1).yvel=(gradient(Q(i,1).n12i))./interptime; % m/day
    Q(i,1).zvel=(gradient(Q(i,1).u12i))./interptime; % m/day
end

% difference velocity of station
for i=1:1:nfQ
    for j=1:1:nfQ
        for t=1:length_xi12_minus1
        xvel_diff(i,j,t)=(Q(i,1).xvel(t))-(Q(j,1).xvel(t));
        yvel_diff(i,j,t)=(Q(i,1).yvel(t))-(Q(j,1).yvel(t));
        zvel_diff(i,j,t)=(Q(i,1).zvel(t))-(Q(j,1).zvel(t));
        end
    end
end

% upper triangle: select only upper triangle of matrix U = triu(X)
lengthvel=length(xvel_diff)
for i=1:lengthvel
    du12(:,:,i)=triu(xvel_diff(:,:,i));
    dv12(:,:,i)=triu(yvel_diff(:,:,i));
    dw12(:,:,i)=triu(zvel_diff(:,:,i));
end

% assign XY map values to movements = movement in Euclidean space
failedFiles={};
for i=1:1:nfQ
    try
    Q(i,:).eloc12=Q(i,1).x+(Q(i,:).e12i)-(Q(i,:).e12i(1,1));
    Q(i,:).nloc12=Q(i,1).y+(Q(i,:).n12i)-(Q(i,:).n12i(1,1));
    Q(i,:).uloc12=(Q(i,1).u12i)-(Q(i,1).u12i(1,1));
        catch 
    failedFiles = [failedFiles; Q(i,1).site];
    continue
    end
end
failedFiles

% difference position of stations
for i=1:1:nfQ
    for j=1:1:nfQ
        for t=1:length_xi12
        xpos_diff(i,j,t)=Q(i,:).eloc12(t)-Q(j,:).eloc12(t);
        ypos_diff(i,j,t)=Q(i,:).nloc12(t)-Q(j,:).nloc12(t);
        zpos_diff(i,j,t)=Q(i,:).uloc12(t)-Q(j,:).uloc12(t);
        end
    end
end

% upper triangle: select only upper triangle of matrix U = triu(X)
lengthpos=length(xpos_diff)
for i=1:lengthpos
    dx12(:,:,i)=triu(xpos_diff(:,:,i));
    dy12(:,:,i)=triu(ypos_diff(:,:,i));
    dz12(:,:,i)=triu(zpos_diff(:,:,i));
end
dx12=dx12(:,:,2:end);dy12=dy12(:,:,2:end);dz12=dz12(:,:,2:end);

%% change to along- and across flowline
flowline = atand(nanmean(dx12(1,2,:))./nanmean(dy12(1,2,:)));
dx12_f = dx12.*sind(flowline) + dy12.*cosd(flowline);
dy12_f = dy12.*sind(flowline) + dx12.*cosd(flowline);
du12_f = du12.*sind(flowline) + dv12.*cosd(flowline);
dv12_f = dv12.*sind(flowline) + du12.*cosd(flowline);

%% STRAIN RATES
% calculate strain rates (velocity over length) day^{-1}
dudx=du12./dx12;
dvdy=dv12./dy12;
dwdz=dw12./dz12;
shear=0.5.*(((dv12)./dx12)+((du12)./dy12)); % shear
length_triu=length(dudx);

% calculate strain rates (velocity over length) year^{-1}
dudx_yr=(du12.*365.25)./dx12; % east west (~longitudinal)
dvdy_yr=(dv12.*365.25)./dy12; % north south (~transverse)
dwdz_yr=(dw12.*365.25)./dz12; % vertical
shear_approx_yr=0.5.*(((dv12.*365.25)./dx12)+((du12.*365.25)./dy12)); % shear 

% calculate strain rates (velocity over length) year^{-1} in direction of
% flowline
du12_f_year = du12_f.*365.25; 
dv12_f_year = dv12_f.*365.25;
lon_yr = (du12_f.*365.25)./dx12_f; % longitudinal
trans_yr = (dv12_f.*365.25)./dy12_f; % transverse
shear_yr=0.5.*(((dv12_f.*365.25)./dx12_f)+((du12_f.*365.25)./dy12_f)); % shear

%% ERRORS IN STRAIN RATES -- correct from Jonny
delta = 0.02; % assumes all Track errors are 2 cm
delta_sqr = delta*delta; % assumes all Track errors are 2 cm

% this is all in units of day^{-1} 
for i=1:1:nfQ
    for j=1:1:nfQ
        for t=2:length_xi12-2
           % delta value, each station
           % delta x
           delta_sqr_x1_Tbehind = Q(i).de12i(t-1)^2;
           delta_sqr_x1_Tforward = Q(i).de12i(t+1)^2;
           delta_sqr_x2_Tbehind = Q(j).de12i(t-1)^2;
           delta_sqr_x2_Tforward = Q(j).de12i(t+1)^2;
           % delta y
           delta_sqr_y1_Tbehind = Q(i).dn12i(t-1)^2;
           delta_sqr_y1_Tforward = Q(i).dn12i(t+1)^2;
           delta_sqr_y2_Tbehind = Q(j).dn12i(t-1)^2;
           delta_sqr_y2_Tforward = Q(j).dn12i(t+1)^2;
           % delta_Exx             
           delta_Exx_yr(i,j,t) = abs(lon_yr(i,j,t)).*sqrt(((1./(du12_f_year(i,j,t).^2)).*...
               (delta_sqr_x2_Tforward+delta_sqr_x2_Tbehind+delta_sqr_x1_Tforward+delta_sqr_x1_Tbehind))+...
               ((1./(dx12_f(i,j,t).^2)).*(delta_sqr_x1_Tforward+delta_sqr_x2_Tforward)));
           % delta_Eyy             
           delta_Eyy_yr(i,j,t) = abs(trans_yr(i,j,t)).*sqrt(((1./(dv12_f_year(i,j,t).^2)).*...
               (delta_sqr_y2_Tforward+delta_sqr_y2_Tbehind+delta_sqr_y1_Tforward+delta_sqr_y1_Tbehind))+...
               ((1./(dy12_f(i,j,t).^2)).*(delta_sqr_y1_Tforward+delta_sqr_y2_Tforward)));
           % delta_Exy
           delta_Exy_yr(i,j,t) = abs(shear_yr(i,j,t)).*(sqrt(((1./(dv12_f_year(i,j,t).^2)).*...
               (delta_sqr_x2_Tforward+delta_sqr_x2_Tbehind+delta_sqr_x1_Tforward+delta_sqr_x1_Tbehind))+...
               ((1./(dy12_f(i,j,t).^2)).*(delta_sqr_x1_Tforward+delta_sqr_x2_Tforward)) + ...
               ((1./(du12_f_year(i,j,t).^2)).*...
               (delta_sqr_y2_Tforward+delta_sqr_y2_Tbehind+delta_sqr_y1_Tforward+delta_sqr_y1_Tbehind))+...
               ((1./(dx12_f(i,j,t).^2)).*(delta_sqr_y1_Tforward+delta_sqr_y2_Tforward))));  
        end
    end
end

%% scrap figure for colors
Fig4 = figure(4); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 19 14]);
axe1 = axes('Position',[0.075 0.70 0.9 0.25],'Box','on','NextPlot','add','XTickLabels',[]);

axes(axe1)
h1=plot([168.85 168.85], [-100 100],'-'); c1 = get(h1,'Color'); hold on;
h2=plot([169.85 168.85], [-100 100],'-'); c2 = get(h2,'Color');
h3=plot([170.85 168.85], [-100 100],'-'); c3 = get(h3,'Color');
h4=plot([171.85 168.85], [-100 100],'-'); c3 = get(h4,'Color');
h5=plot([172.85 168.85], [-100 100],'-'); c4 = get(h5,'Color');
h6=plot([173.85 168.85], [-100 100],'-'); c5 = get(h6,'Color');
h7=plot([174.85 168.85], [-100 100],'-'); c6 = get(h7,'Color');
h8=plot([175.85 168.85], [-100 100],'-'); c7 = get(h8,'Color');
h9=plot([176.85 168.85], [-100 100],'-'); c8 = get(h9,'Color');

%% load L1A geographic files
origin = [68.72, -49.53];
load apcoords_lle
    lats=apcoords_lle(1,:); lons=apcoords_lle(2,:); hs=apcoords_lle(3,:);
    llh=[lats; lons; hs];
    xy_sta_11=llh2localxy(llh,origin);
load lake.mat
    llh_lake = [lake(:,5)'; lake(:,4)'; zeros(337,1)'];
    xy_lake = llh2localxy(llh_lake,origin);
load lil_lake_2013_168.mat
    llh_lil_lake = [lil_lake_2013_168(:,4)'; lil_lake_2013_168(:,3)'; zeros(148,1)'];
    xy_lil_lake = llh2localxy(llh_lil_lake,origin); 
% polarstereo conversion needed values
    radius=6378137.0;    eccen=0.08181919;    lat_true=70;    lon_posy=-45;
    NNL=csvread('NNL_20110617.csv',1,0); 
    [phi,lambda]=polarstereo_inv(NNL(:,1),NNL(:,2),radius,eccen,lat_true,lon_posy);
    llh_nnl = [phi,lambda]; 
    xy_nnl_lake = llh2localxy(llh_nnl',origin);
% NNL crack from 2016/222
    NNL2016=csvread('nnl2016222_latlon.csv',1,0);
    xy_NNL_crack = llh2localxy(horzcat(NNL2016(:,2), NNL2016(:,1))',origin);
    strike_NNL = 90 - atand((xy_NNL_crack(end,2)-xy_NNL_crack(1,2))./...
        (xy_NNL_crack(end,1)-xy_NNL_crack(1,1)) )  %  for plotting purposes
% SNL crack from 2011 June 17
    SNL=csvread('snl_crack2_2011June17.csv',1,0);
    [phi,lambda]=polarstereo_inv(SNL(:,1),SNL(:,2),radius,eccen,lat_true,lon_posy);
    xy_SNL_crack = llh2localxy([phi,lambda]',origin); 
    strike_SNL = atand((xy_SNL_crack(end,2)-xy_SNL_crack(1,2))./...
        (xy_SNL_crack(end,1)-xy_SNL_crack(1,1)) ) 
% strike NL crack in lake    
    load out2011.mat
    patchesC = out2011.patches_C;
    x_crack = patchesC(1:6:end,6); y_crack = patchesC(1:6:end,7);
    strike_NL = atand((y_crack(22,1)-y_crack(8,1))./(x_crack(22,1)-x_crack(8,1))) +180  % +180 for plotting purposes       

%% NL1A—C strain rates
Fig2 = figure(2); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[0 5 18*2 20]);
axe0 = axes('Position',[0.03 0.70 0.15 0.2625],'Box','on','NextPlot','add','XTickLabels',[]);
axe01 = axes('Position',[0.03 0.39 0.15 0.2625],'Box','on','NextPlot','add','XTickLabels',[]);
axe02 = axes('Position',[0.03 0.08 0.15 0.2625],'Box','on','NextPlot','add');
axe1 = axes('Position',[0.24 0.7 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe2 = axes('Position',[0.24 0.39 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe3 = axes('Position',[0.24 0.08 0.35 0.26],'Box','on','NextPlot','add'); 
axe4 = axes('Position',[0.645 0.7 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe5 = axes('Position',[0.645 0.39 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe6 = axes('Position',[0.645 0.08 0.35 0.26],'Box','on','NextPlot','add'); 

sites = {'FL03','FL04','NL01','NL02','NL03','NL04','NL06','NL07','NL08','NL09',...
    'NL10','NL11','NL12','NL13','NLBS'};

% dock at eel pond 
sky_blue = [111, 169, 228]./255; % SNL
metal = [87, 115, 131]./255; % NNL
oar = [251, 219, 154]./255; % NENL
handle = [161, 37, 49]./255; % NL
dark_oar = [164, 114, 63]./255;

axes(axe0);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'a.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, -0.50, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(3,1) xy_sta_11(5,1)],[xy_sta_11(3,2) xy_sta_11(5,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(4,1) xy_sta_11(5,1)],[xy_sta_11(4,2) xy_sta_11(5,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(3,1) xy_sta_11(7,1)],[xy_sta_11(3,2) xy_sta_11(7,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_11(6,1) xy_sta_11(5,1)],[xy_sta_11(6,2) xy_sta_11(5,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_11(3,1) xy_sta_11(6,1)],[xy_sta_11(3,2) xy_sta_11(6,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_11(4,1) xy_sta_11(6,1)],[xy_sta_11(4,2) xy_sta_11(6,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_11(3:4,1)-1.80,xy_sta_11(3:4,2)-0.10,sites(3:4),'FontSize',9,'Rotation',0);
text(xy_sta_11(5,1)+0.40,xy_sta_11(5,2)-0.10,sites(5),'FontSize',9,'Rotation',0);
text(xy_sta_11(6:7,1)-0.70,xy_sta_11(6:7,2)-0.50,sites(6:7),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel(' y [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe01);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'b.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(2,1) xy_sta_11(8,1)],[xy_sta_11(2,2) xy_sta_11(8,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(2,1) xy_sta_11(6,1)],[xy_sta_11(2,2) xy_sta_11(6,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(6,1) xy_sta_11(10,1)],[xy_sta_11(6,2) xy_sta_11(10,2)],'-','Color',c3,'LineWidth',1.5)
plot([xy_sta_11(7,1) xy_sta_11(9,1)],[xy_sta_11(7,2) xy_sta_11(9,2)],'-','Color',c4,'LineWidth',1.5)

% trans
plot([xy_sta_11(6,1) xy_sta_11(9,1)],[xy_sta_11(6,2) xy_sta_11(9,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_11(7,1) xy_sta_11(10,1)],[xy_sta_11(7,2) xy_sta_11(10,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(2,1)-0.55,xy_sta_11(2,2)-0.50,sites(2),'FontSize',9,'Rotation',0);
text(xy_sta_11(8,1)-0.75,xy_sta_11(8,2)-0.5,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_11(6:7,1)-0.65,xy_sta_11(6:7,2)+0.55,sites(6:7),'FontSize',9,'Rotation',0);
text(xy_sta_11(9,1)-0.35,xy_sta_11(9,2)-0.5,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_11(10,1)-0.65,xy_sta_11(10,2)-0.5,sites(10),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe02);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'c.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(8,1) xy_sta_11(12,1)],[xy_sta_11(8,2) xy_sta_11(12,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(11,1) xy_sta_11(12,1)],[xy_sta_11(11,2) xy_sta_11(12,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(10,1) xy_sta_11(11,1)],[xy_sta_11(10,2) xy_sta_11(11,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_11(9,1) xy_sta_11(14,1)],[xy_sta_11(9,2) xy_sta_11(14,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_11(9,1) xy_sta_11(12,1)],[xy_sta_11(9,2) xy_sta_11(12,2)],'-','Color',c5,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(9:10,1)-0.45,xy_sta_11(9:10,2)+0.55,sites(9:10),'FontSize',9,'Rotation',0);
text(xy_sta_11(11,1)-1.5,xy_sta_11(11,2)-0.45,sites(11),'FontSize',9,'Rotation',0);
text(xy_sta_11(8,1)-1.5,xy_sta_11(8,2)-0.5,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_11(12,1)+0.3,xy_sta_11(12,2)-0.25,sites(12),'FontSize',9,'Rotation',0);
text(xy_sta_11(14,1)+0.3,xy_sta_11(14,2)+0.15,sites(14),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]'); xlabel('x [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;
xtickangle(0)

axes(axe1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(5,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(6,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(2,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3) % NL06 NL01
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL01-NL03','NL02-NL03','NL01-NL06'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, 'd. L1C','FontSize',12,'FontWeight','bold'); 
text(168.81,0.06,'Start of Precursor','Rotation',90,'FontSize',10);
text(169.17,0.06,'HF Initiation','Rotation',90,'FontSize',10);
text(169.35,0.11,'Max HF','Rotation',90,'FontSize',10);
text(170.27,-0.18,'+1 day','Rotation',90,'FontSize',10);

text(168.85-0.03,0.44,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.21-0.03,0.44,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.32-0.03,0.44,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(170.32-0.03,0.44,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.8],'FontName','Avenir');  grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');

axes(axe2)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(2,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(9,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1) % NL07 Fl04
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2) % NL04 Fl04
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL04-NL09','NL06-NL08','NL07-FL04','NL04-FL04'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, 'e. L1A','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.4],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');

axes(axe3)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(9,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(4,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(10,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL07-NL11','NL09-NL10','NL10-NL11','Location','NorthEast'); legend boxoff
ylim([-0.22 0.4]);  
xlim([168.32 170.32]);
 
text(168.33, 0.44, 'f. L1B','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.4],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');
xlabel(' Day of Year, 2011  [ UTC ]','FontName','Avenir','FontSize',12);

axes(axe4)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,6,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6)

plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL03-NL04','NL01-NL04','NL02-NL04'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, 'g. L1C','FontSize',12,'FontWeight','bold') 
text(168.81,0.06,'Start of Precursor','Rotation',90,'FontSize',10);
text(169.17,0.11,'HF Initiation','Rotation',90,'FontSize',10);
text(169.35,0.11,'Max HF','Rotation',90,'FontSize',10);
text(170.27,-0.18,'+1 day','Rotation',90,'FontSize',10);

text(168.85-0.03,0.44,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.21-0.03,0.44,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.32-0.03,0.44,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(170.32-0.03,0.44,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.6],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');

axes(axe5)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(2,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(2,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL04-NL08','NL06-NL09','NL04-NL09','NL06-NL08'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);

text(168.33, 0.44, 'h. L1A','FontSize',12,'FontWeight','bold')
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.6],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');

axes(axe6)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(9,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,12,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd = legend('NL07-NL11','NL08-NL11','NL08-NL13','Location','NorthEast'); legend boxoff
set(lgd,'Position',[0.89 0.27 0.1 0.05]);
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, 'i. L1B','FontSize',12,'FontWeight','bold')
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.6:0.1:0.6],'FontName','Avenir');
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');
xlabel(' Day of Year, 2011  [ UTC ]','FontSize',12,'FontName','Avenir');  grid on;

%print(gcf,'-dpng','-r500','paperfig5_strainrate_2011_20240208.png');

%% ERROR IN STRAIN RATE 2011 L1A-L1C
Fig3 = figure(3); clf; 

set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[0 5 18*2 20]);
axe0 = axes('Position',[0.03 0.70 0.15 0.2625],'Box','on','NextPlot','add','XTickLabels',[]);
axe01 = axes('Position',[0.03 0.39 0.15 0.2625],'Box','on','NextPlot','add','XTickLabels',[]);
axe02 = axes('Position',[0.03 0.08 0.15 0.2625],'Box','on','NextPlot','add');
axe1 = axes('Position',[0.24 0.7 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe2 = axes('Position',[0.24 0.39 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe3 = axes('Position',[0.24 0.08 0.35 0.26],'Box','on','NextPlot','add'); 
axe4 = axes('Position',[0.645 0.7 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe5 = axes('Position',[0.645 0.39 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe6 = axes('Position',[0.645 0.08 0.35 0.26],'Box','on','NextPlot','add'); 

axes(axe0);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'a.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% LAKES
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, -0.50, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(3,1) xy_sta_11(5,1)],[xy_sta_11(3,2) xy_sta_11(5,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(4,1) xy_sta_11(5,1)],[xy_sta_11(4,2) xy_sta_11(5,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(3,1) xy_sta_11(7,1)],[xy_sta_11(3,2) xy_sta_11(7,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_11(6,1) xy_sta_11(5,1)],[xy_sta_11(6,2) xy_sta_11(5,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_11(3,1) xy_sta_11(6,1)],[xy_sta_11(3,2) xy_sta_11(6,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_11(4,1) xy_sta_11(6,1)],[xy_sta_11(4,2) xy_sta_11(6,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_11(3:4,1)-1.80,xy_sta_11(3:4,2)-0.10,sites(3:4),'FontSize',9,'Rotation',0);
text(xy_sta_11(5,1)+0.40,xy_sta_11(5,2)-0.10,sites(5),'FontSize',9,'Rotation',0);
text(xy_sta_11(6:7,1)-0.70,xy_sta_11(6:7,2)-0.50,sites(6:7),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel(' y [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;


axes(axe01);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'b.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% LAKES
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(2,1) xy_sta_11(8,1)],[xy_sta_11(2,2) xy_sta_11(8,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(2,1) xy_sta_11(6,1)],[xy_sta_11(2,2) xy_sta_11(6,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(6,1) xy_sta_11(10,1)],[xy_sta_11(6,2) xy_sta_11(10,2)],'-','Color',c3,'LineWidth',1.5)
plot([xy_sta_11(7,1) xy_sta_11(9,1)],[xy_sta_11(7,2) xy_sta_11(9,2)],'-','Color',c4,'LineWidth',1.5)

% trans
plot([xy_sta_11(6,1) xy_sta_11(9,1)],[xy_sta_11(6,2) xy_sta_11(9,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_11(7,1) xy_sta_11(10,1)],[xy_sta_11(7,2) xy_sta_11(10,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(2,1)-0.55,xy_sta_11(2,2)-0.50,sites(2),'FontSize',9,'Rotation',0);
text(xy_sta_11(8,1)-0.75,xy_sta_11(8,2)-0.5,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_11(6:7,1)-0.65,xy_sta_11(6:7,2)+0.55,sites(6:7),'FontSize',9,'Rotation',0);
text(xy_sta_11(9,1)-0.35,xy_sta_11(9,2)-0.5,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_11(10,1)-0.65,xy_sta_11(10,2)-0.5,sites(10),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe02);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'c.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.75, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(8,1) xy_sta_11(12,1)],[xy_sta_11(8,2) xy_sta_11(12,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(11,1) xy_sta_11(12,1)],[xy_sta_11(11,2) xy_sta_11(12,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(10,1) xy_sta_11(11,1)],[xy_sta_11(10,2) xy_sta_11(11,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_11(9,1) xy_sta_11(14,1)],[xy_sta_11(9,2) xy_sta_11(14,2)],'-','Color',c4,'LineWidth',1.5)
%plot([xy_sta_11(8,1) xy_sta_11(12,1)],[xy_sta_11(8,2) xy_sta_11(12,2)],':','Color',oar,'LineWidth',2.5)
plot([xy_sta_11(9,1) xy_sta_11(12,1)],[xy_sta_11(9,2) xy_sta_11(12,2)],'-','Color',c5,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(9:10,1)-0.45,xy_sta_11(9:10,2)+0.55,sites(9:10),'FontSize',9,'Rotation',0);
text(xy_sta_11(11,1)-1.5,xy_sta_11(11,2)-0.45,sites(11),'FontSize',9,'Rotation',0);
text(xy_sta_11(8,1)-1.5,xy_sta_11(8,2)-0.5,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_11(12,1)+0.3,xy_sta_11(12,2)-0.25,sites(12),'FontSize',9,'Rotation',0);
text(xy_sta_11(14,1)+0.3,xy_sta_11(14,2)+0.15,sites(14),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]'); xlabel('x [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;
xtickangle(0)

axes(axe1)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(5,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(6,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(2,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3) % NL06 NL01
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL01-NL03','NL02-NL03','NL02-NL06'); legend boxoff
ylim([0 0.00006]);
xlim([168.32 170.32]);
text(168.35, 5.6*(10^-5), 'd. L1C','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001],'FontName','Avenir'); grid on
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('${\delta}\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);

text(168.85-0.03,0.000065,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.21-0.03,0.000065,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.32-0.03,0.000065,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(170.32-0.03,0.000065,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

axes(axe2)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(2,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(9,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1) % NL07 Fl04
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(1,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2) % NL04 Fl04
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL04-NL09','NL06-NL08','NL07-FL04','NL04-FL04'); legend boxoff
ylim([0 0.00006]);
xlim([168.32 170.32]);
text(168.35, 5.6*(10^-5), 'e. L1A','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.00004],'FontName','Avenir'); grid on
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('${\delta}\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);

axes(axe3)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(9,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(4,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(10,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL07-NL11','NL09-NL10','NL10-NL11','Location','NorthEast'); legend boxoff
ylim([0 0.00006]);  
xlim([168.32 170.32]);
text(168.35, 5.6*(10^-5), 'f. L1B','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001],'FontName','Avenir'); grid on
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('${\delta}\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
xlabel(' Day of Year, 2011  [ UTC ]','FontSize',12); 

axes(axe4)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,6,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL03-NL04','NL01-NL04','NL02-NL04'); legend boxoff
ylim([0 0.00006]);
xlim([168.32 170.32]);
text(168.35, 5.6*(10^-5), 'g. L1C','FontSize',12,'FontWeight','bold')
text(168.85-0.03,0.000065,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.21-0.03,0.000065,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.32-0.03,0.000065,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(170.32-0.03,0.000065,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001],'FontName','Avenir'); grid on
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('${\delta}\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);

axes(axe5)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(2,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(2,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL04-NL08','NL06-NL09','NL04-NL09','NL06-NL08'); legend boxoff
ylim([0 0.00006]); 
xlim([168.32 170.32]);
text(168.35, 5.6*(10^-5), 'h. L1A','FontSize',12,'FontWeight','bold')
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001],'FontName','Avenir'); grid on
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('${\delta}\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);

axes(axe6)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(9,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(3,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(3,12,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd = legend('NL07-NL11','NL08-NL11','NL08-NL13','Location','NorthEast'); legend boxoff
set(lgd,'Position',[0.89 0.27 0.1 0.05]);
ylim([0 0.00006]); 
xlim([168.32 170.32]); 
text(168.35, 5.6*(10^-5), 'i. L1B','FontSize',12,'FontWeight','bold')
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001],'FontName','Avenir'); grid on
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('${\delta}\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
xlabel(' Day of Year, 2011  [ UTC ]','FontSize',12);  grid on;

%print(gcf,'-dpng','-r500','paperfigS2_strainrate_errors_2011_20240208.png');

%% EXAMPLES OF STRAIN RATES CALCULATED OVER JUST NL, NL+NNL, NL+SNL
Fig4 = figure(4); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[0 5 18*2 20]);
axe0 = axes('Position',[0.03 0.70 0.15 0.2625],'Box','on','NextPlot','add','XTickLabels',[]);
axe01 = axes('Position',[0.03 0.39 0.15 0.2625],'Box','on','NextPlot','add','XTickLabels',[]);
axe02 = axes('Position',[0.03 0.08 0.15 0.2625],'Box','on','NextPlot','add');
axe1 = axes('Position',[0.24 0.7 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe2 = axes('Position',[0.24 0.39 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe3 = axes('Position',[0.24 0.08 0.35 0.26],'Box','on','NextPlot','add'); 
axe4 = axes('Position',[0.645 0.7 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe5 = axes('Position',[0.645 0.39 0.35 0.26],'Box','on','NextPlot','add','XTickLabels',[]); 
axe6 = axes('Position',[0.645 0.08 0.35 0.26],'Box','on','NextPlot','add'); 

axes(axe0);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'a.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.55, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(3,1) xy_sta_11(10,1)],[xy_sta_11(3,2) xy_sta_11(10,2)],'-','Color',c1,'LineWidth',1.5) %nl1 nl9
plot([xy_sta_11(3,1) xy_sta_11(2,1)],[xy_sta_11(3,2) xy_sta_11(2,2)],'-','Color',c2,'LineWidth',1.5) %nl1 fl4
plot([xy_sta_11(2,1) xy_sta_11(4,1)],[xy_sta_11(2,2) xy_sta_11(4,2)],'-','Color',c3,'LineWidth',1.5) % nl2 fl04

% trans 
plot([xy_sta_11(4,1) xy_sta_11(10,1)],[xy_sta_11(4,2) xy_sta_11(10,2)],'-','Color',c4,'LineWidth',1.5) % nl2 nl9
plot([xy_sta_11(3,1) xy_sta_11(9,1)],[xy_sta_11(3,2) xy_sta_11(9,2)],'-','Color',c5,'LineWidth',1.5) % nl1 nl8
plot([xy_sta_11(4,1) xy_sta_11(9,1)],[xy_sta_11(4,2) xy_sta_11(9,2)],'-','Color',c6,'LineWidth',1.5) % nl2 nl8

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(2,1)-0.55,xy_sta_11(2,2)-0.50,sites(2),'FontSize',9,'Rotation',0);
text(xy_sta_11(3:4,1)-1.65,xy_sta_11(3:4,2)+0.00,sites(3:4),'FontSize',9,'Rotation',0);
text(xy_sta_11(9,1)-0.35,xy_sta_11(9,2)-0.5,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_11(10,1)-0.65,xy_sta_11(10,2)-0.5,sites(10),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;


axes(axe01);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'b.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.55, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(2,1) xy_sta_11(8,1)],[xy_sta_11(2,2) xy_sta_11(8,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_11(2,1) xy_sta_11(6,1)],[xy_sta_11(2,2) xy_sta_11(6,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_11(6,1) xy_sta_11(10,1)],[xy_sta_11(6,2) xy_sta_11(10,2)],'-','Color',c3,'LineWidth',1.5)
plot([xy_sta_11(7,1) xy_sta_11(9,1)],[xy_sta_11(7,2) xy_sta_11(9,2)],'-','Color',c4,'LineWidth',1.5)

% trans
plot([xy_sta_11(6,1) xy_sta_11(9,1)],[xy_sta_11(6,2) xy_sta_11(9,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_11(7,1) xy_sta_11(10,1)],[xy_sta_11(7,2) xy_sta_11(10,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(2,1)-0.55,xy_sta_11(2,2)-0.50,sites(2),'FontSize',9,'Rotation',0);
text(xy_sta_11(8,1)-0.75,xy_sta_11(8,2)-0.5,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_11(6:7,1)-0.65,xy_sta_11(6:7,2)+0.55,sites(6:7),'FontSize',9,'Rotation',0);
text(xy_sta_11(9,1)-0.35,xy_sta_11(9,2)-0.5,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_11(10,1)-0.65,xy_sta_11(10,2)-0.5,sites(10),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;


axes(axe02);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'c.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
% lakes
scatter(xy_lake(:,1),xy_lake(:,2),4,'s','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_lil_lake(:,1),xy_lil_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
scatter(xy_nnl_lake(:,1),xy_nnl_lake(:,2),4,'bs','filled','MarkerEdgeColor','none','MarkerFaceColor',[0.6 0.6 0.6])
% hydrofractures
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'k','LineWidth',1.5)
plot(xy_SNL_crack(:,1),xy_SNL_crack(:,2),'k-','LineWidth',1.5)
plot(xy_NNL_crack(:,1),xy_NNL_crack(:,2),'k-','LineWidth',1.5)
% lake names
text(1.25, 4.5, 'L1C','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(3.00, 1.55, 'L1A','FontSize',11,'FontName','Avenir','FontAngle','italic')
text(-1.75, -3.5, 'L1B','FontSize',11,'FontName','Avenir','FontAngle','italic')

% lon 
plot([xy_sta_11(7,1) xy_sta_11(11,1)],[xy_sta_11(7,2) xy_sta_11(11,2)],'-','Color',c1,'LineWidth',1.5) % nl6 nl10
plot([xy_sta_11(6,1) xy_sta_11(12,1)],[xy_sta_11(6,2) xy_sta_11(12,2)],'-','Color',c2,'LineWidth',1.5) % nl4 nl11

% trans
plot([xy_sta_11(6,1) xy_sta_11(14,1)],[xy_sta_11(6,2) xy_sta_11(14,2)],'-','Color',c3,'LineWidth',1.5) % nl4 nl13

plot(xy_sta_11(1:14,1),xy_sta_11(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','none');
text(xy_sta_11(6:7,1)-0.65,xy_sta_11(6:7,2)+0.55,sites(6:7),'FontSize',9,'Rotation',0);
text(xy_sta_11(11,1)-1.5,xy_sta_11(11,2)-0.45,sites(11),'FontSize',9,'Rotation',0);
text(xy_sta_11(12,1)+0.3,xy_sta_11(12,2)-0.25,sites(12),'FontSize',9,'Rotation',0);
text(xy_sta_11(14,1)+0.3,xy_sta_11(14,2)+0.15,sites(14),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',8,'FontName','Avenir')

ylabel('y [ km ]'); xlabel('x [ km ]');
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','in','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;
xtickangle(0)

axes(axe1) % L1A and L1C
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(4,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(5,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2) 
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(6,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3) 

plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);

legend('NL01-NL09','NL01-FL04','NL02-FL04'); legend boxoff

ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, '{\bf d. L1A-L1C} (Compare with Fig. S2d.)','FontSize',12); 

text(168.81,0.06,'Start of Precursor','Rotation',90,'FontSize',10);
text(169.17,0.06,'HF Initiation','Rotation',90,'FontSize',10);
text(169.35,0.15,'Max HF','Rotation',90,'FontSize',10);

text(168.85-0.07,-0.18,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.21-0.07,-0.18,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.32-0.07,-0.18,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(170.32-0.08,-0.18,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.8],'FontName','Avenir');  grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');

axes(axe2) % L1A
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3) %nl4 nl9
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(2,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4) %nl6 NL8
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(9,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1) % NL07 Fl04
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2) % NL04 Fl04
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL06-NL08','NL04-NL09','NL04-FL04','NL07-FL04'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);

text(168.33, 0.44, '{\bf e. L1A} (Same as Fig. S2e.)','FontSize',12); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.4],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');


axes(axe3)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(2,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd = legend('NL06-NL10','NL04-NL11','Location','NorthEast'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, '{\bf f. L1A-L1B} (Compare with Fig. S2f)','FontSize',12)
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.6:0.1:0.6],'FontName','Avenir');
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');
xlabel(' Day of Year, 2011  [ UTC ]','FontSize',12,'FontName','Avenir');  grid on;


axes(axe4) % L1A and L1C
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(4,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(4,6,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5) 
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,6,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6) 

plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);

legend('NL01-NL09','NL02-NL09','NL01-NL08','NL02-NL08'); legend boxoff

text(168.81,0.06,'Start of Precursor','Rotation',90,'FontSize',10);
text(169.17,0.11,'HF Initiation','Rotation',90,'FontSize',10);
text(169.35,0.11,'Max HF','Rotation',90,'FontSize',10);
%text(170.27,-0.18,'+1 day','Rotation',90,'FontSize',10);

text(168.85-0.07,-0.18,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.21-0.07,-0.18,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(169.32-0.07,-0.18,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(170.32-0.08,-0.18,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

ylim([-0.22 0.4]); 
xlim([168.32 170.32]);

text(168.33, 0.44, '{\bf g. L1A-L1C} (Compare with Fig. S2g.)','FontSize',12)
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.6],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');


axes(axe5) % L1A 
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5) %nl4 nl8
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(2,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6) %nl6 nl9
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3) % nl4 nl9
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(2,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4) %nl6 NL8
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5) %nl4 nl8
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]); 
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL04-NL08','NL06-NL09','NL04-NL09','NL06-NL08'); legend boxoff
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);

text(168.33, 0.44, '{\bf h. L1A} (Same as Fig. S2h.)','FontSize',12) 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.1:0.6],'FontName','Avenir'); grid on;
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');


axes(axe6) % L1A and L1C
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(2,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,12,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([168.85 168.85], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.21 169.21], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([169.32 169.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([170.32 170.32], [-100 100], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd = legend('NL06-NL10','NL04-NL11','NL04-NL13','Location','NorthEast'); legend boxoff
set(lgd,'Position',[0.89 0.27 0.1 0.05]);
ylim([-0.22 0.4]); 
xlim([168.32 170.32]);
text(168.33, 0.44, '{\bf i. L1A-L1B}  (Compare with Fig. S2i.)','FontSize',12)
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.6:0.1:0.6],'FontName','Avenir');
set(gca,'xtick',[168.5,169,169.5,170]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16,'FontName','Avenir');
xlabel(' Day of Year, 2011  [ UTC ]','FontSize',12,'FontName','Avenir');  grid on;

%print(gcf,'-dpng','-r500','paperfigS4_strainrate_twolakebasins_2011_20240208.png');