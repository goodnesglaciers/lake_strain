%% strain rate plots for 2012 drainage with Track position outputs archived 
% June 2018: reimagined for proposal, added strain rates for L1A, L1C
% May 2019: Added a little map
% Oct 2020: Added duration of elevated pos strain rates
% Oct 2021: Corrected strain rate error estimates based on Jonny's derivation
% June 2023: Figure updates from reviews
% February 2024: replace yellow, 1-sigma errors from TRACK

%% interGPS station distance, station locations in 2012 
clear all; close all;
load('S_neu_err_nif_2012.mat'); % load TRACK NEU position archive

%% interpolate to get consistent time vector, smooth over time window
interptime = 0.0138;      % 20 minutes [decimal day]
xi12=(160.0:interptime:163.0)'; yi12=xi12; zi12=yi12; 
length_xi12=length(xi12);
length_xi12_minus1=length_xi12-1;

Q = SS; % include interpolated values within SS archive
nfQ = 15; % number of stations

timestep = nanmean(diff(Q(1).time_new2(3000:6000)));  % decimal days
timestep_minutes = timestep*24*60; % minutes

span = 23;  % 60-min window width [ points ]
% Remember: points have a 2.6-min window smooth when loaded in from S_2012_* 
% + a 23-point smooth completed here = 60-min window.
timestep_smooth = span.*timestep_minutes % window width [ minutes ]

%return

%% stations 1–15: smooth with 60-min boxcar, then interpret onto 20-min
% time vector
for i=1:15
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


% dx12(5,6,:) = smooth(dx12(5,6,:),5); %1 2
% dy12(5,6,:) = smooth(dy12(5,6,:),5); %1 2
% dv12(5,6,:) = smooth(dv12(5,6,:),5);
% du12(5,6,:) = smooth(du12(5,6,:),5);
% 
% dx12(1,5,:) = smooth(dx12(1,5,:),3); %1 5
% dy12(1,5,:) = smooth(dy12(1,5,:),3); %1 5
% dv12(1,5,:) = smooth(dv12(1,5,:),3);
% du12(1,5,:) = smooth(du12(1,5,:),3);
% 
% dy12(1,7,:) = smooth(dy12(1,7,:),2); %3 5
% dv12(1,7,:) = smooth(dv12(1,7,:),2);
% 
% dy12(6,7,:) = smooth(dy12(6,7,:),3); %3 5
% dv12(6,7,:) = smooth(dv12(6,7,:),3);

%% change to along- and across flowline
flowline =  atand(nanmean(dx12(1,2,:))./nanmean(dy12(1,2,:)));
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

% calculate strin rates (velocity over length) year^{-1} in direction of
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

%% figure for colors
Fig4 = figure(4); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[5 2 19 14]);
axe1 = axes('Position',[0.075 0.70 0.9 0.25],'Box','on','NextPlot','add','XTickLabels',[]);

axes(axe1)
h1=plot([168.85 168.85], [-100 100],'-'); c1 = get(h1,'Color'); hold on;
h2=plot([169.85 168.85], [-100 100],'-'); c2 = get(h2,'Color');
h3=plot([170.85 168.85], [-100 100],'-'); c3 = get(h3,'Color'); % NO YELLOW!!!
h4=plot([171.85 168.85], [-100 100],'-'); c3 = get(h4,'Color');
h5=plot([172.85 168.85], [-100 100],'-'); c4 = get(h5,'Color');
h6=plot([173.85 168.85], [-100 100],'-'); c5 = get(h6,'Color');
h7=plot([174.85 168.85], [-100 100],'-'); c6 = get(h7,'Color');
h8=plot([175.85 168.85], [-100 100],'-'); c7 = get(h8,'Color');
h9=plot([176.85 168.85], [-100 100],'-'); c8 = get(h9,'Color');

%% load scarp position
scarp=load('scarp_20June2011_Euclid.mat');
scarpx=scarp.scarp_20June2011_Euclid(:,1);
scarpy=scarp.scarp_20June2011_Euclid(:,2);
load('lake.mat');
moulin = lake(226,4:5);
radius=6378137.0; eccen=0.08181919;
lat_true=70; lon_posy=-45;
[moulin_x,moulin_y] = polarstereo_fwd(moulin(2),moulin(1),radius,eccen,lat_true,lon_posy);
[lake_x,lake_y] = polarstereo_fwd(lake(:,5),lake(:,4),radius,eccen,lat_true,lon_posy);

%% load north lake geographic files
origin = [68.72, -49.53];
load out2012.mat
    patchesB = out2012.patches;
    patchesC = out2012.patches_C;
load apcoords_lle_2012
    lats=apcoords_lle_2012(1,:); lons=apcoords_lle_2012(2,:); hs=apcoords_lle_2012(3,:);
    llh=[lats; lons; hs];
    xy_sta_12=llh2localxy(llh,origin);
load lake.mat
    llh_lake = [lake(:,5)'; lake(:,4)'; zeros(337,1)'];
    xy_lake = llh2localxy(llh_lake,origin);
load lil_lake_2013_168.mat
    llh_lil_lake = [lil_lake_2013_168(:,4)'; lil_lake_2013_168(:,3)'; zeros(148,1)'];
    xy_lil_lake = llh2localxy(llh_lil_lake,origin);
NNL=csvread('NNL_20110617.csv',1,0);
    radius=6378137.0;
    eccen=0.08181919;
    lat_true=70;
    lon_posy=-45;
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
    x_crack = patchesC(1:6:end,6); y_crack = patchesC(1:6:end,7);
    strike_NL = atand((y_crack(22,1)-y_crack(8,1))./(x_crack(22,1)-x_crack(8,1))) +180  % +180 for plotting purposes    

%% STRAIN RATES AT LAKES 2012
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

sites = {'FL03','FL04','NL01','NL02','NL03','04','NL05','NL06','NL07','NL08','NL09',...
    'NL10','NL11','NL12','NL13','NLBS'};

axes(axe0);
hold on;
TriangleSize=6;
circle_size=8;
text(-7, 5.8, 'a.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_12(1:16,1),xy_sta_12(1:16,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
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
plot([xy_sta_12(3,1) xy_sta_12(5,1)],[xy_sta_12(3,2) xy_sta_12(5,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_12(4,1) xy_sta_12(5,1)],[xy_sta_12(4,2) xy_sta_12(5,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_12(4,1) xy_sta_12(8,1)],[xy_sta_12(4,2) xy_sta_12(8,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_12(3,1) xy_sta_12(7,1)],[xy_sta_12(3,2) xy_sta_12(7,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_12(3,1) xy_sta_12(16,1)],[xy_sta_12(3,2) xy_sta_12(16,2)],'-','Color',c5,'LineWidth',1.5)

plot(xy_sta_12(1:14,1),xy_sta_12(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_12(3:4,1)-1.80,xy_sta_12(3:4,2)-0.10,sites(3:4),'FontSize',9,'Rotation',0);
text(xy_sta_12(5,1)+0.40,xy_sta_12(5,2)-0.10,sites(5),'FontSize',9,'Rotation',0);
text(xy_sta_12(7,1)-0.80,xy_sta_12(7,2)-0.55,sites(7),'FontSize',9,'Rotation',0);
text(xy_sta_12(8,1)-0.20,xy_sta_12(8,2)-0.55,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_12(16,1)-0.20,xy_sta_12(16,2)-0.55,sites(16),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',9,'FontName','Avenir')

ylabel('y [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','out','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe01);
hold on;
text(-7, 5.8, 'b.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_12(1:16,1),xy_sta_12(1:16,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
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
plot([xy_sta_12(7,1) xy_sta_12(11,1)],[xy_sta_12(7,2) xy_sta_12(11,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_12(8,1) xy_sta_12(10,1)],[xy_sta_12(8,2) xy_sta_12(10,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_12(2,1) xy_sta_12(9,1)],[xy_sta_12(2,2) xy_sta_12(9,2)],'-','Color',c3,'LineWidth',1.5)
plot([xy_sta_12(2,1) xy_sta_12(7,1)],[xy_sta_12(2,2) xy_sta_12(7,2)],'-','Color',c4,'LineWidth',1.5)

% trans
plot([xy_sta_12(7,1) xy_sta_12(10,1)],[xy_sta_12(7,2) xy_sta_12(10,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_12(8,1) xy_sta_12(11,1)],[xy_sta_12(8,2) xy_sta_12(11,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_12(1:14,1),xy_sta_12(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_12(2,1)-0.20,xy_sta_12(2,2)+0.55,sites(2),'FontSize',9,'Rotation',0);
text(xy_sta_12(7,1)-0.80,xy_sta_12(7,2)+0.55,sites(7),'FontSize',9,'Rotation',0);
text(xy_sta_12(8,1)-0.20,xy_sta_12(8,2)+0.55,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_12(9,1)-1.00,xy_sta_12(9,2)-0.55,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_12(10,1)-0.60,xy_sta_12(10,2)-0.5,sites(10),'FontSize',9,'Rotation',0);
text(xy_sta_12(11,1)-0.60,xy_sta_12(11,2)-0.55,sites(11),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',9,'FontName','Avenir')

ylabel('y [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','out','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe02);
hold on;
text(-7, 5.8, 'c.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_12(1:16,1),xy_sta_12(1:16,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
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
plot([xy_sta_12(9,1) xy_sta_12(13,1)],[xy_sta_12(9,2) xy_sta_12(13,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_12(11,1) xy_sta_12(12,1)],[xy_sta_12(11,2) xy_sta_12(12,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_12(12,1) xy_sta_12(13,1)],[xy_sta_12(12,2) xy_sta_12(13,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_12(10,1) xy_sta_12(15,1)],[xy_sta_12(10,2) xy_sta_12(15,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_12(10,1) xy_sta_12(13,1)],[xy_sta_12(10,2) xy_sta_12(13,2)],'-','Color',c5,'LineWidth',1.5)

plot(xy_sta_12(1:14,1),xy_sta_12(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_12(9,1)-1.00,xy_sta_12(9,2)-0.55,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_12(10,1)-0.60,xy_sta_12(10,2)+0.5,sites(10),'FontSize',9,'Rotation',0);
text(xy_sta_12(11,1)-0.60,xy_sta_12(11,2)+0.55,sites(11),'FontSize',9,'Rotation',0);
text(xy_sta_12(12,1)-1.00,xy_sta_12(12,2)-0.55,sites(12),'FontSize',9,'Rotation',0);
text(xy_sta_12(13,1)+0.20,xy_sta_12(13,2)-0.5,sites(13),'FontSize',9,'Rotation',0);
text(xy_sta_12(15,1)+0.40,xy_sta_12(15,2)-0.1,sites(15),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',9,'FontName','Avenir')

ylabel('y [ km ]'); xlabel('x [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','out','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;
xtickangle(0)

axes(axe1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(5,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(6,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(6,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL01-NL03','NL02-NL03','NL02-NL06','Location','SouthEast'); legend boxoff
ylim([-0.62 0.3]); 
xlim([160.83 162.83]);
text(160.83, 0.355, 'd. L1C','FontSize',12,'FontWeight','bold'); 
text(161.16,-0.53,'Start of Precursor','Rotation',90,'FontSize',10);
text(161.68,-0.5,'HF Initiation','Rotation',90,'FontSize',10);
text(161.82,-0.5,'Max HF','Rotation',90,'FontSize',10);
text(162.78,0.08,'+1 day','Rotation',90,'FontSize',10);

text(161.17-0.013,0.355,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.69-0.013,0.355,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.83-0.013,0.355,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(162.83-0.016,0.355,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.6:0.2:0.2],'FontName','Avenir');  grid on;
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5]);

axes(axe2)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(3,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(1,9,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(9,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL05-NL09','NL06-NL08','NL05-FL04','NL07-FL04'); legend boxoff
ylim([-0.42 0.82]); 
xlim([160.83 162.83]);
text(160.83, 0.89, 'e. L1A','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.4:0.2:0.8],'FontName','Avenir'); grid on;
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5]);

axes(axe3)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(11,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(4,12,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-1)+(interptime./2), squeeze(lon_yr(12,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL07-NL11','NL09-NL10','NL10-NL11','Location','NorthEast'); legend boxoff
ylim([-0.3 0.61]); 
xlim([160.83 162.83]);
text(160.83, 0.67, 'f. L1B','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.2:0.6],'FontName','Avenir'); grid on;
ylabel('$\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5]);
xlabel(' Day of Year, 2012  [ UTC ]','FontSize',12); 

axes(axe4)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(2,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(6,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
hlegend=legend('NL01-NL05','NL01-NLBS','NL02-NL06');
hlegend.NumColumns = 1; legend boxoff; hlegend.Location = 'SouthEast';
ylim([-0.62 0.3]); 
xlim([160.83 162.83]);
text(160.83, 0.355, 'g. L1C','FontSize',12,'FontWeight','bold')
text(161.16,-0.53,'Start of Precursor','Rotation',90,'FontSize',10);
text(161.68,-0.5,'HF Initiation','Rotation',90,'FontSize',10);
text(161.82,0.05,'Max HF','Rotation',90,'FontSize',10);
text(162.78,0.08,'+1 day','Rotation',90,'FontSize',10);

text(161.17-0.013,0.355,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.69-0.013,0.355,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.83-0.013,0.355,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(162.83-0.025,0.355,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 

set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.8:0.2:0.8],'FontName','Avenir'); grid on;
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',14);
set(gca,'xtick',[161.0,161.5,162.0,162.5]);

axes(axe5)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(4,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
hlegend=legend('NL05-NL08','NL06-NL09','NL05-NL09','NL06-NL08'); 
hlegend.NumColumns = 1; legend boxoff; hlegend.Location = 'NorthEast';
ylim([-0.42 0.82]); 
xlim([160.83 162.83]);
text(160.83, 0.89, 'h. L1A','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.4:0.2:0.8],'FontName','Avenir'); grid on;
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5]);

axes(axe6)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(11,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-1)+(interptime./2), squeeze(trans_yr(3,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd = legend('NL07-NL11','NL08-NL11','NL08-NL13','Location','NorthEast'); legend boxoff
set(lgd,'Position',[0.89 0.23 0.1 0.1]);
ylim([-0.3 0.61]); 
xlim([160.83 162.83]);
text(160.83, 0.67, 'i. L1B','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[-0.2:0.2:0.6]);
ylabel('$\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');
xlabel(' Day of Year, 2012  [ UTC ]','FontSize',12);  grid on;
%
%print(gcf,'-dpng','-r500','paperfig5_strainrate_2012_20240208.png');

%% ERRORS in STRAIN RATES 2012
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
plot(xy_sta_12(1:16,1),xy_sta_12(1:16,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
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
plot([xy_sta_12(3,1) xy_sta_12(5,1)],[xy_sta_12(3,2) xy_sta_12(5,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_12(4,1) xy_sta_12(5,1)],[xy_sta_12(4,2) xy_sta_12(5,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_12(4,1) xy_sta_12(8,1)],[xy_sta_12(4,2) xy_sta_12(8,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_12(3,1) xy_sta_12(7,1)],[xy_sta_12(3,2) xy_sta_12(7,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_12(3,1) xy_sta_12(16,1)],[xy_sta_12(3,2) xy_sta_12(16,2)],'-','Color',c5,'LineWidth',1.5)

plot(xy_sta_12(1:14,1),xy_sta_12(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_12(3:4,1)-1.80,xy_sta_12(3:4,2)-0.10,sites(3:4),'FontSize',9,'Rotation',0);
text(xy_sta_12(5,1)+0.40,xy_sta_12(5,2)-0.10,sites(5),'FontSize',9,'Rotation',0);
text(xy_sta_12(7,1)-0.80,xy_sta_12(7,2)-0.55,sites(7),'FontSize',9,'Rotation',0);
text(xy_sta_12(8,1)-0.20,xy_sta_12(8,2)-0.55,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_12(16,1)-0.20,xy_sta_12(16,2)-0.55,sites(16),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',9,'FontName','Avenir')

ylabel('y [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','out','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe01);
hold on;
text(-7, 5.8, 'b.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_12(1:16,1),xy_sta_12(1:16,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
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
plot([xy_sta_12(7,1) xy_sta_12(11,1)],[xy_sta_12(7,2) xy_sta_12(11,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_12(8,1) xy_sta_12(10,1)],[xy_sta_12(8,2) xy_sta_12(10,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_12(2,1) xy_sta_12(9,1)],[xy_sta_12(2,2) xy_sta_12(9,2)],'-','Color',c3,'LineWidth',1.5)
plot([xy_sta_12(2,1) xy_sta_12(7,1)],[xy_sta_12(2,2) xy_sta_12(7,2)],'-','Color',c4,'LineWidth',1.5)

% trans
plot([xy_sta_12(7,1) xy_sta_12(10,1)],[xy_sta_12(7,2) xy_sta_12(10,2)],'-','Color',c5,'LineWidth',1.5)
plot([xy_sta_12(8,1) xy_sta_12(11,1)],[xy_sta_12(8,2) xy_sta_12(11,2)],'-','Color',c6,'LineWidth',1.5)

plot(xy_sta_12(1:14,1),xy_sta_12(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_12(2,1)-0.20,xy_sta_12(2,2)+0.55,sites(2),'FontSize',9,'Rotation',0);
text(xy_sta_12(7,1)-0.80,xy_sta_12(7,2)+0.55,sites(7),'FontSize',9,'Rotation',0);
text(xy_sta_12(8,1)-0.20,xy_sta_12(8,2)+0.55,sites(8),'FontSize',9,'Rotation',0);
text(xy_sta_12(9,1)-1.00,xy_sta_12(9,2)-0.55,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_12(10,1)-0.60,xy_sta_12(10,2)-0.5,sites(10),'FontSize',9,'Rotation',0);
text(xy_sta_12(11,1)-0.60,xy_sta_12(11,2)-0.55,sites(11),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',9,'FontName','Avenir')

ylabel('y [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]);
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','out','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
grid on;

axes(axe02);
hold on;
text(-7, 5.8, 'c.','FontWeight','bold','FontSize',13);
hold on;
plot(xy_sta_12(1:16,1),xy_sta_12(1:16,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
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
plot([xy_sta_12(9,1) xy_sta_12(13,1)],[xy_sta_12(9,2) xy_sta_12(13,2)],'-','Color',c1,'LineWidth',1.5)
plot([xy_sta_12(11,1) xy_sta_12(12,1)],[xy_sta_12(11,2) xy_sta_12(12,2)],'-','Color',c2,'LineWidth',1.5)
plot([xy_sta_12(12,1) xy_sta_12(13,1)],[xy_sta_12(12,2) xy_sta_12(13,2)],'-','Color',c3,'LineWidth',1.5)

% trans
plot([xy_sta_12(10,1) xy_sta_12(15,1)],[xy_sta_12(10,2) xy_sta_12(15,2)],'-','Color',c4,'LineWidth',1.5)
plot([xy_sta_12(10,1) xy_sta_12(13,1)],[xy_sta_12(10,2) xy_sta_12(13,2)],'-','Color',c5,'LineWidth',1.5)

plot(xy_sta_12(1:14,1),xy_sta_12(1:14,2),'k^','MarkerSize',TriangleSize-2,'MarkerFaceColor','k');
text(xy_sta_12(9,1)-1.00,xy_sta_12(9,2)-0.55,sites(9),'FontSize',9,'Rotation',0);
text(xy_sta_12(10,1)-0.60,xy_sta_12(10,2)+0.5,sites(10),'FontSize',9,'Rotation',0);
text(xy_sta_12(11,1)-0.60,xy_sta_12(11,2)+0.55,sites(11),'FontSize',9,'Rotation',0);
text(xy_sta_12(12,1)-1.00,xy_sta_12(12,2)-0.55,sites(12),'FontSize',9,'Rotation',0);
text(xy_sta_12(13,1)+0.20,xy_sta_12(13,2)-0.5,sites(13),'FontSize',9,'Rotation',0);
text(xy_sta_12(15,1)+0.40,xy_sta_12(15,2)-0.1,sites(15),'FontSize',9,'Rotation',0);

x_north = [-4.166666666]; y_north = [3]; u_north = [0.1577]; v_north = [1.5]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.466666666,2.7,'N','FontSize',12)
x_north = [-2.5]; y_north = [-4.166666667]; u_north = [-2]; v_north = [0]; 
quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.5,'LineWidth',1.1)
text(-4.8,-4.7,'Flow Direction','FontSize',9,'FontName','Avenir')

ylabel('y [ km ]'); xlabel('x [ km ]'); 
xlim([-5 5]); ylim([-5.3 5.3]); grid on;
set(gca,'xtick',[-5:1:5],'ytick',[-5:1:5],'tickdir','out','LineWidth',1.1,'FontSize',10,'FontName','Avenir'); 
xtickangle(0)

axes(axe1)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(5,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(6,7,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(6,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL01-NL03','NL02-NL03','NL02-NL06','Location','NorthEast'); legend boxoff
ylim([0 0.00006]); 
xlim([160.83 162.83]);
text(160.87, 5.6*(10^-5), 'd. L1C','FontSize',12,'FontWeight','bold'); 
text(161.17-0.013,0.000065,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.69-0.013,0.000065,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.83-0.013,0.000065,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(162.83-0.015,0.000065,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001]);  grid on;
ylabel('$\delta\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',14);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');

axes(axe2)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(3,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(1,9,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(9,11,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd=legend('NL05-NL09','NL06-NL08','NL05-FL04','NL07-FL04'); legend boxoff
set(lgd,'FontSize',9)
ylim([0 0.00006]); 
xlim([160.83 162.83]);
text(160.87, 5.6*(10^-5), 'e. L1A','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001]); grid on;
ylabel('$\delta\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');

axes(axe3)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(11,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(4,12,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Exx_yr(12,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
legend('NL07-NL11','NL09-NL10','NL10-NL11','Location','NorthEast'); legend boxoff
ylim([0 0.00006]); 
xlim([160.83 162.83]);
text(160.87, 5.6*(10^-5), 'f. L1B','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001]); grid on;
ylabel('$\delta\dot{\epsilon}_{xx}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');
xlabel(' Day of Year, 2012  [ UTC ]','FontSize',12); 

axes(axe4)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(2,5,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(6,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c3)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
hlegend=legend('NL01-NL05','NL01-NLBS','NL02-NL06');
hlegend.NumColumns = 1; legend boxoff; hlegend.Location = 'NorthEast';
ylim([0 0.00006]); 
xlim([160.83 162.83]);
text(160.87, 5.6*(10^-5), 'g. L1C','FontSize',12,'FontWeight','bold')
text(161.17-0.013,0.000065,'\it{t}_{1}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.69-0.013,0.000065,'\it{t}_{2}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
text(161.83-0.013,0.000065,'\it{t}_{3}','FontSize',11,'FontName','Avenir'); 
text(162.83-0.025,0.000065,'\it{t}_{4}','Rotation',0,'FontSize',11,'FontName','Avenir'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001]); grid on;
ylabel('$\delta\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');

axes(axe5)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,3,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(4,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c6)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(1,4,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(3,10,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c2)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
hlegend=legend('NL05-NL08','NL06-NL09','NL05-NL09','NL06-NL08'); 
hlegend.NumColumns = 1; legend boxoff; hlegend.Location = 'NorthEast';
ylim([0 0.00006]); 
xlim([160.83 162.83]);
text(160.87, 5.6*(10^-5), 'h. L1A','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001]); grid on;
ylabel('$\delta\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');

axes(axe6)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(11,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c1)
hold on
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(3,13,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c5)
plot(xi12(4:end-2)+(interptime./2), squeeze(delta_Eyy_yr(3,14,4:end)),'.-','LineWidth',1.1,'MarkerSize',circle_size,'Color',c4)
plot([161.2 161.2], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.72 161.72], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
plot([161.85 161.85], [-10 10], 'k','LineWidth',1,'Color',[0.4 0.4 0.4]);
lgd = legend('NL07-NL11','NL08-NL11','NL08-NL13','Location','NorthEast'); legend boxoff
set(lgd,'Position',[0.89 0.23 0.1 0.1]);
ylim([0 0.00006]); 
xlim([160.83 162.83]);
text(160.87, 5.6*(10^-5), 'i. L1B','FontSize',12,'FontWeight','bold'); 
set(gca,'tickdir','in','LineWidth',1.1,'FontSize',12,'ytick',[0:0.00002:0.0001]);
xlabel(' Day of Year, 2012  [ UTC ]','FontSize',12);  grid on;
ylabel('$\delta\dot{\epsilon}_{yy}\  [\ yr^{-1}\ ]$', 'Interpreter','latex','FontSize',16);
set(gca,'xtick',[161.0,161.5,162.0,162.5],'FontName','Avenir');
%
%print(gcf,'-dpng','-r500','paperfigS3_strainrate_errors_2012_20240208.png');