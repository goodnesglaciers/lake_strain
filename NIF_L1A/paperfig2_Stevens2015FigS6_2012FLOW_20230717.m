%% Stevens et al. [2015] Supplementary Info FIGURE 6 --> 2012 displacements
% GPS vs NIF: all stations data in 3 long columns
% 11 June 2021: change for flowline stations in 2012 drainage, remove NIF
% outputs
% 13 Oct 2021: remove 3\sigma outliers for the rough day boundary
% 2023 Junly 17: LAS revisions for reviews
clear all; close all

%% 2012 load displ timeseries
load('MLE2012_47_range7_gc_20km_FLOW_plotting.mat') % saved 2012 NIF displ timeseries
t=ts_enu.epochs-ref_epoch;
Vhat=X0_enu.'*ones(1,Nepochs)+V0_enu'*t';
time_2012 = ts_enu.epochs; % time vector
open_time_2012 = 161.85; % time of HF open, 2012

%% load north lake geographic files and morlighem bed relative to North Lake
origin = [68.72, -49.53];
load apcoords_lle_2012_FLOW.mat
    lats=apcoords_lle_2012(1,:); lons=apcoords_lle_2012(2,:); hs=apcoords_lle_2012(3,:);
    llh=[lats; lons; hs];
    xy_sta_12=llh2localxy(llh,origin);
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
% load bedmap
load('morlighem_for_nevis_100km.mat'); % origin = M1 moulin 
load cmapland2.mat % bed topo colormap
% dock at eel pond 
sky_blue = [111, 169, 228]./255; % SNL
metal = [87, 115, 131]./255; % NNL
oar = [251, 219, 154]./255; % NENL
handle = [161, 37, 49]./255; % NL
dark_oar = [164, 114, 63]./255;

%% begin figure
fig21 = figure('Units','centimeters','Position',1.2.*[2 1 18 22]);
    clf
    
    axe0 = axes('Position',[0.055 0.79 0.8750 0.18],'Box','on','NextPlot','add');
    xlabel('x [ km ]','FontSize',9);
    ylabel('y [ km ]','FontSize',9);
    set(gca,'xaxislocation','top','FontSize',8)
    
    axe1 = axes('Position',[0.08 0.04 0.27 0.72],'Box','on','NextPlot','add');
    ylabel('GPS Station','FontSize',11);
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [-1:1:5]);
    station_names = flipud(ts_enu.stations(1:19));
    set(gca,'YTick',-4:0.5:5,'YTickLabel',...
        vertcat(station_names(1),station_names(18:19),station_names(2:3),station_names(4:17)))
    
    axe2 = axes('Position',[0.37 0.04 0.27 0.72],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [-1:1:5]);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    xlabel('Days since 2012 L1A maximum hydro-fracture opening','FontSize',11,'Fontname','Avenir')
    
    axe3 = axes('Position',[0.66 0.04 0.27 0.72],'Box','on','NextPlot','add');
    ylabel('Displacement  [ m ]','FontSize',11,'Fontname','Avenir');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [-1:1:5]);
    set(gca,'yaxislocation','right','ytick',[-4:0.5:5.5])
 
    TriangleSize = 7;
    offset = 0.5; % displacements vertical offset
    offset_v = 0.5*(-10:1:11); % displacements y-intercept
   
axes(axe0)
hold on;
surf(morlighem_for_nevis_100km.X_km./1e3,morlighem_for_nevis_100km.Y_km./1e3,...
    morlighem_for_nevis_100km.B_km-400,'EdgeColor','None','Facealpha',0.5); % plot surface below station points
view([0 90])
caxis([-800 0]); colormap(cmap)
text(-39, 10.75, 'a','FontWeight','bold','FontSize',10);

t2=colorbar('EastOutside');
set(t2,'YTick',[-400,-200,0,200,400],'TickDirection','out'); 
hold all
set(get(t2,'xlabel'),'String',{'Bed Elevation  [ m a.s.l. ]'},'FontSize',8,'Fontname','Avenir');
set(t2, 'Position', [0.915 0.80 .006 0.16],'YTick',[-800,-600,-400,-200,0],...
    'YTickLabel',[-400,-200,0,200,400]);

hold all

% surface ice sheet contours
surf_contours = [0:100:1600]; % m a.s.l.
[C, h] = contour(morlighem_for_nevis_100km.X_km./1e3, morlighem_for_nevis_100km.Y_km./1e3, ...
    morlighem_for_nevis_100km.S_km, surf_contours,'Color',[0.6 0.6 0.6],'LineWidth',0.8);
clabel(C,h,'FontSize',7,'Color',[0.4 0.4 0.4],'FontName','Avenir','LabelSpacing',300)

% main array
plot3(xy_sta_12(3:16,1),xy_sta_12(3:16,2),300.*ones(14,1),'k^','MarkerSize',TriangleSize-4,'MarkerFaceColor',metal,'markerEdgeColor',metal); hold on
% flowline
plot3(xy_sta_12(1:2,1),xy_sta_12(1:2,2),300.*ones(2,1),'^','MarkerSize',TriangleSize-4,'MarkerFaceColor',handle,'markerEdgeColor',handle); 
plot3(xy_sta_12(17:19,1),xy_sta_12(17:19,2),300.*ones(3,1),'^','MarkerSize',TriangleSize-4,'MarkerFaceColor',handle,'markerEdgeColor',handle); 
% lakes
plot3(xy_lake(:,1),xy_lake(:,2),300.*ones(length(xy_lake(:,1)),1),'o','MarkerSize',TriangleSize-6,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot3(xy_lil_lake(:,1),xy_lil_lake(:,2),300.*ones(length(xy_lil_lake(:,1)),1),'o','MarkerSize',TriangleSize-6,'MarkerFaceColor','b','MarkerEdgeColor','none')
plot3(xy_nnl_lake(:,1),xy_nnl_lake(:,2),300.*ones(length(xy_nnl_lake(:,1)),1),'o','MarkerSize',TriangleSize-6,'MarkerFaceColor','b','MarkerEdgeColor','none')

% site names
text(xy_sta_12(1,1)-3,xy_sta_12(1,2)-0.8,'FL03','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(2,1)+0.5,xy_sta_12(2,2)-0.5,'FL04','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(17,1)-3,xy_sta_12(17,2)-1,'FL01','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(18,1)-3,xy_sta_12(18,2)-1,'FL02','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(19,1)+0.5,xy_sta_12(19,2)-1,'FL06','Fontname','AvenirBold','FontSize',6);

text(xy_sta_12(3,1)+0.5,xy_sta_12(3,2)+1,'NL01','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(4,1)-3,xy_sta_12(4,2)+1,'NL02','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(5,1)+0.5,xy_sta_12(5,2)+0.75,'NL03','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(6,1)-3,xy_sta_12(6,2)+0.8,'NL04','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(8,1)+0.5,xy_sta_12(8,2)+0.5,'NL06','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(16,1)+0.2,xy_sta_12(16,2)+0.75,'BS','Fontname','AvenirBold','FontSize',5);
text(xy_sta_12(7,1)-0.95,xy_sta_12(7,2)+0.8,'05','Fontname','AvenirBold','FontSize',5);

text(xy_sta_12(9,1)-0.75,xy_sta_12(9,2)+0.9,'07','Fontname','AvenirBold','FontSize',5);
text(xy_sta_12(10,1)+0.45,xy_sta_12(10,2)+0.0,'08','Fontname','AvenirBold','FontSize',5);
text(xy_sta_12(11,1)+0.25,xy_sta_12(11,2)+0.6,'09','Fontname','AvenirBold','FontSize',5);

text(xy_sta_12(12,1)-3,xy_sta_12(12,2)-1,'NL10','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(13,1)-1.5,xy_sta_12(13,2)-0.9,'NL11','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(14,1)+0.5,xy_sta_12(14,2)-0.9,'NL12','Fontname','AvenirBold','FontSize',6);
text(xy_sta_12(15,1)-1.5,xy_sta_12(15,2)-1,'NL13','Fontname','AvenirBold','FontSize',6);

axis equal; grid on;
set(gca,'FontName','Avenir','FontSize',7,'xtick',[-40:5:50],'ytick',[-12:6:12])
xlim([-40 50]); ylim([-12 12]);

% begin GPS displacements
for i=3:16
 
 II=unique(ts_enu(i,:).epochmindex);

 axes(axe1);% E
 title('Flowline displacement [ m ]','Fontname','Avenir','FontSize',11);
  ind=(i-1)*3+1;
  
% For 30-sec GPS data, use a moving window of 2 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = -1*(ts_enu(i,:).d(1:3:end)-Vhat(ind,II)-F(1,II))-offset_v(i-2);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',metal,'markerSize',2.5);
  hold on;
  xlim([-2 5]);ylim([-4.25 5.5])
  box on

 axes(axe2) %N
  title('Across-flow displacement [ m ]','Fontname','Avenir','FontSize',11)
  ind=(i-1)*3+2;
% For 30-sec GPS data, use a moving window of 2 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(2:3:end)-Vhat(ind,II)-F(2,II)-offset_v(i-2);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',metal,'markerSize',2.5);  
  hold on; 
  xlim([-2 5]); ylim([-4.25 5.5])
  box on
  
 axes(axe3) %U
  title('Vertical displacement [ m ]','Fontname','Avenir','FontSize',11)
  ind=(i-1)*3+3;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(3:3:end)-Vhat(ind,II)-F(3,II)-offset_v(i-2);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',metal,'markerSize',2.5); 
  hold on;
  xlim([-2 5]); ylim([-4.25 5.5])
  
end

for i=1:2 % FL03 FL04
 
 II=unique(ts_enu(i,:).epochmindex);

 axes(axe1);% E
  ind=(i-1)*3+1;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = -1*(ts_enu(i,:).d(1:3:end)-Vhat(ind,II)-F(1,II))-offset_v(i+16);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
  hold on;
  xlim([-1 4]);ylim([-4.25 5.5])
  box on

 axes(axe2) %N
  ind=(i-1)*3+2;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(2:3:end)-Vhat(ind,II)-F(2,II)-offset_v(i+16);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
hold on; 
  xlim([-1 4]); ylim([-4.25 5.5])
  box on

  
 axes(axe3) %U
  ind=(i-1)*3+3;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(3:3:end)-Vhat(ind,II)-F(3,II)-offset_v(i+16);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
hold on;
  xlim([-1 4]); ylim([-4.25 5.5])
end
  
for i=17:18 % FL01 FL02
 
 II=unique(ts_enu(i,:).epochmindex);

 axes(axe1);% E
  ind=(i-1)*3+1;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = -1*(ts_enu(i,:).d(1:3:end)-Vhat(ind,II)-F(1,II))-offset_v(i-2);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
  hold on;
  xlim([-1 4]);ylim([-4.25 5.5])
  box on

 axes(axe2) %N
  ind=(i-1)*3+2;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(2:3:end)-Vhat(ind,II)-F(2,II)-offset_v(i-2);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
  hold on; 
  xlim([-1 4]); ylim([-4.25 5.5])
  box on
  
 axes(axe3) %U
  ind=(i-1)*3+3;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(3:3:end)-Vhat(ind,II)-F(3,II)-offset_v(i-2);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
  hold on;
  xlim([-1 4]); ylim([-4.25 5.5])

end
    
      
for i=19 % FL06
 
 II=unique(ts_enu(i,:).epochmindex);

 axes(axe1);% E
  ind=(i-1)*3+1;
% For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = -1*(ts_enu(i,:).d(1:3:end)-Vhat(ind,II)-F(1,II))-offset_v(i);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
  hold on;
  xlim([-1 4]);ylim([-4.25 5.5])
  box on

 axes(axe2) %N
  ind=(i-1)*3+2;
    % For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(2:3:end)-Vhat(ind,II)-F(2,II)-offset_v(i);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
    hold on; 
  xlim([-1 4]); ylim([-4.25 5.5])
  box on

 axes(axe3) %U
  ind=(i-1)*3+3;
    % For 30-sec GPS data, use a moving window of 4 days in length:
    window = round((4)/((nanmean(diff(ts_enu(i,:).epochs))))); % window =  number of data points
    plot_data = ts_enu(i,:).d(3:3:end)-Vhat(ind,II)-F(3,II)-offset_v(i);
    [TF] = isoutlier(plot_data,'movmean',window); 
    plot(ts_enu(i,:).epochs(~TF,1)-open_time_2012,plot_data(~TF),'.','Color',handle,'markerSize',2.5); 
  hold on;
  xlim([-1 4]); ylim([-4.25 5.5])

end
    
  axes(axe3)
  grid on
  set(gca,'FontSize',8,'tickdir','in','FontName','Avenir')
  text(-0.8, 5.35, 'd','FontWeight','bold','FontSize',10);
  
  axes(axe2)
  grid on
  set(gca,'FontSize',8,'tickdir','in','FontName','Avenir')
  text(-0.8, 5.35, 'c','FontWeight','bold','FontSize',10);

  axes(axe1)
  grid on
  set(gca,'FontSize',8,'tickdir','in','FontName','Avenir')
  text(-0.8, 5.35, 'b','FontWeight','bold','FontSize',10);
%% print figure
figurename=sprintf('paperfig2_flowline_20230717.png');
print(gcf,'-dpng','-r600',figurename);  