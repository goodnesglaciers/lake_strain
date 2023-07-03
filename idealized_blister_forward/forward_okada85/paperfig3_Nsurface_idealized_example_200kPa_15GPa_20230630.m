%% Nsurface principle stress -- JGR:ES paper figure for principal stresses 
% for idealized case (example)
% 2023 June 21 LAS -- idealized example figure for revisions (1250 m;
% V=0.01 km^3)
close all; clear all;

%% load saved displacement, strain, stress, and idealized patches
load surface_disp_1250m.mat
load Gsurface1250m_strain.mat
load sigma_ij_1250m.mat
load('../make_idealized_blisters/patches_idealized.mat')

%% Gsurface (Green function matrix) 
% surface positions as defined in NIF runfile + makegeom -- do not change !!!
xy_surf = Gsurface.xy_surf; xx_vec = xy_surf(:,1)'; yy_vec = xy_surf(:,2)'; % [ km ]
xx = Gsurface.xx; yy = Gsurface.yy; Nsurface = Gsurface.Nsurface; % [ km ] 121x121

%% calculate principle stresses
for i = 1:5
     for j = 1:14641         
% open and slip 
        princ_sigma1_1250m(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx(i,j)+sigma_ij_1250m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx(i,j)-sigma_ij_1250m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_1250m_SLder.sigma_xy(i,j).^2)));            
% open 
        princ_sigma1_1250m_open(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_open(i,j)+sigma_ij_1250m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_open(i,j)-sigma_ij_1250m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_1250m_SLder.sigma_xy_open(i,j).^2)));   
% slip 
        princ_sigma1_1250m_slip(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_slip(i,j)+sigma_ij_1250m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_slip(i,j)-sigma_ij_1250m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_1250m_SLder.sigma_xy_slip(i,j).^2)));   
     end
end

%% Figure: H = 1250 m; V = 0.1 km^3
v=[-1500:5:1500]; % sigma contours [ kPa ]
v200 = [200 200];  % sigma threshold [ kPa ]
v_disp = [-3:0.025:3]; % surface displacement contours [ m ]
m = 10; % font size
load BWR.mat % colormap

fig1=figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[1 1 150 150]./5.5);
tt=0.21; ss = 28/30; scale=0.65;
tt2=tt.*1.00;

axe1 = axes('Position',[0.07 0.78 tt tt*ss],'Box','on','xticklabel',[]); 
axe2 = axes('Position',[0.30 0.78 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);

axe12 = axes('Position',[0.07 0.545 tt tt*ss],'Box','on','xticklabel',[]); 
axe22 = axes('Position',[0.30 0.545 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);
axe32 = axes('Position',[0.53 0.545 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 

axe13 = axes('Position',[0.07 0.31 tt tt*ss],'Box','on'); 
axe23 = axes('Position',[0.30 0.31 tt tt*ss],'Box','on','yticklabel',[]);
axe33 = axes('Position',[0.53 0.31 tt tt*ss],'Box','on','yticklabel',[]); 

axes(axe1)
text(-17.5,14,'a','FontWeight','bold','FontSize',m+1);
hold on;
% Basal Opening (idealized)
z=patches_idealized.opening(:,5);
Nx=40; Ny=40; width = 0.5;
 isf=0;
 for i=1:Nx
  for j=1:Ny
    isf=isf+1;
    
    x1 = Gsurface.patchesB(isf,6);
    x2 = x1 + width; % patch width [km]
    y1 = Gsurface.patchesB(isf,7)-(Gsurface.patchesB(isf,2)/2);
    y2 = y1 + width; % patch width [km]
    x=[x1, x2, x2, x1];
    y=[y1, y1, y2, y2];
    patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6]); hold on;
  
  end
 end
ylabel('y  [ km ]'); 
text(-11, 11, 'V=0.1 km^{3}','FontName','Avenir','FontSize',m-1)
title('Basal Cavity Opening');
caxis([-3 3]); % m
xlim([-12 12]); ylim([-12 12]);
colormap(axe1,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on

%%case 2 -- slip
axes(axe2)
text(-13.5,14,'b','FontWeight','bold','FontSize',m+1);
hold on;
% Basal Slip (idealized)
z=-0.5.*patches_idealized.slip(:,5);
Nx=40; Ny=40; width = 0.5;
 isf=0;
 for i=1:Nx
  for j=1:Ny
    isf=isf+1;
    
    x1 = Gsurface.patchesB(isf,6);
    x2 = x1 + width; % patch width [km]
    y1 = Gsurface.patchesB(isf,7)-(Gsurface.patchesB(isf,2)/2);
    y2 = y1 + width; % patch width [km]
    x=[x1, x2, x2, x1];
    y=[y1, y1, y2, y2];
    patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6]); hold on;
    
  end
 end

caxis([-3 3]); % m
text(-11, 11, 'V=0.1 km^{3}','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
colormap(axe2,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('Basal Slip'); %xlabel('x  [ km ]');

cc2 = colorbar('EastOutside'); 
set(cc2,'Position',[0.52 0.78 .005 (tt*ss)],'xtick',[-3 -2 -1 0 1 2 3],'tickdir','in');
text(17.0,-11.0, 'Imposed Opening or Slip [ m ]','Rotation',90,'FontSize',m,'FontName','Avenir');

% surface displacements

[xx_init,yy_init] = meshgrid(-15:0.25:15); % initial points
step = 1; 
[xx_plot,yy_plot]=meshgrid(-12:step:12); % interpolated, lower-res grid

axes(axe12) % open only
text(-17.5,14,'c','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U_open(5,:))),121,121),v_disp); 
set(h3,'LineColor','none');

FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E_open(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N_open(5,:),121,121),xx_plot,yy_plot);
quiver_scale = max(step/max(max(FE)), step/max(max(FN))); % scale by maximum length
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-8.75, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');

ylabel('y  [ km ]');
text(-11, -10.8, 'V=0.1 km^{3}   H=1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('U_{bedopen}');
caxis([-3 3]); % m
xlim([-12 12]); ylim([-12 12]);
colormap(axe12,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on 

axes(axe22) % slip only
text(-13.5,14,'d','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U_slip(5,:))),121,121),v_disp); 
set(h3,'LineColor','none'); 

FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E_slip(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N_slip(5,:),121,121),xx_plot,yy_plot);
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-8.75, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');

caxis([-3 3]); % m
text(-11,  -10.8, 'V=0.1 km^{3}   H=1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
colormap(axe22,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('U_{bedslip}'); 

axes(axe32) % open and slip
text(-13.5,14,'e','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U(5,:))),121,121),v_disp); 
set(h3,'LineColor','none'); 

FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N(5,:),121,121),xx_plot,yy_plot);
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-8.75, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');

caxis([-3 3]); % m
text(-11, -10.8, 'V=0.1 km^{3}   H=1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
colormap(axe32,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('U_{bedopen+bedslip}'); 

cc2 = colorbar('EastOutside'); 
set(cc2,'Position',[0.75 0.545 .005 (tt*ss)],'xtick',[-3 -2 -1 0 1 2 3],'tickdir','in');
text(17.6,-3.8, 'Uplift [ m ]','Rotation',90,'FontSize',m,'FontName','Avenir');

% sigma1
axes(axe13) % open
text(-17.5,14,'f','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m_open(5,:))./1e3),121,121),v); set(h3,'LineColor','none'); 
[C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);

ylabel('y  [ km ]'); xlabel('x  [ km ]');
text(-11, -10.8, 'V=0.1 km^{3}   H=1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('\sigma_{1,bedopen}');
caxis([-1000 1000]); colormap(axe13,BWR); 
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on 


axes(axe23) % slip
text(-13.5,14,'g','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m_slip(5,:))./1e3),121,121),v);  set(h3,'LineColor','none'); 
[C2,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);

caxis([-1000 1000]); colormap(axe23,BWR); 
text(-11, -10.8, 'V=0.1 km^{3}   H=1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('\sigma_{1,bedslip}'); xlabel('x  [ km ]');

axes(axe33) % slip and open
text(-13.5,14,'h','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m(5,:))./1e3),121,121),v); set(h3,'LineColor','none'); 
[C3,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);

caxis([-1000 1000]); colormap(axe33,BWR); 
text(-11, -10.8, 'V=0.1 km^{3}   H=1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('\sigma_{1,bedopen+bedslip}'); xlabel('x  [ km ]');

cc43 = colorbar('EastOutside'); 
set(cc43,'Position',[0.75 0.31 .005 (tt*ss)],'xtick',[-1000 -500 0 500 1000],'tickdir','in');
text(17.8,-3.0, '\sigma_{1} [ kPa ]','Rotation',90,'FontSize',m,'FontName','Avenir');

%% print figure
% figurename=sprintf('paperfig3_idealized_example_15GPa_20230630.png');
% print(gcf,'-dpng','-r600',figurename);
