%% Nsurface principle stress -- JGR:ES paper figure for principal stresses 
% for idealized case (example)
% 2023 June 21 LAS -- idealized example figure for revisions (1250 m;
% V=0.01 km^3)
% 2023 July 06 LAS -- include 4km E-W vertical HF opening of 0.5 m 
close all; clear all;

%% load saved displacement, strain, stress, and idealized patches
load surface_disp_1250m_EWvert_surfC_4km.mat
load Gsurface1250m_strain_EWvert_surfC_4km.mat
load sigma_ij_1250m_EWvert_surfC_4km.mat
load('../make_idealized_blisters/patches_idealized.mat')

%% Gsurface (Green function matrix) 
% surface positions as defined in NIF runfile + makegeom -- do not change !!!
xy_surf = Gsurface.xy_surf; xx_vec = xy_surf(:,1)'; yy_vec = xy_surf(:,2)'; % [ km ]
xx = Gsurface.xx; yy = Gsurface.yy; Nsurface = Gsurface.Nsurface; % [ km ] 121x121

%% calculate principle stresses
for i = 1:5 
     for j = 1:14641         
% bed open and slip 
        princ_sigma1_1250m(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx(i,j)+sigma_ij_1250m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx(i,j)-sigma_ij_1250m_SLder.sigma_yy(i,j))).^2) + ...
            (sigma_ij_1250m_SLder.sigma_xy(i,j).^2)));            
% bed open 
        princ_sigma1_1250m_open(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_open(i,j)+sigma_ij_1250m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_open(i,j)-sigma_ij_1250m_SLder.sigma_yy_open(i,j))).^2) + ...
            (sigma_ij_1250m_SLder.sigma_xy_open(i,j).^2)));   
% bed slip 
        princ_sigma1_1250m_slip(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_slip(i,j)+sigma_ij_1250m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_slip(i,j)-sigma_ij_1250m_SLder.sigma_yy_slip(i,j))).^2) + ...
            (sigma_ij_1250m_SLder.sigma_xy_slip(i,j).^2)));   
% vertical crack open 
        princ_sigma1_1250m_vert(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_vert(i,j)+sigma_ij_1250m_SLder.sigma_yy_vert(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_vert(i,j)-sigma_ij_1250m_SLder.sigma_yy_vert(i,j))).^2) + ...
            (sigma_ij_1250m_SLder.sigma_xy_vert(i,j).^2)));   
% bed open, bed slip, vertical crack open 
        princ_sigma1_1250m_vertopenslip(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_vertopenslip(i,j)+sigma_ij_1250m_SLder.sigma_yy_vertopenslip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_vertopenslip(i,j)-sigma_ij_1250m_SLder.sigma_yy_vertopenslip(i,j))).^2) + ...
            (sigma_ij_1250m_SLder.sigma_xy_vertopenslip(i,j).^2)));  
     end
end

%% Figure: H = 1250 m; V = 0.1 km^3
v=[-1500:5:1500]; % sigma contours [ kPa ]
v200 = [200 200];  % sigma threshold [ kPa ]
v_disp = [-3:0.025:3]; % surface displacement contours [ m ]
m = 10; % font size
load BWR.mat % colormap
[cbar_def] = [-3 3]; % [ m ]

fig1=figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[1 1 200 120]./5.5);
tt=0.17; ss = 20/12; scale=0.65;

axe1 = axes('Position',[0.04 0.69 tt tt*ss],'Box','on','xticklabel',[]); 
axe2 = axes('Position',[0.22 0.69 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);
axe3 = axes('Position',[0.40 0.69 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);

axe12 = axes('Position',[0.04 0.3725 tt tt*ss],'Box','on','xticklabel',[]); 
axe22 = axes('Position',[0.22 0.3725 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);
axe32 = axes('Position',[0.40 0.3725 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 
axe42 = axes('Position',[0.58 0.3725 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 
axe52 = axes('Position',[0.76 0.3725 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 

axe13 = axes('Position',[0.04 0.055 tt tt*ss],'Box','on'); 
axe23 = axes('Position',[0.22 0.055 tt tt*ss],'Box','on','yticklabel',[]);
axe33 = axes('Position',[0.40 0.055 tt tt*ss],'Box','on','yticklabel',[]); 
axe43 = axes('Position',[0.58 0.055 tt tt*ss],'Box','on','yticklabel',[]); 
axe53 = axes('Position',[0.76 0.055 tt tt*ss],'Box','on','yticklabel',[]); 

axes(axe1)
text(-16.25,13.5,'a','FontWeight','bold','FontSize',m+1);
hold on;
% Vertical Opening (idealize
z=0;
Nx=40; Ny=40; width = Gsurface.patchesB(1,1);
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
    patch(x,y,z,'EdgeColor',[0.6 0.6 0.6]); hold on;
  
  end
 end
% when plotting HF crack, x- of first supfault patch is 1 patch length up-strike 
% from Gsurface.patchesC(1,6). the final subfault patch extends along-strike 
% for one additional patch length (patch length along-strike =
% Gsurface.patchesC(:,1)). ** check visually **
plot([Gsurface.patchesC(1,6)-(Gsurface.patchesC(1,1)) ...
    Gsurface.patchesC(end,6)+(Gsurface.patchesC(1,1))],...
    [Gsurface.patchesC(1,7) Gsurface.patchesC(end,7)],...
    '-k','LineWidth',3);
ylabel('y  [ km ]'); 
text(-11, 11, 'Crack Opening = 0.5 m','FontName','Avenir','FontSize',m-1)
title('Vertical Crack Opening');
caxis([cbar_def(1) cbar_def(2)]);
xlim([-12 12]); ylim([-12 12]);
colormap(axe1,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on

axes(axe2)
text(-12,13.5,'b','FontWeight','bold','FontSize',m+1);
hold on;
% Basal Opening (idealized)
z=patches_idealized.opening(:,5);
Nx=40; Ny=40; width = Gsurface.patchesB(1,1);
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
text(-11, 11, 'V = 0.1 km^{3}','FontName','Avenir','FontSize',m-1)
title('Basal Cavity Opening');
caxis([cbar_def(1) cbar_def(2)]);
xlim([-12 12]); ylim([-12 12]);
colormap(axe2,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on

axes(axe3)
text(-12,13.5,'c','FontWeight','bold','FontSize',m+1);
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

caxis([cbar_def(1) cbar_def(2)]);
text(-11, 11, 'V = 0.1 km^{3}','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
colormap(axe3,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('Basal Slip'); %xlabel('x  [ km ]');

cc2 = colorbar('EastOutside'); 
set(cc2,'Position',[0.58 0.69 .004 (tt*ss)],'xtick',[-3 -2 -1 0 1 2 3],'tickdir','in');
text(17.5,-9.5, 'Imposed Opening or Slip [ m ]','Rotation',90,'FontSize',m,'FontName','Avenir');


% surface displacements

[xx_init,yy_init] = meshgrid(-15:0.25:15); % initial points
step = 1; % [ km ] increase if coarser vector field is desired
[xx_plot,yy_plot]=meshgrid(-12:step:12); % interpolated, 1-km lower-res grid

axes(axe22) % U_{bedopen}
text(-12,13.5,'e','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U_open(5,:))),121,121),v_disp); 
set(h3,'LineColor','none');
FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E_open(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N_open(5,:),121,121),xx_plot,yy_plot);
% scale 150% of max length of max bedopen vector
quiver_scale = 1.50*(max(step/max(max(FE)), step/max(max(FN)))); 
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-9.00, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('u_{bedopen}');
caxis([cbar_def(1) cbar_def(2)]); colormap(axe22,BWR); 
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on 


axes(axe12) % U_{HF}
text(-16.25,13.5,'d','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U_vert(5,:))),121,121),v_disp); 
set(h3,'LineColor','none');

FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E_vert(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N_vert(5,:),121,121),xx_plot,yy_plot);
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-9.00, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');
text(-11, -10.8, 'H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('u_{HF}'); ylabel('y  [ km ]');
caxis([cbar_def(1) cbar_def(2)]); colormap(axe12,BWR); 
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on 


axes(axe32) % U_{bedslip}
text(-12,13.5,'f','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U_slip(5,:))),121,121),v_disp); 
set(h3,'LineColor','none'); 

FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E_slip(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N_slip(5,:),121,121),xx_plot,yy_plot);
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-9.00, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');
caxis([cbar_def(1) cbar_def(2)]); colormap(axe32,BWR); 
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on
title('u_{bedslip}'); 


axes(axe42) % U_{HF+bedopen+bedslip}
text(-12,13.5,'g','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U(5,:))),121,121),v_disp); 
set(h3,'LineColor','none'); 
FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N(5,:),121,121),xx_plot,yy_plot);
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-9.00, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');
caxis([cbar_def(1) cbar_def(2)]); colormap(axe42,BWR); 
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on
title('u_{bedopen+bedslip}'); 


axes(axe52) % U_{bedopen+bedslip}
text(-12,13.5,'h','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((surface_disp_1250m.Disp_U_vertopenslip(5,:))),121,121),v_disp); 
set(h3,'LineColor','none'); 

FE = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_E_vertopenslip(5,:),121,121),xx_plot,yy_plot);
FN = interp2(xx_init,yy_init,reshape(surface_disp_1250m.Disp_N_vertopenslip(5,:),121,121),xx_plot,yy_plot);
quiver(xx_plot,yy_plot,quiver_scale.*FE,quiver_scale.*FN,'k','autoscale','off','LineWidth',0.6,'Color',BWR(10,:));
quiver(-9.00, 11.1, quiver_scale.*1, quiver_scale.*0, 'k','autoscale','off','LineWidth',1.0,...
    'Color',BWR(10,:),'MaxHeadSize',0.8);
text(-9, 10, '1 m','FontSize',m-1,'Color',BWR(10,:),'FontName','Avenir');
caxis([cbar_def(1) cbar_def(2)]); colormap(axe52,BWR); 
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on
title('u_{HF+bedopen+bedslip}'); 

cc2 = colorbar('EastOutside'); 
set(cc2,'Position',[0.94 0.3725 .004 (tt*ss)],'xtick',[-3 -2 -1 0 1 2 3],'tickdir','in');
text(17.4,-8.8, 'Vertical Displacement [ m ]','Rotation',90,'FontSize',m,'FontName','Avenir');


% \sigma_{1} row
cbar_sigma = [-1e3 1e3]; % [ kPa ]

axes(axe13) % \sigma_{1,HF} 
text(-16.25,13.5,'i','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m_vert(5,:))./1e3),121,121),v); set(h3,'LineColor','none'); 
[C3,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_vert(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);
caxis([cbar_sigma(1) cbar_sigma(2)]); colormap(axe13,BWR); 
text(-11, -10.8, 'H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on
title('\sigma_{1,HF}'); xlabel('x  [ km ]'); ylabel('y  [ km ]'); 


axes(axe23) % \sigma_{1,bedopen}
text(-12,13.5,'j','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m_open(5,:))./1e3),121,121),v); set(h3,'LineColor','none'); 
[C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('\sigma_{1,bedopen}'); xlabel('x  [ km ]');
caxis([cbar_sigma(1) cbar_sigma(2)]);  colormap(axe23,BWR); 
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on 


axes(axe33) % \sigma_{1,bedslip}
text(-12,13.5,'k','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m_slip(5,:))./1e3),121,121),v);  set(h3,'LineColor','none'); 
[C2,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);
caxis([cbar_sigma(1) cbar_sigma(2)]);  colormap(axe33,BWR); 
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('\sigma_{1,bedslip}'); xlabel('x  [ km ]');
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on


axes(axe43) % \sigma_{1,bedopen+bedslip}
text(-12,13.5,'l','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m(5,:))./1e3),121,121),v); set(h3,'LineColor','none'); 
[C3,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);
caxis([cbar_sigma(1) cbar_sigma(2)]);  colormap(axe43,BWR); 
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('\sigma_{1,bedopen+bedslip}'); xlabel('x  [ km ]');
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top'); grid on


axes(axe53) % \sigma_{1,HF+bedopen+bedslip}
text(-12,13.5,'m','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1250m_vertopenslip(5,:))./1e3),121,121),v); set(h3,'LineColor','none'); 
[C3,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_vertopenslip(5,:))./1e3),121,121),[v200 v200],'Color',BWR(10,:),'LineWidth',1.1);
caxis([cbar_sigma(1) cbar_sigma(2)]);  colormap(axe53,BWR); 
text(-11, -10.8, 'V = 0.1 km^{3}    H = 1250 m   \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-12 12]); ylim([-12 12]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('\sigma_{1,HF+bedopen+bedslip}'); xlabel('x  [ km ]');


cc43 = colorbar('EastOutside'); 
set(cc43,'Position',[0.94 0.055 .004 (tt*ss)],'xtick',[-1000 -500 0 500 1000],'tickdir','in');
text(17.8,-3.0, '\sigma_{1} [ kPa ]','Rotation',90,'FontSize',m,'FontName','Avenir');

% %% print figure
%  figurename=sprintf('paperfig4_idealized_example_4kmHF_surfC_20230724.png');
%  print(gcf,'-dpng','-r600',figurename);
