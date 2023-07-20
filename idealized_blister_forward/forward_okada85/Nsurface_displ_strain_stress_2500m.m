%% Forward model surface displacements and stress over 20 km x 20 km region 
% at the half-space surface centered on the lake 
% Laura A. Stevens Oct 26, 2015 (while in Oxford)
% 25 Feb 2021 LAS
% 4 March 2021 LAS -- now with strains! for 2012!
% 22 March 2021 LAS -- idealized simulations with H=500m
% 07 April 2023 LAS -- S. Larochelle finds issue with strain sign conventions 
% 30 June 2023 LAS -- paper revisions, idealized patch locations fixed
clear all; close all

%% load Green function matrix, idealized opening and slip
load Gsurface2500m_strain.mat % 40 x 40 patch bed
load('../make_idealized_blisters/patches_idealized.mat') % 40 x 40 patch bed
% extract opening and slip along planes
bed_crack_open_m =patches_idealized.opening'; %  [ m ] from CYL derivation
bed_crack_slip_m =-0.5.*patches_idealized.slip'; % 0.5 m down-dip (-x direction = west = down flowline)
Vtot = [0.001, 0.005, 0.01, 0.05, 0.1]; % [ km^3 ] idealized lake volumes

%% calculate surface displacements (forward problem)
% surface positions as defined in NIF runfile + makegeom -- do not change !!!
xy_surf = Gsurface.xy_surf; xx = Gsurface.xx; yy = Gsurface.yy; 
xx_vec = xy_surf(:,1)'; yy_vec = xy_surf(:,2)'; % [ km ]
Nsurface = Gsurface.Nsurface; % length(xy_surf);

for j=1:length(bed_crack_open_m(:,1))
    % bed opening
    B_G3bv = Gsurface.G3bv*bed_crack_open_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;

        B_G3bv_E(j,i) = B_G3bv(ind1); %E
        B_G3bv_N(j,i) = B_G3bv(ind2); %N
        B_G3bv_U(j,i) = B_G3bv(ind3); %U
    end
end

for j=1:length(bed_crack_slip_m(:,1))
    % bed slip
    B_G2bv = Gsurface.G2bv*bed_crack_slip_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;

        B_G2bv_E(j,i) = B_G2bv(ind1); %E
        B_G2bv_N(j,i) = B_G2bv(ind2); %N
        B_G2bv_U(j,i) = B_G2bv(ind3); %U
    end
end

% linear: add all sources of displacement
Disp_E =  B_G2bv_E + B_G3bv_E ;
Disp_N =  B_G2bv_N + B_G3bv_N ;
Disp_U =  B_G2bv_U + B_G3bv_U ;

Disp_E_open =   B_G3bv_E ;
Disp_N_open =   B_G3bv_N ;
Disp_U_open =   B_G3bv_U ;

Disp_E_slip =   B_G2bv_E ;
Disp_N_slip =   B_G2bv_N ;
Disp_U_slip =   B_G2bv_U ;

%% save displacement timeseries
surface_disp_2500m.Disp_E = Disp_E;
surface_disp_2500m.Disp_N = Disp_N;
surface_disp_2500m.Disp_U = Disp_U;
surface_disp_2500m.Disp_E_slip = Disp_E_slip;
surface_disp_2500m.Disp_N_slip = Disp_N_slip;
surface_disp_2500m.Disp_U_slip = Disp_U_slip;
surface_disp_2500m.Disp_E_open = Disp_E_open;
surface_disp_2500m.Disp_N_open = Disp_N_open;
surface_disp_2500m.Disp_U_open = Disp_U_open;
save surface_disp_2500m.mat surface_disp_2500m
%load surface_disp_2500m.mat

%% calculate surface strains (forward problem) from SL
% GsurfaceNUMBER (Green function matrix) --> strains and "tilts"

for j=1:length(bed_crack_open_m(:,1))
    % bed opening (multiply by -1 to get sign convention to positive strain = TENSION)
    B_G3bv = -1.*Gsurface.G3bv_strain*bed_crack_open_m(j,:)'; 
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;
        ind4=Nsurface*3+i; ind5=Nsurface*4+i; ind6=Nsurface*5+i;
        
        B_G3bv_uZE(j,i) = B_G3bv(ind1); %uZE = duz/dx
        B_G3bv_uZN(j,i) = B_G3bv(ind2); %uZN = duz/dy
        B_G3bv_uNN(j,i) = B_G3bv(ind3); %uNN = duy/dy
        B_G3bv_uNE(j,i) = B_G3bv(ind4); %uNE = duy/dx
        B_G3bv_uEN(j,i) = B_G3bv(ind5); %uEN = dux/dy
        B_G3bv_uEE(j,i) = B_G3bv(ind6); %uEE = dux/dx
    end
end

for j=1:length(bed_crack_slip_m(:,1))
    % bed slip (multiply by -1 to get sign convention to positive strain = TENSION)
    B_G2bv = -1.*Gsurface.G2bv_strain*bed_crack_slip_m(j,:)'; 
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;
        ind4=Nsurface*3+i; ind5=Nsurface*4+i; ind6=Nsurface*5+i;
        
        B_G2bv_uZE(j,i) = B_G2bv(ind1);     %uZE = duz/dx
        B_G2bv_uZN(j,i) = B_G2bv(ind2);     %uZN = duz/dy
        B_G2bv_uNN(j,i) = B_G2bv(ind3);     %uNN = duy/dy
        B_G2bv_uNE(j,i) = B_G2bv(ind4);     %uNE = duy/dx
        B_G2bv_uEN(j,i) = B_G2bv(ind5);     %uEN = dux/dy
        B_G2bv_uEE(j,i) = B_G2bv(ind6);     %uEE = dux/dx
    end
end

% linear elastic: add all sources of strain 
NU = 0.3; % Poisson ratio

% open and slip
duz_dx =  (B_G3bv_uZE + B_G2bv_uZE)./1e3; % u in m but dx is in km
duz_dy =  (B_G3bv_uZN + B_G2bv_uZN)./1e3;
dux_dx =  (B_G3bv_uEE + B_G2bv_uEE)./1e3;
duy_dx =  (B_G3bv_uNE + B_G2bv_uNE)./1e3;
duy_dy =  (B_G3bv_uNN + B_G2bv_uNN)./1e3;
dux_dy =  (B_G3bv_uEN + B_G2bv_uEN)./1e3;
dux_dz =  (-1.*duz_dx); % From free surface boundary condition (strain_xz = 0 at the surface)
duy_dz =  (-1.*duz_dy); % From free surface boundary condition (strain_yz = 0 at the surface)
duz_dz =  (-1.*(dux_dx + duy_dy).*(NU/(1-NU))); % From free surface boundary condition (stress_zz = 0 at the surface)

% slip only
duz_dx_slip =   (B_G2bv_uZE)./1e3; % u in m but dx is in km
duz_dy_slip =   (B_G2bv_uZN)./1e3;
dux_dx_slip =   (B_G2bv_uEE)./1e3;
duy_dx_slip =   (B_G2bv_uNE)./1e3;
duy_dy_slip =   (B_G2bv_uNN)./1e3;
dux_dy_slip =   (B_G2bv_uEN)./1e3;
dux_dz_slip =  (-1.*duz_dx_slip); % From free surface boundary condition (strain_xz = 0 at the surface)
duy_dz_slip =  (-1.*duz_dy_slip); % From free surface boundary condition (strain_yz = 0 at the surface)
duz_dz_slip =  (-1.*(dux_dx_slip + duy_dy_slip).*(NU/(1-NU))); % From free surface boundary condition (stress_zz = 0 at the surface)

% open only
duz_dx_open =  (B_G3bv_uZE)./1e3; % u in m but dx is in km
duz_dy_open =  (B_G3bv_uZN)./1e3;
dux_dx_open =  (B_G3bv_uEE)./1e3;
duy_dx_open =  (B_G3bv_uNE)./1e3;
duy_dy_open =  (B_G3bv_uNN)./1e3;
dux_dy_open =  (B_G3bv_uEN)./1e3;
dux_dz_open =  (-1.*duz_dx_open); % From free surface boundary condition (strain_xz = 0 at the surface)
duy_dz_open =  (-1.*duz_dy_open); % from tilts to strains; see okada85.m header
duz_dz_open =  (-1.*(dux_dx_open + duy_dy_open).*(NU/(1-NU))); % From free surface boundary condition (stress_zz = 0 at the surface)

%% (Infinitesimal) strain tensor at the surface (SL)
% eps_ij = 0.5*(dui_dj + duj_di)
% eps is short for epsilon

eps_xx = dux_dx; eps_xy = 0.5*(dux_dy+duy_dx); eps_xz = 0.5*(dux_dz+duz_dx);
eps_yx = 0.5*(duy_dx+dux_dy); eps_yy = duy_dy; eps_yz = 0.5*(duy_dz+duz_dy);
eps_zx = 0.5*(duz_dx+dux_dz); eps_zy = 0.5*(duz_dy+duy_dz); eps_zz = duz_dz;

%% Stress tensor at the surface
% sigma_ij = lambda * delta_ij * eps_kk + 2*mu*eps_ij
% where eps_kk = eps_xx + eps_yy + eps_zz is the volumetric strain

% shear modulus (mu) middle estimate
lambda = 2.25e9; 
mu = 1.5e9;
% 
% % shear modulus (mu) lower estimate
% lambda = 0.48e9; 
% mu = 0.32e9;
% 
% % shear modulus (mu) upper estimate
% lambda = 5.85e9; 
% mu = 3.9e9;

eps_kk = eps_xx + eps_yy + eps_zz; 

sigma_xx = lambda*eps_kk+2*mu*eps_xx; sigma_xy = 2*mu*eps_xy; sigma_xz = 2*mu*eps_xz;
sigma_yx = 2*mu*eps_yx; sigma_yy = lambda*eps_kk+2*mu*eps_yy; sigma_yz = 2*mu*eps_yz;
sigma_zx = 2*mu*eps_zx; sigma_zy = 2*mu*eps_zy; sigma_zz = lambda*eps_kk+2*mu*eps_zz;

%% Slip only 
% (Infinitesimal) strain tensor at the surface (SL)
% eps_ij = 0.5*(dui_dj + duj_di)
% eps is short for epsilon
eps_xx_slip = dux_dx_slip; eps_xy_slip = 0.5*(dux_dy_slip+duy_dx_slip); eps_xz_slip = 0.5*(dux_dz_slip+duz_dx_slip);
eps_yx_slip = 0.5*(duy_dx_slip+dux_dy_slip); eps_yy_slip = duy_dy_slip; eps_yz_slip = 0.5*(duy_dz_slip+duz_dy_slip);
eps_zx_slip = 0.5*(duz_dx_slip+dux_dz_slip); eps_zy_slip = 0.5*(duz_dy_slip+duy_dz_slip); eps_zz_slip = duz_dz_slip;

% Stress tensor at the surface
% sigma_ij = lambda * delta_ij * eps_kk + 2*mu*eps_ij
% where eps_kk = eps_xx + eps_yy + eps_zz is the volumetric strain
eps_kk_slip = eps_xx_slip + eps_yy_slip + eps_zz_slip; 

sigma_xx_slip = lambda*eps_kk_slip+2*mu*eps_xx_slip; sigma_xy_slip = 2*mu*eps_xy_slip; sigma_xz_slip = 2*mu*eps_xz_slip;
sigma_yx_slip = 2*mu*eps_yx_slip; sigma_yy_slip = lambda*eps_kk_slip+2*mu*eps_yy_slip; sigma_yz_slip = 2*mu*eps_yz_slip;
sigma_zx_slip = 2*mu*eps_zx_slip; sigma_zy_slip = 2*mu*eps_zy_slip; sigma_zz_slip = lambda*eps_kk_slip+2*mu*eps_zz_slip;

%% Open only 
% (Infinitesimal) strain tensor at the surface (SL)
% eps_ij = 0.5*(dui_dj + duj_di)
% eps is short for epsilon
eps_xx_open = dux_dx_open; eps_xy_open = 0.5*(dux_dy_open+duy_dx_open); eps_xz_open = 0.5*(dux_dz_open+duz_dx_open);
eps_yx_open = 0.5*(duy_dx_open+dux_dy_open); eps_yy_open = duy_dy_open; eps_yz_open = 0.5*(duy_dz_open+duz_dy_open);
eps_zx_open = 0.5*(duz_dx_open+dux_dz_open); eps_zy_open = 0.5*(duz_dy_open+duy_dz_open); eps_zz_open = duz_dz_open;

% Stress tensor at the surface
% sigma_ij = lambda * delta_ij * eps_kk + 2*mu*eps_ij
% where eps_kk = eps_xx + eps_yy + eps_zz is the volumetric strain
eps_kk_open = eps_xx_open + eps_yy_open + eps_zz_open; 

sigma_xx_open = lambda*eps_kk_open+2*mu*eps_xx_open; sigma_xy_open = 2*mu*eps_xy_open; sigma_xz_open = 2*mu*eps_xz_open;
sigma_yx_open = 2*mu*eps_yx_open; sigma_yy_open = lambda*eps_kk_open+2*mu*eps_yy_open; sigma_yz_open = 2*mu*eps_yz_open;
sigma_zx_open = 2*mu*eps_zx_open; sigma_zy_open = 2*mu*eps_zy_open; sigma_zz_open = lambda*eps_kk_open+2*mu*eps_zz_open;

%% save for plotting 
sigma_ij_2500m_SLder.sigma_xx = sigma_xx;
sigma_ij_2500m_SLder.sigma_xy = sigma_xy;
sigma_ij_2500m_SLder.sigma_xz = sigma_xz;
sigma_ij_2500m_SLder.sigma_yx = sigma_yx;
sigma_ij_2500m_SLder.sigma_yy = sigma_yy;
sigma_ij_2500m_SLder.sigma_yz = sigma_yz;
sigma_ij_2500m_SLder.sigma_zx = sigma_zx;
sigma_ij_2500m_SLder.sigma_zy = sigma_zy;
sigma_ij_2500m_SLder.sigma_zz = sigma_zz;

sigma_ij_2500m_SLder.sigma_xx_slip = sigma_xx_slip;
sigma_ij_2500m_SLder.sigma_xy_slip = sigma_xy_slip;
sigma_ij_2500m_SLder.sigma_xz_slip = sigma_xz_slip;
sigma_ij_2500m_SLder.sigma_yx_slip = sigma_yx_slip;
sigma_ij_2500m_SLder.sigma_yy_slip = sigma_yy_slip;
sigma_ij_2500m_SLder.sigma_yz_slip = sigma_yz_slip;
sigma_ij_2500m_SLder.sigma_zx_slip = sigma_zx_slip;
sigma_ij_2500m_SLder.sigma_zy_slip = sigma_zy_slip;
sigma_ij_2500m_SLder.sigma_zz_slip = sigma_zz_slip;

sigma_ij_2500m_SLder.sigma_xx_open = sigma_xx_open;
sigma_ij_2500m_SLder.sigma_xy_open = sigma_xy_open;
sigma_ij_2500m_SLder.sigma_xz_open = sigma_xz_open;
sigma_ij_2500m_SLder.sigma_yx_open = sigma_yx_open;
sigma_ij_2500m_SLder.sigma_yy_open = sigma_yy_open;
sigma_ij_2500m_SLder.sigma_yz_open = sigma_yz_open;
sigma_ij_2500m_SLder.sigma_zx_open = sigma_zx_open;
sigma_ij_2500m_SLder.sigma_zy_open = sigma_zy_open;
sigma_ij_2500m_SLder.sigma_zz_open = sigma_zz_open;

save sigma_ij_2500m.mat sigma_ij_2500m_SLder

%% pre-figure load datasets
load apcoords_lle
    lats=apcoords_lle(1,:);
    lons=apcoords_lle(2,:);
    hs=apcoords_lle(3,:); % UNITS don't matter, not used
    llh=[lats; lons; hs];
    origin=[68.72, -49.53];
    xy_sta_11=llh2localxy(llh,origin);
    load moulin_2011.mat
    llh_moulin = [moulin_2011(:,4)'; moulin_2011(:,3)'; zeros(37,1)'];
    xy_moulin = llh2localxy(llh_moulin,origin);
    load lake.mat
    llh_lake = [lake(:,5)'; lake(:,4)'; zeros(337,1)'];
    xy_lake = llh2localxy(llh_lake,origin);    
    load BWR.mat
    load out2012.mat
    patchesB = out2012.patches;
    patchesC = out2012.patches_C;
    cmin_vert = -0.15; cmax_vert = 0.15;
    cmin_bup = -1.0; cmax_bup = 1.0;
    cmin_bslip = -0.5; cmax_bslip = 0.5;
    TriangleSize = 7;
    scale_factor = 10;  scale_factor2 = 10;

%% %%%%%% plot forward problem surface displacements
%% bed open + bed slip
for i=1:1:5
fig1 = figure('Units','centimeters','Position',[0.5 0.5 25.5 22]);%,...
                 % 'Name',[titlestring ' 2'],'NumberTitle','off');
    clf
    axe1 = axes('Position',[0.059 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe2 = axes('Position',[0.059 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe3 = axes('Position',[0.059 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe4 = axes('Position',[0.379 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe5 = axes('Position',[0.379 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe6 = axes('Position',[0.379 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
    
    axe7 = axes('Position',[0.699 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe8 = axes('Position',[0.699 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe9 = axes('Position',[0.699 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
 %   
% bed open, bed slip --> surface displacement
axes(axe1)
vv=[-1.5025:0.005:1.5025]; 
[C1,h1]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_E(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-0.5 0.5]);
colormap(axe1,BWR); 
t1=colorbar('EastOutside'); set(t1,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t1,'ylabel'),'String','x-displacement  [ m ]','FontSize',10);
set(t1, 'Position', [.295 0.70 .005 (0.27)-0.02]);

text(-29, 20.2, 'a.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Open and Slip: Displacement, U','FontSize',10)
text(-21,19,sprintf('H = 2500 m   V = %4.3f km^{3}',Vtot(i)),'FontWeight','bold','FontSize',10);
% axis equal

axes(axe2)
[C2,h2]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_N(i,:),121,121),vv);
set(h2,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar;  
caxis([-0.5 0.5]);
colormap(axe2,BWR); 
t2=colorbar('EastOutside'); set(t2,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t2,'ylabel'),'String','y-displacement  [ m ]','FontSize',10);
set(t2, 'Position', [.295 0.3825 .005 (0.27)-0.02]);
text(-29, 20.2, 'b.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

axes(axe3)
vvv=[-1.0025:0.005:1.0025];
[C3,h3]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_U(i,:),121,121),vvv);
set(h3,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-1 1]);
colormap(axe3,BWR); 
t3=colorbar('EastOutside'); set(t3,'YTick',[-1,-0.5,0,0.5,1]); hold all
set(get(t3,'ylabel'),'String','z-displacement  [ m ]','FontSize',10);
set(t3, 'Position', [.295 0.065 .005 (0.27)-0.02]);
% axis equal  

text(-29, 20.2, 'c.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10); 
xlabel(' x [ km ]','FontSize',10);
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
  
% bed open, bed slip --> surface strain
axes(axe4)
vv=(10^(-4)).*[-1.5025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_xx(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe4,BWR); 
t4=colorbar('EastOutside'); set(t4,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t4,'ylabel'),'String','\epsilon_{xx}  [ m/m ]','FontSize',10);
set(t4, 'Position', [.615 0.70 .005 (0.27)-0.02]);

% x_north = [-8.8]; y_north = [8]; u_north = [0]; v_north = [1.5]; 
% quiver(x_north, y_north, u_north, v_north, 'k','AutoScale','off','MaxHeadSize',1.2,'LineWidth',1.2)
% text(-9.366666666,7.7,'N','FontSize',11,'FontWeight','bold')

text(-24, 20.2, 'd.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Open and Slip: Strain, \epsilon','FontSize',10)
% axis equal

axes(axe5)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_yy(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe5,BWR); 
t5=colorbar('EastOutside'); set(t5,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t5,'ylabel'),'String','\epsilon_{yy}  [ m/m ]','FontSize',10);
set(t5, 'Position', [.615 0.3825 .005 (0.27)-0.02]);

text(-24, 20.2, 'e.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
% axis equal

axes(axe6)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_zz(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe6,BWR); 
t6=colorbar('EastOutside'); set(t6,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t6,'ylabel'),'String','\epsilon_{zz}  [ m/m ]','FontSize',10);
set(t6, 'Position', [.615 0.065 .005 (0.27)-0.02]);
text(-24, 20.2, 'f.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlabel(' x [ km ]','FontSize',10);
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
% axis equal

% STRESSES sigma_xx
axes(axe7) 
vvKPA=[-8005:10:8005]; v150=150;
[C4,h4]=contourf(xx, yy, reshape(sigma_xx(i,:),121,121)./1e3,vvKPA);
set(h4,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_xx(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-600 600]);
colormap(axe7,BWR); 
t7=colorbar('EastOutside'); set(t7,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t7,'ylabel'),'String','\sigma_{xx}  [ kPa ]','FontSize',10);
set(t7,'Position', [.935 0.70 .005 (0.27)-0.02]);
text(-24, 20.2, 'g.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('Open and Slip: Stress, \sigma','FontSize',10)


% sigma_yy
axes(axe8)
[C5,h5]=contourf(xx, yy, reshape(sigma_yy(i,:),121,121)./1e3,vvKPA);
set(h5,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_yy(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-600 600]);
colormap(axe8,BWR); 
t8=colorbar('EastOutside'); set(t8,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t8,'ylabel'),'String','\sigma_{yy}  [ kPa ]','FontSize',10);
set(t8,'Position', [.935 0.3825 .005 (0.27)-0.02]);
text(-24, 20.2, 'h.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% sigma_zz
axes(axe9)
[C5,h6]=contourf(xx, yy, reshape(sigma_zz(i,:),121,121)./1e3,vvKPA);
set(h6,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_zz(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-600 600]);
colormap(axe9,BWR); 
t9=colorbar('EastOutside'); set(t9,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t9,'ylabel'),'String','\sigma_{zz}  [ kPa ]','FontSize',10);
set(t9, 'Position', [.935 0.065 .005 (0.27)-0.02]);
text(-24, 20.2, 'i.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
xlabel(' x [ km ]','FontSize',10);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% print movie frame
rez=300; %resolution (dpi) of final graphic
figpos=getpixelposition(fig1); %dont need to change anything here
resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
set(fig1,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
path='Forward_2500m_check_plots'; %the folder where you want to put the file
name=sprintf('%04d_2500m_bedopen_bedslip',i); %what you want the file to be called
print(fig1,fullfile(path,name),'-dpng',['-r',num2str(rez)],'-opengl') %save file 
  
close(gcf)

end

%% bed open ONLY
for i=1:1:5
fig1 = figure('Units','centimeters','Position',[0.5 0.5 25.5 22]);%,...
                 % 'Name',[titlestring ' 2'],'NumberTitle','off');
    clf
    axe1 = axes('Position',[0.059 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe2 = axes('Position',[0.059 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe3 = axes('Position',[0.059 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe4 = axes('Position',[0.379 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe5 = axes('Position',[0.379 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe6 = axes('Position',[0.379 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
    
    axe7 = axes('Position',[0.699 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe8 = axes('Position',[0.699 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe9 = axes('Position',[0.699 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
 %   
% bed open --> surface displacement
axes(axe1)
vv=[-1.5025:0.005:1.5025]; 
[C1,h1]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_E_open(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-0.5 0.5]);
colormap(axe1,BWR); 
t1=colorbar('EastOutside'); set(t1,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t1,'ylabel'),'String','x-displacement  [ m ]','FontSize',10);
set(t1, 'Position', [.295 0.70 .005 (0.27)-0.02]);

text(-29, 20.2, 'a.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Bed Open: Displacement, U','FontSize',10)
text(-21,19,sprintf('H = 2500 m   V = %4.3f km^{3}',Vtot(i)),'FontWeight','bold','FontSize',10);
% axis equal

axes(axe2)
[C2,h2]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_N_open(i,:),121,121),vv);
set(h2,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar;  
caxis([-0.5 0.5]);
colormap(axe2,BWR); 
t2=colorbar('EastOutside'); set(t2,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t2,'ylabel'),'String','y-displacement  [ m ]','FontSize',10);
set(t2, 'Position', [.295 0.3825 .005 (0.27)-0.02]);

text(-29, 20.2, 'b.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% axis equal

axes(axe3)
vvv=[-1.0025:0.005:1.0025];
[C3,h3]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_U_open(i,:),121,121),vvv);
set(h3,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-1 1]);
colormap(axe3,BWR); 
t3=colorbar('EastOutside'); set(t3,'YTick',[-1,-0.5,0,0.5,1]); hold all
set(get(t3,'ylabel'),'String','z-displacement  [ m ]','FontSize',10);
set(t3, 'Position', [.295 0.065 .005 (0.27)-0.02]);
% axis equal  

text(-29, 20.2, 'c.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10); 
xlabel(' x [ km ]','FontSize',10);
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
  
% bed open --> surface strain
axes(axe4)
vv=(10^(-4)).*[-1.5025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_xx_open(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe4,BWR); 
t4=colorbar('EastOutside'); set(t4,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t4,'ylabel'),'String','\epsilon_{xx}  [ m/m ]','FontSize',10);
set(t4, 'Position', [.615 0.70 .005 (0.27)-0.02]);
text(-24, 20.2, 'd.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Bed Open: Strain, \epsilon','FontSize',10)
% axis equal

axes(axe5)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_yy_open(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe5,BWR); 
t5=colorbar('EastOutside'); set(t5,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t5,'ylabel'),'String','\epsilon_{yy}  [ m/m ]','FontSize',10);
set(t5, 'Position', [.615 0.3825 .005 (0.27)-0.02]);
text(-24, 20.2, 'e.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
% axis equal

axes(axe6)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_zz_open(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe6,BWR); 
t6=colorbar('EastOutside'); set(t6,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t6,'ylabel'),'String','\epsilon_{zz}  [ m/m ]','FontSize',10);
set(t6, 'Position', [.615 0.065 .005 (0.27)-0.02]);
text(-24, 20.2, 'f.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlabel(' x [ km ]','FontSize',10);
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
% axis equal

% STRESSES sigma_xx
axes(axe7)
vvKPA=[-8005:10:8005]; v150=150;
[C4,h4]=contourf(xx, yy, reshape(sigma_xx_open(i,:),121,121)./1e3,vvKPA);
set(h4,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_xx_open(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-600 600]);
colormap(axe7,BWR); 
t7=colorbar('EastOutside'); set(t7,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t7,'ylabel'),'String','\sigma_{xx}  [ kPa ]','FontSize',10);
set(t7,'Position', [.935 0.70 .005 (0.27)-0.02]);

text(-24, 20.2, 'g.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Bed Open: Stress, \sigma','FontSize',10)

% sigma_yy
axes(axe8)
[C5,h5]=contourf(xx, yy, reshape(sigma_yy_open(i,:),121,121)./1e3,vvKPA);
set(h5,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_yy_open(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-600 600]);
colormap(axe8,BWR); 
t8=colorbar('EastOutside'); set(t8,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t8,'ylabel'),'String','\sigma_{yy}  [ kPa ]','FontSize',10);
set(t8,'Position', [.935 0.3825 .005 (0.27)-0.02]);

text(-24, 20.2, 'h.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% sigma_zz
axes(axe9)
[C5,h6]=contourf(xx, yy, reshape(sigma_zz_open(i,:),121,121)./1e3,vvKPA);
set(h6,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_zz_open(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-600 600]);
colormap(axe9,BWR); 
t9=colorbar('EastOutside'); set(t9,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t9,'ylabel'),'String','\sigma_{zz}  [ kPa ]','FontSize',10);
set(t9, 'Position', [.935 0.065 .005 (0.27)-0.02]);

text(-24, 20.2, 'i.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
xlabel(' x [ km ]','FontSize',10);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% print movie frame
rez=300; %resolution (dpi) of final graphic
figpos=getpixelposition(fig1); %dont need to change anything here
resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
set(fig1,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
path='Forward_2500m_check_plots'; %the folder where you want to put the file
name=sprintf('%04d_2500m_bedopen',i); %what you want the file to be called
print(fig1,fullfile(path,name),'-dpng',['-r',num2str(rez)],'-opengl') %save file 
  
close(gcf)

end


%% bed slip ONLY
for i=1:1:5
fig1 = figure('Units','centimeters','Position',[0.5 0.5 25.5 22]);%,...
                 % 'Name',[titlestring ' 2'],'NumberTitle','off');
    clf
    axe1 = axes('Position',[0.059 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe2 = axes('Position',[0.059 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe3 = axes('Position',[0.059 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    %set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe4 = axes('Position',[0.379 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe5 = axes('Position',[0.379 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe6 = axes('Position',[0.379 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
    
    axe7 = axes('Position',[0.699 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe8 = axes('Position',[0.699 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe9 = axes('Position',[0.699 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    %set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [],'FontSize',12);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
 %   
% bed slip --> surface displacement
axes(axe1)
vv=[-1.5025:0.005:1.5025]; 
[C1,h1]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_E_slip(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-0.5 0.5]);
colormap(axe1,BWR); 
t1=colorbar('EastOutside'); set(t1,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t1,'ylabel'),'String','x-displacement  [ m ]','FontSize',10);
set(t1, 'Position', [.295 0.70 .005 (0.27)-0.02]);

text(-29, 20.2, 'a.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Bed Slip: Displacement, U','FontSize',10)
text(-21,19,sprintf('H = 2500 m   V = %4.3f km^{3}',Vtot(i)),'FontWeight','bold','FontSize',10);
% axis equal

axes(axe2)
[C2,h2]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_N_slip(i,:),121,121),vv);
set(h2,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar;  
caxis([-0.5 0.5]);
colormap(axe2,BWR); 
t2=colorbar('EastOutside'); set(t2,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t2,'ylabel'),'String','y-displacement  [ m ]','FontSize',10);
set(t2, 'Position', [.295 0.3825 .005 (0.27)-0.02]);

text(-29, 20.2, 'b.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% axis equal

axes(axe3)
vvv=[-1.0025:0.005:1.0025];
[C3,h3]=contourf(xx, yy, reshape(surface_disp_2500m.Disp_U_slip(i,:),121,121),vvv);
set(h3,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-1 1]);
colormap(axe3,BWR); 
t3=colorbar('EastOutside'); set(t3,'YTick',[-1,-0.5,0,0.5,1]); hold all
set(get(t3,'ylabel'),'String','z-displacement  [ m ]','FontSize',10);
set(t3, 'Position', [.295 0.065 .005 (0.27)-0.02]);
% axis equal  

text(-29, 20.2, 'c.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10); 
xlabel(' x [ km ]','FontSize',10);
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
  
% bed slip --> surface strain
axes(axe4)
vv=(10^(-4)).*[-1.5025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_xx_slip(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe4,BWR); 
t4=colorbar('EastOutside'); set(t4,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t4,'ylabel'),'String','\epsilon_{xx}  [ m/m ]','FontSize',10);
set(t4, 'Position', [.615 0.70 .005 (0.27)-0.02]);
text(-24, 20.2, 'd.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

title('Bed Slip: Strain, \epsilon','FontSize',10)
% axis equal

axes(axe5)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_yy_slip(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe5,BWR); 
t5=colorbar('EastOutside'); set(t5,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t5,'ylabel'),'String','\epsilon_{yy}  [ m/m ]','FontSize',10);
set(t5, 'Position', [.615 0.3825 .005 (0.27)-0.02]);
text(-24, 20.2, 'e.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
% axis equal

axes(axe6)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_zz_slip(i,:),121,121),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-0.0001 0.0001]);
colormap(axe6,BWR); 
t6=colorbar('EastOutside'); set(t6,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t6,'ylabel'),'String','\epsilon_{zz}  [ m/m ]','FontSize',10);
set(t6, 'Position', [.615 0.065 .005 (0.27)-0.02]);
text(-24, 20.2, 'f.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlabel(' x [ km ]','FontSize',10);
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
% axis equal

% STRESSES sigma_xx
axes(axe7)
vvKPA=[-8005:10:8005]; v150=150;
[C4,h4]=contourf(xx, yy, reshape(sigma_xx_slip(i,:),121,121)./1e3,vvKPA);
set(h4,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_xx_slip(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-600 600]);
colormap(axe7,BWR); 
t7=colorbar('EastOutside'); set(t7,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t7,'ylabel'),'String','\sigma_{xx}  [ kPa ]','FontSize',10);
set(t7,'Position', [.935 0.70 .005 (0.27)-0.02]);

text(-24, 20.2, 'g.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('Bed Slip: Stress, \sigma','FontSize',10)


% sigma_yy
axes(axe8)
[C5,h5]=contourf(xx, yy, reshape(sigma_yy_slip(i,:),121,121)./1e3,vvKPA);
set(h5,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_yy_slip(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%scatter(xy_moulin(1,1),xy_moulin(1,2),10,'yo','filled','MarkerEdgeColor','k')

%colorbar; 
caxis([-600 600]);
colormap(axe8,BWR); 
t8=colorbar('EastOutside'); set(t8,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t8,'ylabel'),'String','\sigma_{yy}  [ kPa ]','FontSize',10);
set(t8,'Position', [.935 0.3825 .005 (0.27)-0.02]);

text(-24, 20.2, 'h.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% sigma_zz
axes(axe9)
[C5,h6]=contourf(xx, yy, reshape(sigma_zz_slip(i,:),121,121)./1e3,vvKPA);
set(h6,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_zz_slip(i,:),121,121)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
%scatter(xy_lake(:,1),xy_lake(:,2),2,'bo','filled','MarkerEdgeColor','none')
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
%colorbar; 
caxis([-600 600]);
colormap(axe9,BWR); 
t9=colorbar('EastOutside'); set(t9,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t9,'ylabel'),'String','\sigma_{zz}  [ kPa ]','FontSize',10);
set(t9, 'Position', [.935 0.065 .005 (0.27)-0.02]);

text(-24, 20.2, 'i.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-15 15]); ylim([-15 15]);
xlabel(' x [ km ]','FontSize',10);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% print movie frame
rez=300; %resolution (dpi) of final graphic
figpos=getpixelposition(fig1); %dont need to change anything here
resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
set(fig1,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
path='Forward_2500m_check_plots'; %the folder where you want to put the file
name=sprintf('%04d_2500m_bedslip',i); %what you want the file to be called
print(fig1,fullfile(path,name),'-dpng',['-r',num2str(rez)],'-opengl') %save file 
  
close(gcf)

end
