%% Forwards model surface displacements and stress over 20 km x 20 km region 
% at the half-space surface (used saved NIF L1A outputs for 20 km x 20 km 
% bed plane and vertical crack along L1A hydro-fracture scarp).
% Laura A. Stevens Oct 26, 2015 (Oxford visit)
% 2021 Mar 04 LAS -- now with strains!
% 2023 May 31 LAS -- updated with Stacy Larochelle infinitesimal strain 
% calcs; finalized for revisions; vary \lambda and \mu to range over 
% material properties. 

clear all; close all;
path(path,'software/general')
path(path,'software/functions/dataIO')
path(path,'software/objects')
path(path,'software/functions/elmat')
path(path,'software/functions/time')
path(path,'software/filter_tools')
path(path,'software/matools')
path(path,'software/okada85')

%% 2011 NIF L1A outputs
load MLE_47_all_6z_2_20kmB.mat % NIF outputs for 20 km bed; L1A set-up.
load Gsurface500_strain.mat % Green functions for 20 km X 20 km surface using
% 0.5-km spacing across surface and the NL set-up

%% extract opening and slip along planes
time_2011 = ts_enu.epochs; % time vector
vert_crack_open_m = slip_b(1:Nsubfa,:)'; % m
bed_crack_open_m = slip_b(Nsubfa+1:Nsubfa+Nsubfb,:)'; % m
bed_crack_slip_m = slip_b(Nsubfa+Nsubfb+1:Nsubfa+Nsubfb+Nsubfc,:)'; % m

%% extract deformation at GPS station locations (not needed here)
% for i=1:4
% II=unique(ts_enu(i,:).epochindex);
% 
% inde=(i-1)*3+1;
% east_2011(i,:).site = ts_enu(i,:).d(1:3:end)-Vhat(inde,II)-F(1,II);
% east_2011(i,:).time = ts2(i,:).epochs;
% indn=(i-1)*3+2;
% north_2011(i,:).site = ts_enu(i,:).d(2:3:end)-Vhat(indn,II)-F(2,II);
% north_2011(i,:).time = ts2(i,:).epochs;
% indu=(i-1)*3+3;
% up_2011(i,:).site = ts_enu(i,:).d(3:3:end)-Vhat(indu,II)-F(3,II);
% up_2011(i,:).time = ts2(i,:).epochs;
% end

%% calculate surface displacements (forward problem)
% Gsurface500 (green's function matrix) --> displacement
% xy for the 20km x 20 km surface using meshgrid(-20:0.5:20, -20:0.5:20); [ km ]
xy_surf = Gsurface500.xy_surf; %horzcat(reshape(xx,41*41,1), reshape(yy,41*41,1));
Nsurface = Gsurface500.Nsurface; %length(xy_surf);

for j=1:length(vert_crack_open_m(:,1))
    % crack opening
    B_G3v = Gsurface500.G3v*vert_crack_open_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;

        B_G3v_E(j,i) = B_G3v(ind1); %E
        B_G3v_N(j,i) = B_G3v(ind2); %N
        B_G3v_U(j,i) = B_G3v(ind3); %U
    end
end

for j=1:length(bed_crack_open_m(:,1))
    % bed opening
    B_G3bv = Gsurface500.G3bv*bed_crack_open_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;

        B_G3bv_E(j,i) = B_G3bv(ind1); %E
        B_G3bv_N(j,i) = B_G3bv(ind2); %N
        B_G3bv_U(j,i) = B_G3bv(ind3); %U
    end
end

for j=1:length(bed_crack_slip_m(:,1))
    % bed slip
    B_G2bv = Gsurface500.G2bv*bed_crack_slip_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;

        B_G2bv_E(j,i) = B_G2bv(ind1); %E
        B_G2bv_N(j,i) = B_G2bv(ind2); %N
        B_G2bv_U(j,i) = B_G2bv(ind3); %U
    end
end
 
%% calculate surface strain components (forward problem)
% Gsurface500 (green's function matrix) --> strains and "tilts" (SL)
for j=1:length(vert_crack_open_m(:,1))
    % crack opening (multiply by -1 to get sign convention to positive strain = TENSION)
    B_G3v = -1.*Gsurface500.G3v_strain*vert_crack_open_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;
        ind4=Nsurface*3+i; ind5=Nsurface*4+i; ind6=Nsurface*5+i;
        
        B_G3v_uZE(j,i) = B_G3v(ind1); %uZE = duz/dx
        B_G3v_uZN(j,i) = B_G3v(ind2); %uZN = duz/dy
        B_G3v_uNN(j,i) = B_G3v(ind3); %uNN = duy/dy
        B_G3v_uNE(j,i) = B_G3v(ind4); %uNE = duy/dx
        B_G3v_uEN(j,i) = B_G3v(ind5); %uEN = dux/dy
        B_G3v_uEE(j,i) = B_G3v(ind6); %uEE = dux/dx
    end
end

for j=1:length(bed_crack_open_m(:,1))
    % bed opening (multiply by -1 to get sign convention to positive strain = TENSION)
    B_G3bv = -1.*Gsurface500.G3bv_strain*bed_crack_open_m(j,:)';
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
    B_G2bv = -1.*Gsurface500.G2bv_strain*bed_crack_slip_m(j,:)';
    for i=1:Nsurface
        ind1=i; ind2=Nsurface+i; ind3=Nsurface*2+i;
        ind4=Nsurface*3+i; ind5=Nsurface*4+i; ind6=Nsurface*5+i;
        
        B_G2bv_uZE(j,i) = B_G2bv(ind1); %uZE = duz/dx
        B_G2bv_uZN(j,i) = B_G2bv(ind2); %uZN = duz/dy
        B_G2bv_uNN(j,i) = B_G2bv(ind3); %uNN = duy/dy
        B_G2bv_uNE(j,i) = B_G2bv(ind4); %uNE = duy/dx
        B_G2bv_uEN(j,i) = B_G2bv(ind5); %uEN = dux/dy
        B_G2bv_uEE(j,i) = B_G2bv(ind6); %uEE = dux/dx
    end
end

% linear elastic: add all sources of strain 
NU = 0.3; % Poisson ratio

% vert open (SL)
duz_dx_vert =  (B_G3v_uZE)./1e3; % u in m but dx is in km
duz_dy_vert =  (B_G3v_uZN)./1e3;
dux_dx_vert =  (B_G3v_uEE)./1e3;
duy_dx_vert =  (B_G3v_uNE)./1e3;
duy_dy_vert =  (B_G3v_uNN)./1e3;
dux_dy_vert =  (B_G3v_uEN)./1e3;
dux_dz_vert =  -1.*duz_dx_vert; % From free surface boundary condition (strain_xz = 0 at the surface)
duy_dz_vert =  -1.*duz_dy_vert; % From free surface boundary condition (strain_yz = 0 at the surface)
duz_dz_vert =  -1.*(dux_dx_vert + duy_dy_vert).*(NU/(1-NU)); % From free surface boundary condition (stress_zz = 0 at the surface)

% bed open (SL)
duz_dx_open =  (B_G3bv_uZE)./1e3; % u in m but dx is in km
duz_dy_open =  (B_G3bv_uZN)./1e3;
dux_dx_open =  (B_G3bv_uEE)./1e3;
duy_dx_open =  (B_G3bv_uNE)./1e3;
duy_dy_open =  (B_G3bv_uNN)./1e3;
dux_dy_open =  (B_G3bv_uEN)./1e3;
dux_dz_open =  -1.*duz_dx_open; % From free surface boundary condition (strain_xz = 0 at the surface)
duy_dz_open =  -1.*duz_dy_open; % From free surface boundary condition (strain_yz = 0 at the surface)
duz_dz_open =  -1.*(dux_dx_open + duy_dy_open).*(NU/(1-NU)); % From free surface boundary condition (stress_zz = 0 at the surface)

% bed slip (SL)
duz_dx_slip =   (B_G2bv_uZE)./1e3; % u in m but dx is in km
duz_dy_slip =   (B_G2bv_uZN)./1e3;
dux_dx_slip =   (B_G2bv_uEE)./1e3;
duy_dx_slip =   (B_G2bv_uNE)./1e3;
duy_dy_slip =   (B_G2bv_uNN)./1e3;
dux_dy_slip =   (B_G2bv_uEN)./1e3;
dux_dz_slip =  -1.*duz_dx_slip; % From free surface boundary condition (strain_xz = 0 at the surface)
duy_dz_slip =  -1.*duz_dy_slip; % From free surface boundary condition (strain_yz = 0 at the surface)
duz_dz_slip =  -1.*(dux_dx_slip + duy_dy_slip).*(NU/(1-NU)); % From free surface boundary condition (stress_zz = 0 at the surface)

%% Vertical opening
% (Infinitesimal) strain tensor at the surface (SL)
% eps_ij = 0.5*(dui_dj + duj_di)
% eps is short for epsilon

eps_xx_vert = dux_dx_vert; eps_xy_vert = 0.5*(dux_dy_vert+duy_dx_vert); eps_xz_vert = 0.5*(dux_dz_vert+duz_dx_vert);
eps_yx_vert = 0.5*(duy_dx_vert+dux_dy_vert); eps_yy_vert = duy_dy_vert; eps_yz_vert = 0.5*(duy_dz_vert+duz_dy_vert);
eps_zx_vert = 0.5*(duz_dx_vert+dux_dz_vert); eps_zy_vert = 0.5*(duz_dy_vert+duy_dz_vert); eps_zz_vert = duz_dz_vert;

% Stress tensor at the surface (SL)
% sigma_ij = lambda * delta_ij * eps_kk + 2*mu*eps_ij
% where eps_kk = eps_xx + eps_yy + eps_zz is the volumetric strain
eps_kk_vert = eps_xx_vert + eps_yy_vert + eps_zz_vert; 

% shear modulus (mu) middle estimate
lambda = 2.25e9; 
mu = 1.5e9;

% shear modulus (mu) lower estimate
% lambda = 0.48e9; 
% mu = 0.32e9;

% shear modulus (mu) upper estimate
% lambda = 5.85e9; 
% mu = 3.9e9;

sigma_xx_vert = lambda*eps_kk_vert+2*mu*eps_xx_vert; sigma_xy_vert = 2*mu*eps_xy_vert; sigma_xz_vert = 2*mu*eps_xz_vert;
sigma_yx_vert = 2*mu*eps_yx_vert; sigma_yy_vert = lambda*eps_kk_vert+2*mu*eps_yy_vert; sigma_yz_vert = 2*mu*eps_yz_vert;
sigma_zx_vert = 2*mu*eps_zx_vert; sigma_zy_vert = 2*mu*eps_zy_vert; sigma_zz_vert = lambda*eps_kk_vert+2*mu*eps_zz_vert;

%% Slip bed
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

%% Open bed
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

%% save sigma_ij for plotting 

% HF (just the vertical hydro-fracture)
sigma_ij_2011_SLder_HF.sigma_xx_vert = sigma_xx_vert;
sigma_ij_2011_SLder_HF.sigma_xy_vert = sigma_xy_vert;
sigma_ij_2011_SLder_HF.sigma_xz_vert = sigma_xz_vert;
sigma_ij_2011_SLder_HF.sigma_yx_vert = sigma_yx_vert;
sigma_ij_2011_SLder_HF.sigma_yy_vert = sigma_yy_vert;
sigma_ij_2011_SLder_HF.sigma_yz_vert = sigma_yz_vert;
sigma_ij_2011_SLder_HF.sigma_zx_vert = sigma_zx_vert;
sigma_ij_2011_SLder_HF.sigma_zy_vert = sigma_zy_vert;
sigma_ij_2011_SLder_HF.sigma_zz_vert = sigma_zz_vert;
sigma_ij_2011_SLder_HF.sigma1 = (0.5.*(sigma_xx_vert+sigma_yy_vert)) + ...
    (sqrt(((0.5.*(sigma_xx_vert-sigma_yy_vert)).^2) + (sigma_xy_vert.^2)));
sigma_ij_2011_SLder_HF.sigma2 = (0.5.*(sigma_xx_vert+sigma_yy_vert)) - ...
    (sqrt(((0.5.*(sigma_xx_vert-sigma_yy_vert)).^2) + (sigma_xy_vert.^2)));
sigma_ij_2011_SLder_HF.theta = 0.5.*(atan2((2.*sigma_xy_vert),(sigma_xx_vert-sigma_yy_vert)));
%save sigma_ij_2011_SLder_HF.mat sigma_ij_2011_SLder_HF

% bedslip (just bedslip)
sigma_ij_2011_SLder_bedslip.sigma_xx_slip = sigma_xx_slip;
sigma_ij_2011_SLder_bedslip.sigma_xy_slip = sigma_xy_slip;
sigma_ij_2011_SLder_bedslip.sigma_xz_slip = sigma_xz_slip;
sigma_ij_2011_SLder_bedslip.sigma_yx_slip = sigma_yx_slip;
sigma_ij_2011_SLder_bedslip.sigma_yy_slip = sigma_yy_slip;
sigma_ij_2011_SLder_bedslip.sigma_yz_slip = sigma_yz_slip;
sigma_ij_2011_SLder_bedslip.sigma_zx_slip = sigma_zx_slip;
sigma_ij_2011_SLder_bedslip.sigma_zy_slip = sigma_zy_slip;
sigma_ij_2011_SLder_bedslip.sigma_zz_slip = sigma_zz_slip;
sigma_ij_2011_SLder_bedslip.sigma1 = (0.5.*(sigma_xx_slip+sigma_yy_slip)) + ...
    (sqrt(((0.5.*(sigma_xx_slip-sigma_yy_slip)).^2) + (sigma_xy_slip.^2)));
sigma_ij_2011_SLder_bedslip.sigma2 = (0.5.*(sigma_xx_slip+sigma_yy_slip)) - ...
    (sqrt(((0.5.*(sigma_xx_slip-sigma_yy_slip)).^2) + (sigma_xy_slip.^2)));
sigma_ij_2011_SLder_bedslip.theta = 0.5.*(atan2((2.*sigma_xy_slip),(sigma_xx_slip-sigma_yy_slip)));
%save sigma_ij_2011_SLder_bedslip.mat sigma_ij_2011_SLder_bedslip

% bedopen (just bedopen)
sigma_ij_2011_SLder_bedopen.sigma_xx_open = sigma_xx_open;
sigma_ij_2011_SLder_bedopen.sigma_xy_open = sigma_xy_open;
sigma_ij_2011_SLder_bedopen.sigma_xz_open = sigma_xz_open;
sigma_ij_2011_SLder_bedopen.sigma_yx_open = sigma_yx_open;
sigma_ij_2011_SLder_bedopen.sigma_yy_open = sigma_yy_open;
sigma_ij_2011_SLder_bedopen.sigma_yz_open = sigma_yz_open;
sigma_ij_2011_SLder_bedopen.sigma_zx_open = sigma_zx_open;
sigma_ij_2011_SLder_bedopen.sigma_zy_open = sigma_zy_open;
sigma_ij_2011_SLder_bedopen.sigma_zz_open = sigma_zz_open;
sigma_ij_2011_SLder_bedopen.sigma1 = (0.5.*(sigma_xx_open+sigma_yy_open)) + ...
    (sqrt(((0.5.*(sigma_xx_open-sigma_yy_open)).^2) + (sigma_xy_open.^2)));
sigma_ij_2011_SLder_bedopen.sigma2 = (0.5.*(sigma_xx_open+sigma_yy_open)) - ...
    (sqrt(((0.5.*(sigma_xx_open-sigma_yy_open)).^2) + (sigma_xy_open.^2)));
sigma_ij_2011_SLder_bedopen.theta = 0.5.*(atan2((2.*sigma_xy_open),(sigma_xx_open-sigma_yy_open)));
%save sigma_ij_2011_SLder_bedopen.mat sigma_ij_2011_SLder_bedopen

% %% load pre-saved stress calculations
% load sigma_ij_2011_SLder_HF.mat
% load sigma_ij_2011_SLder_bedopen.mat
% load sigma_ij_2011_SLder_bedslip.mat

%% displacement, strains, stresses check figure
load apcoords_lle
    lats=apcoords_lle(1,:); lons=apcoords_lle(2,:); hs=apcoords_lle(3,:); 
    llh=[lats; lons; hs];
    xy_sta_11=llh2localxy(llh,origin);
    load moulin_2011.mat % 2011 L1A moulin location
    llh_moulin = [moulin_2011(:,4)'; moulin_2011(:,3)'; zeros(37,1)'];
    xy_moulin = llh2localxy(llh_moulin,origin);
    load lake.mat; llh_lake = [lake(:,5)'; lake(:,4)'; zeros(337,1)'];
    xy_lake = llh2localxy(llh_lake,origin);
    [xx,yy] = meshgrid(-20:0.5:20, -20:0.5:20); % [ km ]
    load BWR.mat
    TriangleSize = 7;

%%%%%%%% plot forward problem surface displacements as a check
for i=1280 % 2011 t_{3} index
fig1 = figure('Units','centimeters','Position',[0.5 0.5 25.5 22]);
    clf
    axe1 = axes('Position',[0.059 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    
    axe2 = axes('Position',[0.059 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    
    axe3 = axes('Position',[0.059 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    
    axe4 = axes('Position',[0.379 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe5 = axes('Position',[0.379 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe6 = axes('Position',[0.379 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
    
    axe7 = axes('Position',[0.699 0.69 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);  
    
    axe8 = axes('Position',[0.699 0.3775 0.23 0.27],'Box','on','NextPlot','add');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', [], 'XTick', []);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []);
    
    axe9 = axes('Position',[0.699 0.06 0.23 0.27],'Box','on','NextPlot','add'); 
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', [], 'YTick', []); 
   
% crack open, bed open, bed slip --> surface displacement
axes(axe1)
vv=[-0.5025:0.005:0.5025]; 
[C1,h1]=contourf(xx, yy, reshape(B_G3v_E(i,:)+B_G2bv_E(i,:)+B_G3bv_E(i,:),81,81),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-0.5 0.5]);
colormap(axe1,BWR); 
t1=colorbar('EastOutside'); set(t1,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t1,'ylabel'),'String','x-displacement  [ m ]','FontSize',10);
set(t1, 'Position', [.295 0.70 .005 (0.27)-0.02]);
text(-29, 20.2, 'a.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('Displacements','FontSize',10)
text(-29,24.5,sprintf('2011/%6.3f',time_2011(i)),'FontWeight','bold','FontSize',10);

axes(axe2)
[C2,h2]=contourf(xx, yy, reshape(B_G3v_N(i,:)+B_G2bv_N(i,:)+B_G3bv_N(i,:),81,81),vv);
set(h2,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)  
caxis([-0.5 0.5]);
colormap(axe2,BWR); 
t2=colorbar('EastOutside'); set(t2,'YTick',[-0.5,-0.25,0,0.25,0.5]); hold all
set(get(t2,'ylabel'),'String','y-displacement  [ m ]','FontSize',10);
set(t2, 'Position', [.295 0.3825 .005 (0.27)-0.02]);
text(-29, 20.2, 'b.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10);% xlabel(' x [ km ]');
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

axes(axe3)
vvv=[-1.0025:0.005:1.0025];
[C3,h3]=contourf(xx, yy, reshape(B_G3v_U(i,:)+B_G2bv_U(i,:)+B_G3bv_U(i,:),81,81),vvv);
set(h3,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-1 1]);
colormap(axe3,BWR); 
t3=colorbar('EastOutside'); set(t3,'YTick',[-1,-0.5,0,0.5,1]); hold all
set(get(t3,'ylabel'),'String','z-displacement  [ m ]','FontSize',10);
set(t3, 'Position', [.295 0.065 .005 (0.27)-0.02]);
text(-29, 20.2, 'c.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
ylabel(' y [ km ]','FontSize',10); 
xlabel(' x [ km ]','FontSize',10);
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
  
% crack open, bed open, bed slip --> surface strain
axes(axe4)
vv=(10^(-4)).*[-1.5025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_xx_vert(i,:)+eps_xx_open(i,:)+eps_xx_slip(i,:),81,81),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-0.0001 0.0001]);
colormap(axe4,BWR); 
t4=colorbar('EastOutside'); set(t4,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t4,'ylabel'),'String','\epsilon_{xx}  [ m/m ]','FontSize',10);
set(t4, 'Position', [.615 0.70 .005 (0.27)-0.02]);
text(-24, 20.2, 'd.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('Strains, \epsilon','FontSize',10)

axes(axe5)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_yy_vert(i,:)+eps_yy_open(i,:)+eps_yy_slip(i,:),81,81),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-0.0001 0.0001]);
colormap(axe5,BWR); 
t5=colorbar('EastOutside'); set(t5,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t5,'ylabel'),'String','\epsilon_{yy}  [ m/m ]','FontSize',10);
set(t5, 'Position', [.615 0.3825 .005 (0.27)-0.02]);
text(-24, 20.2, 'e.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

axes(axe6)
vv=(10^(-4)).*[-0.8025:0.005:0.8025]; 
[C1,h1]=contourf(xx, yy, reshape(eps_zz_vert(i,:)+eps_zz_open(i,:)+eps_zz_slip(i,:),81,81),vv);
set(h1,'LineColor','none'); 
hold on;
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5) 
caxis([-0.0001 0.0001]);
colormap(axe6,BWR); 
t6=colorbar('EastOutside'); set(t6,'YTick',[-0.0001,-0.00005,0,0.00005,0.0001]); hold all
set(get(t6,'ylabel'),'String','\epsilon_{zz}  [ m/m ]','FontSize',10);
set(t6, 'Position', [.615 0.065 .005 (0.27)-0.02]);
text(-24, 20.2, 'f.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlabel(' x [ km ]','FontSize',10);
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% STRESSES sigma_xx
axes(axe7)
vvKPA=[-8005:10:8005]; v150=150;
[C4,h4]=contourf(xx, yy, reshape(sigma_xx_vert(i,:)+sigma_xx_open(i,:)+sigma_xx_slip(i,:),81,81)./1e3,vvKPA);
set(h4,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_xx_vert(i,:)+sigma_xx_open(i,:)+sigma_xx_slip(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-600 600]);
colormap(axe7,BWR); 
t7=colorbar('EastOutside'); set(t7,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t7,'ylabel'),'String','\sigma_{xx}  [ kPa ]','FontSize',10);
set(t7,'Position', [.935 0.70 .005 (0.27)-0.02]);
text(-24, 20.2, 'g.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;
title('Stresses, \sigma','FontSize',10)

% sigma_yy
axes(axe8)
[C5,h5]=contourf(xx, yy, reshape(sigma_yy_vert(i,:)+sigma_yy_open(i,:)+sigma_yy_slip(i,:),81,81)./1e3,vvKPA);
set(h5,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_yy_vert(i,:)+sigma_yy_open(i,:)+sigma_yy_slip(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-600 600]);
colormap(axe8,BWR); 
t8=colorbar('EastOutside'); set(t8,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t8,'ylabel'),'String','\sigma_{yy}  [ kPa ]','FontSize',10);
set(t8,'Position', [.935 0.3825 .005 (0.27)-0.02]);
text(-24, 20.2, 'h.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-20 20]); ylim([-20 20]);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% sigma_zz
axes(axe9)
[C5,h6]=contourf(xx, yy, reshape(sigma_zz_vert(i,:)+sigma_zz_open(i,:)+sigma_zz_slip(i,:),81,81)./1e3,vvKPA);
set(h6,'LineColor','none'); 
hold on;
contour(xx, yy,reshape(sigma_zz_vert(i,:)+sigma_zz_open(i,:)+sigma_zz_slip(i,:),81,81)./1e3,[v150 v150],'k','LineWidth',1.1);
plot(xy_sta_11(:,1),xy_sta_11(:,2),'k^','MarkerSize',TriangleSize-5,'MarkerFaceColor','none'); hold on;
plot(patchesC(1:6:end,6),patchesC(1:6:end,7),'y','LineWidth',1.5)
caxis([-600 600]);
colormap(axe9,BWR); 
t9=colorbar('EastOutside'); set(t9,'YTick',[-600,-450,-300,-150,0,150,300,450,600]); 
hold all
set(get(t9,'ylabel'),'String','\sigma_{zz}  [ kPa ]','FontSize',10);
set(t9, 'Position', [.935 0.065 .005 (0.27)-0.02]);
text(-24, 20.2, 'i.','FontSize',10,'FontWeight','bold');
set(gca,'FontName','Avenir');
xlim([-20 20]); ylim([-20 20]);
xlabel(' x [ km ]','FontSize',10);
set(gca,'xtick',[-20:10:20],'ytick',[-20:10:20],'tickdir','out','LineWidth',1.1,'FontSize',10); 
grid on;

% print movie frame
rez=300; %resolution (dpi) of final graphic
figpos=getpixelposition(fig1); %dont need to change anything here
resolution=get(0,'ScreenPixelsPerInch'); %dont need to change anything here
set(fig1,'paperunits','inches','papersize',figpos(3:4)/resolution,'paperposition',[0 0 figpos(3:4)/resolution]); %dont need to change anything here
name=sprintf('check_fig_2011_%04d',i); 
print(fig1,fullfile(name),'-dpng',['-r',num2str(rez)],'-opengl') %save file 
print(gcf,'-dpng','-r600',name);  
close(gcf)

end
