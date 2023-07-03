%% Idealized blister simulations with derivation from CYL
% LAS 22 March 2021
% LAS 24 April 2021 -- basal slip of 1 m when r<R+1
% LAS 06 March 2023 -- review for sharing with S. Larochelle
% LAS 30 June 2023 -- edits for review, patch center location
clear all; close all;

%% Lai et al 2021 for L1A 2011 values:
alpha = 0.32; % [ ] from Lai 2021
hmax2011 = 1.03; % [ m ] from Lai 2021
Vi_2011 = 0.008*(1000^3); % [ m^3 ] from Lai 2021

%% Volume range for idealized forward experiments
Vi = [0.001, 0.005, 0.01, 0.05, 0.1].*(1e3^3); % [ m^3 ] 
% hmax, R, h from Lai (2021) equations over range in lake volume
hmax = hmax2011.*((Vi./Vi_2011).^(1/3)); % h_max [ m ]
R = sqrt(Vi./(2*pi*alpha.*hmax)); % Radius [ m ]
r = linspace(-10000,10000,401); % r [ m ] 
for i=1:401; h_r(i,:) = hmax.*(sqrt(1-((r(i)./R).^2))); end % [ m ] 
h_r = real(h_r); % h_r [ m ]

%% visualize cavity opening
figure(1); clf
hold on; plot(r./1e3,h_r(:,1),'.-');
plot(r./1e3,h_r(:,2),'.-');
plot(r./1e3,h_r(:,3),'.-');
plot(r./1e3,h_r(:,4),'.-');
plot(r./1e3,h_r(:,5),'.-');
xlabel('r [ km ]'); ylabel('h [ m ]');
legend('0.001 km3','0.005 km3','0.01 km3','0.05 km3','0.1 km3');
title(' h(r,t=0) ' );
% print figure
figurename=sprintf('blister_along_r.png');
print(gcf,'-dpng','-r600',figurename);

%% convert blister opening to NIF input geometry 
% what we would be able to resolve with an equivalent GPS array
load('../forward_okada85/Gsurface500m_strain.mat') % load to get patch locations
patchesB = Gsurface500.patchesB;
spacing = (patchesB(1,1)).*1e3 % [ +500 m ] 
half_spacing = spacing./2 % [ +250 m ] 
patches_center(:,1) = (patchesB(:,6).*1e3) + half_spacing; % [ m ] patch extends in -x direction
patches_center(:,2) = (patchesB(:,7).*1e3); % [ m ] y-location already at lower midpoint

% X and Y for h_r
[X,Y]=meshgrid(-1e5:50:1e5); % [ m ]
RR = sqrt(X.^2 + Y.^2); % [ m ]

% opening - interpolate to patch centers from h_r
for i=1:5
        vq1 = interp1(r, h_r(:,i), RR); % interp to R in XY
        vq2 = interp2(X, Y, vq1,patches_center(:,1), patches_center(:,2));
        opening(:,i) = vq2;
end
opening(isnan(opening))=0; 

% slip - set as 1 m slip for all uplift patches 
slip_logical = gt(opening,0); 
slip = double(slip_logical);

%% visualize slip and cavity opening
figure(2); clf;
scatter(patches_center(:,1),patches_center(:,2),10,opening(:,5),'.');
hold on;
scatter(patches_center(:,1),patches_center(:,2),50,slip(:,5),'o');
colorbar;
axis equal
% print figure
figurename=sprintf('blister_open_slip_NIF_geometry.png');
print(gcf,'-dpng','-r600',figurename);

 %% save idealized patch opening and slip for NIF forward calculation
patches_idealized.opening = opening;
patches_idealized.slip = slip;
patches_idealized.patchesB = patchesB;
patches_idealized.patchesB_center = patches_center;
%save patches_idealized.mat patches_idealized