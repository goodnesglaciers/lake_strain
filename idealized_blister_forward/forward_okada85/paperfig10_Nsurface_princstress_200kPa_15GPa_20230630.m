%% Nsurface principle stress -- JGR:ES paper figure for principal stresses 
% across idealized tested range in ice thickness and lake volume
% 26 Feb 2021 LAS
% 4 March 2021 LAS -- now with a better understanding of strain
% 22 March 2021 LAS -- principle surface stresses across thickness and volume
% 06 March 2023 LAS -- review for sharing with S. Larochelle
% 07 April 2023 LAS -- check with revised strain sign conventions 
% 30 June 2023 LAS -- finalized for revisions
close all;

%% opening and slip along idealized patches
load('../make_idealized_blisters/patches_idealized.mat') % 40 x 40 patch bed

%% surface positions
% surface positions as defined in NIF runfile + makegeom -- do not change !!!
load Gsurface1500m_strain.mat
xy_surf = Gsurface.xy_surf; xx_vec = xy_surf(:,1)'; yy_vec = xy_surf(:,2)'; % [ km ]
xx = Gsurface.xx; yy = Gsurface.yy; Nsurface = Gsurface.Nsurface; % [ km ] 121x121

%% calculate principle stresses from sigma_ij
v200 = [200 200]; % stress threshold

% load saved stresses
load sigma_ij_500m.mat
load sigma_ij_750m.mat
load sigma_ij_1000m.mat
load sigma_ij_1250m.mat
load sigma_ij_1500m.mat
load sigma_ij_1750m.mat
load sigma_ij_2000m.mat
load sigma_ij_2250m.mat
load sigma_ij_2500m.mat
load sigma_ij_2750m.mat
load sigma_ij_3000m.mat

for i = 1:5
     for j = 1:14641
         
        % open and slip 
%500m 
        princ_sigma1_500m(i,j) = (0.5.*(sigma_ij_500m_SLder.sigma_xx(i,j)+sigma_ij_500m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_500m_SLder.sigma_xx(i,j)-sigma_ij_500m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_500m_SLder.sigma_xy(i,j).^2)));   
%750m 
        princ_sigma1_750m(i,j) = (0.5.*(sigma_ij_750m_SLder.sigma_xx(i,j)+sigma_ij_750m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_750m_SLder.sigma_xx(i,j)-sigma_ij_750m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_750m_SLder.sigma_xy(i,j).^2))); 
%1000m 
        princ_sigma1_1000m(i,j) = (0.5.*(sigma_ij_1000m_SLder.sigma_xx(i,j)+sigma_ij_1000m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1000m_SLder.sigma_xx(i,j)-sigma_ij_1000m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_1000m_SLder.sigma_xy(i,j).^2))); 
%1250m 
        princ_sigma1_1250m(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx(i,j)+sigma_ij_1250m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx(i,j)-sigma_ij_1250m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_1250m_SLder.sigma_xy(i,j).^2))); 
%1500m 
        princ_sigma1_1500m(i,j) = (0.5.*(sigma_ij_1500m_SLder.sigma_xx(i,j)+sigma_ij_1500m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1500m_SLder.sigma_xx(i,j)-sigma_ij_1500m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_1500m_SLder.sigma_xy(i,j).^2)));   
%1750m 
        princ_sigma1_1750m(i,j) = (0.5.*(sigma_ij_1750m_SLder.sigma_xx(i,j)+sigma_ij_1750m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1750m_SLder.sigma_xx(i,j)-sigma_ij_1750m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_1750m_SLder.sigma_xy(i,j).^2))); 
%2000m 
        princ_sigma1_2000m(i,j) = (0.5.*(sigma_ij_2000m_SLder.sigma_xx(i,j)+sigma_ij_2000m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2000m_SLder.sigma_xx(i,j)-sigma_ij_2000m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_2000m_SLder.sigma_xy(i,j).^2))); 
%2250m 
        princ_sigma1_2250m(i,j) = (0.5.*(sigma_ij_2250m_SLder.sigma_xx(i,j)+sigma_ij_2250m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2250m_SLder.sigma_xx(i,j)-sigma_ij_2250m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_2250m_SLder.sigma_xy(i,j).^2)));        
%2500m 
        princ_sigma1_2500m(i,j) = (0.5.*(sigma_ij_2500m_SLder.sigma_xx(i,j)+sigma_ij_2500m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2500m_SLder.sigma_xx(i,j)-sigma_ij_2500m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_2500m_SLder.sigma_xy(i,j).^2)));   
%2750m 
        princ_sigma1_2750m(i,j) = (0.5.*(sigma_ij_2750m_SLder.sigma_xx(i,j)+sigma_ij_2750m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2750m_SLder.sigma_xx(i,j)-sigma_ij_2750m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_2750m_SLder.sigma_xy(i,j).^2)));             
%3000m 
        princ_sigma1_3000m(i,j) = (0.5.*(sigma_ij_3000m_SLder.sigma_xx(i,j)+sigma_ij_3000m_SLder.sigma_yy(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_3000m_SLder.sigma_xx(i,j)-sigma_ij_3000m_SLder.sigma_yy(i,j))).^2) + (sigma_ij_3000m_SLder.sigma_xy(i,j).^2)));        
        
   
        % open 
%500m 
        princ_sigma1_500m_open(i,j) = (0.5.*(sigma_ij_500m_SLder.sigma_xx_open(i,j)+sigma_ij_500m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_500m_SLder.sigma_xx_open(i,j)-sigma_ij_500m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_500m_SLder.sigma_xy_open(i,j).^2)));   
%750m 
        princ_sigma1_750m_open(i,j) = (0.5.*(sigma_ij_750m_SLder.sigma_xx_open(i,j)+sigma_ij_750m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_750m_SLder.sigma_xx_open(i,j)-sigma_ij_750m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_750m_SLder.sigma_xy_open(i,j).^2))); 
%1000m 
        princ_sigma1_1000m_open(i,j) = (0.5.*(sigma_ij_1000m_SLder.sigma_xx_open(i,j)+sigma_ij_1000m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1000m_SLder.sigma_xx_open(i,j)-sigma_ij_1000m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_1000m_SLder.sigma_xy_open(i,j).^2))); 
%1250m 
        princ_sigma1_1250m_open(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_open(i,j)+sigma_ij_1250m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_open(i,j)-sigma_ij_1250m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_1250m_SLder.sigma_xy_open(i,j).^2))); 
%1500m 
        princ_sigma1_1500m_open(i,j) = (0.5.*(sigma_ij_1500m_SLder.sigma_xx_open(i,j)+sigma_ij_1500m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1500m_SLder.sigma_xx_open(i,j)-sigma_ij_1500m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_1500m_SLder.sigma_xy_open(i,j).^2)));   
%1750m 
        princ_sigma1_1750m_open(i,j) = (0.5.*(sigma_ij_1750m_SLder.sigma_xx_open(i,j)+sigma_ij_1750m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1750m_SLder.sigma_xx_open(i,j)-sigma_ij_1750m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_1750m_SLder.sigma_xy_open(i,j).^2))); 
%2000m 
        princ_sigma1_2000m_open(i,j) = (0.5.*(sigma_ij_2000m_SLder.sigma_xx_open(i,j)+sigma_ij_2000m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2000m_SLder.sigma_xx_open(i,j)-sigma_ij_2000m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_2000m_SLder.sigma_xy_open(i,j).^2))); 
%2250m 
        princ_sigma1_2250m_open(i,j) = (0.5.*(sigma_ij_2250m_SLder.sigma_xx_open(i,j)+sigma_ij_2250m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2250m_SLder.sigma_xx_open(i,j)-sigma_ij_2250m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_2250m_SLder.sigma_xy_open(i,j).^2)));        
%2500m 
        princ_sigma1_2500m_open(i,j) = (0.5.*(sigma_ij_2500m_SLder.sigma_xx_open(i,j)+sigma_ij_2500m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2500m_SLder.sigma_xx_open(i,j)-sigma_ij_2500m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_2500m_SLder.sigma_xy_open(i,j).^2)));   
%2750m 
        princ_sigma1_2750m_open(i,j) = (0.5.*(sigma_ij_2750m_SLder.sigma_xx_open(i,j)+sigma_ij_2750m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2750m_SLder.sigma_xx_open(i,j)-sigma_ij_2750m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_2750m_SLder.sigma_xy_open(i,j).^2)));             
%3000m 
        princ_sigma1_3000m_open(i,j) = (0.5.*(sigma_ij_3000m_SLder.sigma_xx_open(i,j)+sigma_ij_3000m_SLder.sigma_yy_open(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_3000m_SLder.sigma_xx_open(i,j)-sigma_ij_3000m_SLder.sigma_yy_open(i,j))).^2) + (sigma_ij_3000m_SLder.sigma_xy_open(i,j).^2)));    
        
        % slip 
%500m 
        princ_sigma1_500m_slip(i,j) = (0.5.*(sigma_ij_500m_SLder.sigma_xx_slip(i,j)+sigma_ij_500m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_500m_SLder.sigma_xx_slip(i,j)-sigma_ij_500m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_500m_SLder.sigma_xy_slip(i,j).^2)));   
%750m 
        princ_sigma1_750m_slip(i,j) = (0.5.*(sigma_ij_750m_SLder.sigma_xx_slip(i,j)+sigma_ij_750m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_750m_SLder.sigma_xx_slip(i,j)-sigma_ij_750m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_750m_SLder.sigma_xy_slip(i,j).^2))); 
%1000m 
        princ_sigma1_1000m_slip(i,j) = (0.5.*(sigma_ij_1000m_SLder.sigma_xx_slip(i,j)+sigma_ij_1000m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1000m_SLder.sigma_xx_slip(i,j)-sigma_ij_1000m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_1000m_SLder.sigma_xy_slip(i,j).^2))); 
%1250m 
        princ_sigma1_1250m_slip(i,j) = (0.5.*(sigma_ij_1250m_SLder.sigma_xx_slip(i,j)+sigma_ij_1250m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1250m_SLder.sigma_xx_slip(i,j)-sigma_ij_1250m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_1250m_SLder.sigma_xy_slip(i,j).^2))); 
%1500m 
        princ_sigma1_1500m_slip(i,j) = (0.5.*(sigma_ij_1500m_SLder.sigma_xx_slip(i,j)+sigma_ij_1500m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1500m_SLder.sigma_xx_slip(i,j)-sigma_ij_1500m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_1500m_SLder.sigma_xy_slip(i,j).^2)));   
%1750m 
        princ_sigma1_1750m_slip(i,j) = (0.5.*(sigma_ij_1750m_SLder.sigma_xx_slip(i,j)+sigma_ij_1750m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_1750m_SLder.sigma_xx_slip(i,j)-sigma_ij_1750m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_1750m_SLder.sigma_xy_slip(i,j).^2))); 
%2000m 
        princ_sigma1_2000m_slip(i,j) = (0.5.*(sigma_ij_2000m_SLder.sigma_xx_slip(i,j)+sigma_ij_2000m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2000m_SLder.sigma_xx_slip(i,j)-sigma_ij_2000m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_2000m_SLder.sigma_xy_slip(i,j).^2))); 
%2250m 
        princ_sigma1_2250m_slip(i,j) = (0.5.*(sigma_ij_2250m_SLder.sigma_xx_slip(i,j)+sigma_ij_2250m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2250m_SLder.sigma_xx_slip(i,j)-sigma_ij_2250m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_2250m_SLder.sigma_xy_slip(i,j).^2)));        
%2500m 
        princ_sigma1_2500m_slip(i,j) = (0.5.*(sigma_ij_2500m_SLder.sigma_xx_slip(i,j)+sigma_ij_2500m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2500m_SLder.sigma_xx_slip(i,j)-sigma_ij_2500m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_2500m_SLder.sigma_xy_slip(i,j).^2)));   
%2750m 
        princ_sigma1_2750m_slip(i,j) = (0.5.*(sigma_ij_2750m_SLder.sigma_xx_slip(i,j)+sigma_ij_2750m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_2750m_SLder.sigma_xx_slip(i,j)-sigma_ij_2750m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_2750m_SLder.sigma_xy_slip(i,j).^2)));             
%3000m 
        princ_sigma1_3000m_slip(i,j) = (0.5.*(sigma_ij_3000m_SLder.sigma_xx_slip(i,j)+sigma_ij_3000m_SLder.sigma_yy_slip(i,j))) + ...
            (sqrt(((0.5.*(sigma_ij_3000m_SLder.sigma_xx_slip(i,j)-sigma_ij_3000m_SLder.sigma_yy_slip(i,j))).^2) + (sigma_ij_3000m_SLder.sigma_xy_slip(i,j).^2)));    

     end
end

%% find r_m: fails when there is no contour at 200 kPa --> 0
for i=1
     % dmax: open + slip
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(9,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(10,i) = 0 % nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: open 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(9,i) =  0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(10,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: slip 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(4,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(5,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(6,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(7,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(9,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(10,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));      
end

for i=2
     % dmax: open + slip
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: open 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: slip 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(4,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(5,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(6,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(7,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(9,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(10,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));      
end

for i=3
     % dmax: open + slip
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: open 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: slip 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(4,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(5,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(6,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(7,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(9,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(10,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));      
end

for i=4
     % dmax: open + slip
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: open 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: slip 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(5,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(6,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(7,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(9,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(10,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));      
end

for i=5
     % dmax: open + slip
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: open 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(6,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(7,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(8,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(9,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(10,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_open(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_open(11,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));     
        
        % dmax: slip 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(1,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(2,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));        
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(3,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(4,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2))) 
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(5,i) = nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_1750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(6,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)))
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(7,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2250m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(8,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2500m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(9,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_2750m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(10,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));
        [C,~]=contour(xx,yy,reshape(((princ_sigma1_3000m_slip(i,:))./1e3),121,121),[v200 v200]); C(C > 15) = NaN;
        idx = C(1,:)<=0; C(1,idx)=NaN;
        dmax_slip(11,i) = 0 %nanmax(sqrt((C(1,2:end).^2)+(C(2,2:end).^2)));      
end

%% rows = volume; columns = thickness
save dmax_idealized_200kPa_20230630.mat dmax
save dmax_open_idealized_200kPa_20230630.mat dmax_open
save dmax_slip_idealized_200kPa_20230630.mat dmax_slip

%% H = 1500 m; V = 0.05 km3; Idealized forward paper figure
v = [-1500:5:1500]; % sigma contours [ kPa ]
v200 = [200 200]; % sigma threshold [ kPa ]
m = 10; % font size 
load BWR.mat % colormap

H = [500,1000,1500,2000,2500,3000]; % Ice thickness range [ m ]
V = [0.001,0.005,0.01,0.05,0.1].*(1e9); % Volume range [ m^3 ]
dmax_x = repmat(H',1,11); % parameter space H
dmax_y = repmat(linspace(500,3000,11),5,1); % parameter space V

% Lai et al [2021] for L1A 2011 for cavity opening height
alpha = 0.32; % [ ] from Lai 2021
hmax2011 = 1.03; % [ m ]
Vi_2011 = 0.008*(1e3^3); % [ m^3 ]
% Volume range
Vi = [0.001, 0.005, 0.01, 0.05, 0.1].*(1e9); % [ m^3 ]
% hmax, R, h over range
hmax = hmax2011.*((Vi./Vi_2011).^(1/3)); % [ m ]
R = sqrt(Vi./(2*pi*alpha.*hmax)); % [ m ]
r = linspace(-15000,15000,800); % [ m ] 
for i=1:800; h_r(i,:) = hmax.*(sqrt(1-((r(i)./R).^2))); end
h_r = real(h_r);

fig1=figure(1); clf; 
set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',[1 1 150 150]./5.5);
tt=0.21; ss = 28/30; scale=0.65;

axe0 = axes('Position',[0.825 0.924 0.16 0.05],'Box','on'); 
axe1 = axes('Position',[0.07 0.78 tt tt*ss],'Box','on','xticklabel',[]); 
axe2 = axes('Position',[0.32 0.78 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);
axe3 = axes('Position',[0.57 0.78 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 

axe12 = axes('Position',[0.07 0.555 tt tt*ss],'Box','on','xticklabel',[]); 
axe22 = axes('Position',[0.32 0.555 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]);
axe32 = axes('Position',[0.57 0.555 tt tt*ss],'Box','on','yticklabel',[],'xticklabel',[]); 

axe13 = axes('Position',[0.07 0.33 tt tt*ss],'Box','on'); 
axe23 = axes('Position',[0.32 0.33 tt tt*ss],'Box','on','yticklabel',[]);
axe33 = axes('Position',[0.57 0.33 tt tt*ss],'Box','on','yticklabel',[]); 

axe5 = axes('Position',[0.07 0.05 tt tt*ss],'Box','on'); 
axe6 = axes('Position',[0.32 0.05 tt tt*ss],'Box','on','yticklabel',[]); 
axe7 = axes('Position',[0.57 0.05 tt tt*ss],'Box','on','yticklabel',[]); 

axes(axe0) % visualize cavity opening
hold on; 
plot(r./1e3,h_r(:,1),'-','LineWidth',1.2,'Color',[217 240 163]./250);
plot(r./1e3,h_r(:,2),'-','LineWidth',1.2,'Color',[173 221 142]./250);
plot(r./1e3,h_r(:,3),'-','LineWidth',1.2,'Color',[120 198 121]./250);
plot(r./1e3,h_r(:,4),'-','LineWidth',1.2,'Color',[49 163 84]./250);
plot(r./1e3,h_r(:,5),'-','LineWidth',1.2,'Color',[0 104 55]./250);
xlabel('x  [ km ]'); ylabel('h  [ m ]');
title('Cavity Opening Height','FontSize',m-1,'FontName','Avenir');
xlim([-15 15]); ylim([0 3]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[0:1:3],'FontSize',m-1,'FontName','Avenir',...
    'LineWidth',0.8);
grid on 

axes(axe1) % case 1 -- open
text(-23.5,14,'a','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1500m_open(4,:))./1e3),121,121),v); 
set(h3,'LineColor','none'); 

% idealized slip patches with opening and/or slip
z = double(gt(patches_idealized.opening(:,4),0)); % patches with opening and/or slip
Nx=40; Ny=40; width = 0.5; % patch width [km]
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
    
    if z(isf)==1
        patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6],'FaceColor','none','LineWidth',0.5);
    else end
  end
 end

 % r_m
[C,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_open(4,:))./1e3),121,121),[v200 v200],'k','LineWidth',1.2);
C(C > 15) = NaN; idx = C(1,:)<=0; C(1,idx)=NaN;
d_plot = sqrt(((C(1,1:end).^2)+(C(2,1:end).^2))); [I,J] = max(d_plot);
quiver([0], [0], [C(1,J)], [C(2,J)],'k','autoscale','off','LineWidth',1.4,'AutoScaleFactor',0.75)
text(7.5, -7.2, 'r_{m}','FontName','Avenir')
quiver([6.7],[-6.9],[-5.8],[+3.1],'LineWidth',1.0,'Color',[0.5 0.5 0.5],'MaxHeadSize',0.7);

ylabel('y  [ km ]'); %xlabel('x  [ km ]');
text(-14, 13.3, 'V=0.05 km^{3}  H = 1500 m  \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
title('\sigma_{1,bedopen}');
caxis([-1000 1000]);
xlim([-15 15]); ylim([-15 15]);
colormap(axe1,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on 

axes(axe2) % case 2 -- slip
text(-18.5,14,'b','FontWeight','bold','FontSize',m+1);
hold on;
[~,h3]=contourf(xx,yy,reshape(((princ_sigma1_1500m_slip(4,:))./1e3),121,121),v);  set(h3,'LineColor','none'); 

% idealized slip patches with opening and/or slip
z = double(gt(patches_idealized.opening(:,4),0)); % patches with opening and/or slip
Nx=40; Ny=40; width = 0.5; % patch width [km]
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
    
    if z(isf)==1
        patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6],'FaceColor','none','LineWidth',0.5);
    else end
  end
 end

% r_m
[C2,~]=contour(xx,yy,reshape(((princ_sigma1_1500m_slip(4,:))./1e3),121,121),[v200 v200],'k','LineWidth',1.2);
% C2(C2 > 15) = NaN;
% d_plot2 = sqrt(((C2(1,1:end).^2)+(C2(2,1:end).^2)));
% [I,J] = max(d_plot2); idx = C2(1,:)<=0; C2(1,idx)=NaN;
% quiver([0], [0], [C2(1,J)], [C2(2,J)],'k','autoscale','off','LineWidth',1.4,'AutoScaleFactor',0.75)
% text(10, 5, 'r_{m}','FontName','Avenir')
% quiver([9.7],[4.3],[-3.5],[-3.3],'LineWidth',1.0,'Color',[0.5 0.5 0.5],'MaxHeadSize',0.7);

caxis([-1000 1000]);
text(-14, 13.3, 'V=0.05 km^{3}  H = 1500 m  \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-15 15]); ylim([-15 15]);
colormap(axe2,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('\sigma_{1,bedslip}'); 

axes(axe3) % case 3 -- open and slip
text(-18.5,14,'c','FontWeight','bold','FontSize',m+1);
hold on;
[c3,h3]=contourf(xx,yy,reshape(((princ_sigma1_1500m(4,:))./1e3),121,121),v); set(h3,'LineColor','none'); 

% idealized slip patches with opening and/or slip
z = double(gt(patches_idealized.opening(:,4),0)); % patches with opening and/or slip
Nx=40; Ny=40; width = 0.5; % patch width [km]
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
    
    if z(isf)==1
        patch(x,y,z(isf),'EdgeColor',[0.6 0.6 0.6],'FaceColor','none','LineWidth',0.5);
    else end
  end
 end

% r_m
[C3,h]=contour(xx,yy,reshape(((princ_sigma1_1500m(4,:))./1e3),121,121),[v200 v200],'k','LineWidth',1.2);
C3(C3 > 15) = NaN; idx = C3(1,:)<=0; C3(1,idx)=NaN;
d_plot = sqrt(((C3(1,1:end).^2)+(C3(2,1:end).^2))); [I,J] = max(d_plot);
quiver([0], [0], [C3(1,J)], [C3(2,J)],'k','autoscale','off','LineWidth',1.4,'AutoScaleFactor',0.75)
text(7.6, 7.2, 'r_{m}','FontName','Avenir')
quiver([7.2],[7.2],[-5.9],[-3.2],'LineWidth',1.0,'Color',[0.5 0.5 0.5],'MaxHeadSize',0.7);

caxis([-1000 1000]);
text(-14, 13.3, 'V=0.05 km^{3}  H = 1500 m  \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-15 15]); ylim([-15 15]);
colormap(axe3,BWR); 
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-15:5:15],'FontSize',m,'FontName','Avenir',...
    'LineWidth',0.8,'Layer','top');
grid on
title('\sigma_{1,bedopen+bedslip}'); %xlabel('x  [ km ]');

cc43 = colorbar('EastOutside'); 
set(cc43,'Position',[0.789 0.785 .005 (tt*ss)-0.11],'xtick',[-1000 -500 0 500 1000],'tickdir','in');
text(23.75,-11.0, '\sigma_{1} [ kPa ]','Rotation',90,'FontSize',m,'FontName','Avenir');


axes(axe12) % range over lake volume
text(-23.5,1300,'d','FontWeight','bold','FontSize',m+1);
hold on;
for i=1:5
    this = reshape(princ_sigma1_1500m_open(i,:),121,121);
    transect_open(i,:) = this(62,:); % transect at y = 0 m 
end

plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);
plot(xx(42,:), transect_open(5,:)./1e3,'LineWidth',1.2,'Color',[0 104 55]./250); % 0.001
plot(xx(42,:), transect_open(4,:)./1e3,'LineWidth',1.2,'Color',[49 163 84]./250);
plot(xx(42,:), transect_open(3,:)./1e3,'LineWidth',1.2,'Color',[120 198 121]./250);
plot(xx(42,:), transect_open(2,:)./1e3,'LineWidth',1.2,'Color',[173 221 142]./250);
plot(xx(42,:), transect_open(1,:)./1e3,'LineWidth',1.2,'Color',[217 240 163]./250); % 0.1

ylabel('\sigma_{1} [ kPa ]');
text(-14, 1280, 'H = 1500 m  \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
text(-14, 275, '200 kPa','FontName','Avenir')
xlim([-15 15]); ylim([-200 1400]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-200:200:2000],'FontSize',m,'FontName','Avenir','LineWidth',0.8);
grid on 

axes(axe22)
text(-18.5,1300,'e','FontWeight','bold','FontSize',m+1);
hold on;
for i=1:5
    this = reshape(princ_sigma1_1500m_slip(i,:),121,121);
    transect_slip(i,:) = this(62,:); % transect at y = 0 m 
end

plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);
plot(xx(42,:), transect_slip(5,:)./1e3,'LineWidth',1.2,'Color',[0 104 55]./250); % 0.001
plot(xx(42,:), transect_slip(4,:)./1e3,'LineWidth',1.2,'Color',[49 163 84]./250);
plot(xx(42,:), transect_slip(3,:)./1e3,'LineWidth',1.2,'Color',[120 198 121]./250);
plot(xx(42,:), transect_slip(2,:)./1e3,'LineWidth',1.2,'Color',[173 221 142]./250);
plot(xx(42,:), transect_slip(1,:)./1e3,'LineWidth',1.2,'Color',[217 240 163]./250); % 0.1

text(-14, 1280, 'H = 1500 m  \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
text(-14, 275, '200 kPa','FontName','Avenir')
xlim([-15 15]); ylim([-200 1400]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-200:200:2000],'FontSize',m,'FontName','Avenir','LineWidth',0.8);
grid on

axes(axe32)
text(-18,1300,'f','FontWeight','bold','FontSize',m+1);
hold on;
for i=1:5
    this = reshape(princ_sigma1_1500m(i,:),121,121);
    transect(i,:) = this(62,:); % transect at y = 0 m
end

plot(xx(42,:), transect(5,:)./1e3,'LineWidth',1.2,'Color',[0 104 55]./250); % 0.001
plot(xx(42,:), transect(4,:)./1e3,'LineWidth',1.2,'Color',[49 163 84]./250);
plot(xx(42,:), transect(3,:)./1e3,'LineWidth',1.2,'Color',[120 198 121]./250);
plot(xx(42,:), transect(2,:)./1e3,'LineWidth',1.2,'Color',[173 221 142]./250);
plot(xx(42,:), transect(1,:)./1e3,'LineWidth',1.2,'Color',[217 240 163]./250); % 0.1
plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);

plot(xx(42,:), transect(5,:)./1e3,'LineWidth',1.2,'Color',[0 104 55]./250); % 0.001
plot(xx(42,:), transect(4,:)./1e3,'LineWidth',1.2,'Color',[49 163 84]./250);
plot(xx(42,:), transect(3,:)./1e3,'LineWidth',1.2,'Color',[120 198 121]./250);
plot(xx(42,:), transect(2,:)./1e3,'LineWidth',1.2,'Color',[173 221 142]./250);
plot(xx(42,:), transect(1,:)./1e3,'LineWidth',1.2,'Color',[217 240 163]./250); % 0.1

text(-14, 1280, 'H = 1500 m  \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
text(-14, 275, '200 kPa','FontName','Avenir')
xlim([-15 15]); ylim([-200 1400]);
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-200:200:2000],'FontSize',m,'FontName','Avenir','LineWidth',0.8);
grid on 

lgd = legend('0.1 km^{3}','0.05 km^{3} (Panels a-c)','0.01 km^{3}','0.005 km^{3}','0.001 km^{3}','200 kPa');
legend boxoff;
set(lgd,'Position',[0.83 0.56 0.1 0.1],'FontSize',m-1);
text(18,825,'Lake Volume:','FontSize',m,'FontName','Avenir')


axes(axe13) % range over ice thickness
text(-23.5,1400,'g','FontWeight','bold','FontSize',m+1);
hold on;
for i=4
    this051 = reshape(princ_sigma1_500m_open(i,:),121,121);
    this1 = reshape(princ_sigma1_1000m_open(i,:),121,121);
    this15 = reshape(princ_sigma1_1500m_open(i,:),121,121);
    this2 = reshape(princ_sigma1_2000m_open(i,:),121,121);
    this25 = reshape(princ_sigma1_2500m_open(i,:),121,121);
    this3 = reshape(princ_sigma1_3000m_open(i,:),121,121);
    transect_thick_open = vertcat(this1(62,:),this15(62,:),this2(62,:),this25(62,:),this3(62,:)); % transect at y = 0 m 
end

plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);
plot(xx(52,:), this051(62,:)./1e3,'LineWidth',1.2,'Color',[8 81 156]./250); % 500
plot(xx(42,:), transect_thick_open(1,:)./1e3,'LineWidth',1.2,'Color',[49 130 189]./250);
plot(xx(42,:), transect_thick_open(2,:)./1e3,'LineWidth',1.2,'Color',[107 174 214]./250);
plot(xx(42,:), transect_thick_open(3,:)./1e3,'LineWidth',1.2,'Color',[158 202 225]./250);
plot(xx(42,:), transect_thick_open(4,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); 
plot(xx(42,:), transect_thick_open(5,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); % 3000

% plot bounds of blister in x
bound_pos = max(r(h_r(:,4)~=0)); % maximum x value for non-zero h_r
bound_neg = min(r(h_r(:,4)~=0)); % minimum x value for non-zero h_r
ha = area([bound_pos bound_neg]./1e3, [2000 2000],'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.2);
ha2 = area([bound_pos bound_neg]./1e3, [-2000 -2000],'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.2); 

plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);
plot(xx(52,:), this051(62,:)./1e3,'LineWidth',1.2,'Color',[8 81 156]./250); % 500
plot(xx(42,:), transect_thick_open(1,:)./1e3,'LineWidth',1.2,'Color',[49 130 189]./250);
plot(xx(42,:), transect_thick_open(2,:)./1e3,'LineWidth',1.2,'Color',[107 174 214]./250);
plot(xx(42,:), transect_thick_open(3,:)./1e3,'LineWidth',1.2,'Color',[158 202 225]./250);
plot(xx(42,:), transect_thick_open(4,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); 
plot(xx(42,:), transect_thick_open(5,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); % 3000

text(-4, 1600, 'Blister Width','FontSize',m-2,'FontName','Avenir')
plot([bound_neg bound_pos]./1e3,[1450 1450],'k.-','Linewidth',1.2,'MarkerSize',10)

ylabel('\sigma_{1} [ kPa ]','FontName','Avenir');
text(-14, 275, '200 kPa','FontName','Avenir')
text(-14, -350, 'V=0.05 km^{3}                      \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
xlim([-15 15]); ylim([-500 1500]); xlabel('x  [ km ]');
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-500:500:2000],'FontSize',m,'FontName','Avenir','LineWidth',0.8);
grid on 

axes(axe23)
text(-18.5,1400,'h','FontWeight','bold','FontSize',m+1);
hold on;
for i=4
    this052 = reshape(princ_sigma1_500m_slip(i,:),121,121);
    this1 = reshape(princ_sigma1_1000m_slip(i,:),121,121);
    this15 = reshape(princ_sigma1_1500m_slip(i,:),121,121);
    this2 = reshape(princ_sigma1_2000m_slip(i,:),121,121);
    this25 = reshape(princ_sigma1_2500m_slip(i,:),121,121);
    this3 = reshape(princ_sigma1_3000m_slip(i,:),121,121);
    transect_thick_slip = vertcat(this1(62,:),this15(62,:),this2(62,:),this25(62,:),this3(62,:)); % transect at y = 0 m 
end

plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);
plot(xx(52,:), this052(62,:)./1e3,'LineWidth',1.2,'Color',[8 81 156]./250);  % 500
plot(xx(42,:), transect_thick_slip(1,:)./1e3,'LineWidth',1.2,'Color',[49 130 189]./250);
plot(xx(42,:), transect_thick_slip(2,:)./1e3,'LineWidth',1.2,'Color',[107 174 214]./250);
plot(xx(42,:), transect_thick_slip(3,:)./1e3,'LineWidth',1.2,'Color',[158 202 225]./250);
plot(xx(42,:), transect_thick_slip(4,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); 
plot(xx(42,:), transect_thick_slip(5,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); % 3000

% plot bounds of blister in x
ha = area([bound_pos bound_neg]./1e3, [2000 2000],'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.2);
ha2 = area([bound_pos bound_neg]./1e3, [-2000 -2000],'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.2);

plot([-20 20],[200 200],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);
plot(xx(52,:), this052(62,:)./1e3,'LineWidth',1.2,'Color',[8 81 156]./250);  % 500
plot(xx(42,:), transect_thick_slip(1,:)./1e3,'LineWidth',1.2,'Color',[49 130 189]./250);
plot(xx(42,:), transect_thick_slip(2,:)./1e3,'LineWidth',1.2,'Color',[107 174 214]./250);
plot(xx(42,:), transect_thick_slip(3,:)./1e3,'LineWidth',1.2,'Color',[158 202 225]./250);
plot(xx(42,:), transect_thick_slip(4,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); 
plot(xx(42,:), transect_thick_slip(5,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); % 3000

text(-5, 1600, 'Slip Patch Width','FontSize',m-2,'FontName','Avenir')
plot([bound_neg bound_pos]./1e3,[1450 1450],'k.-','Linewidth',1.2,'MarkerSize',10) 

text(-14, -350, 'V=0.05 km^{3}                      \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
text(-14, 275, '200 kPa','FontName','Avenir')
xlim([-15 15]); ylim([-500 1500]); xlabel('x  [ km ]');
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-500:500:2000],'FontSize',m,'FontName','Avenir','LineWidth',0.8);
grid on 

axes(axe33)
text(-18,1400,'i','FontWeight','bold','FontSize',m+1);
hold on;
for i=4
    this053 = reshape(princ_sigma1_500m(i,:),121,121);
    this1 = reshape(princ_sigma1_1000m(i,:),121,121);
    this15 = reshape(princ_sigma1_1500m(i,:),121,121);
    this2 = reshape(princ_sigma1_2000m(i,:),121,121);
    this25 = reshape(princ_sigma1_2500m(i,:),121,121);
    this3 = reshape(princ_sigma1_3000m(i,:),121,121);
    transect_thick = vertcat(this1(62,:),this15(62,:),this2(62,:),this25(62,:),this3(62,:)); % transect at y = 0 m 
end

plot(xx(52,:), this053(62,:)./1e3,'LineWidth',1.2,'Color',[8 81 156]./250);  % 500
plot(xx(42,:), transect_thick(1,:)./1e3,'LineWidth',1.2,'Color',[49 130 189]./250);
plot(xx(42,:), transect_thick(2,:)./1e3,'LineWidth',1.2,'Color',[107 174 214]./250);
plot(xx(42,:), transect_thick(3,:)./1e3,'LineWidth',1.2,'Color',[158 202 225]./250);
plot(xx(42,:), transect_thick(4,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); 
plot(xx(42,:), transect_thick(5,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250);  % 3000
plot([-20 20],[150 150],'-','LineWidth',1.2,'Color',[0.6 0.6 0.6]);

% plot bounds of blister in x
ha = area([bound_pos bound_neg]./1e3, [2000 2000],'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.2);
ha2 = area([bound_pos bound_neg]./1e3, [-2000 -2000],'FaceColor',[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.2);

plot(xx(52,:), this053(62,:)./1e3,'LineWidth',1.2,'Color',[8 81 156]./250);  % 500
plot(xx(42,:), transect_thick(1,:)./1e3,'LineWidth',1.2,'Color',[49 130 189]./250);
plot(xx(42,:), transect_thick(2,:)./1e3,'LineWidth',1.2,'Color',[107 174 214]./250);
plot(xx(42,:), transect_thick(3,:)./1e3,'LineWidth',1.2,'Color',[158 202 225]./250);
plot(xx(42,:), transect_thick(4,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); 
plot(xx(42,:), transect_thick(5,:)./1e3,'LineWidth',1.2,'Color',[198 219 239]./250); % 3000

text(-9, 1600, 'Blister and Slip Patch Width','FontSize',m-2,'FontName','Avenir')
plot([bound_neg bound_pos]./1e3,[1450 1450],'k.-','Linewidth',1.2,'MarkerSize',10)

text(-14, -350, 'V=0.05 km^{3}                      \mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
text(-14, 275, '200 kPa','FontName','Avenir')
xlim([-15 15]); ylim([-500 1500]); xlabel('x  [ km ]');
set(gca,'tickdir','in','xtick',[-15:5:15],'ytick',[-500:500:2000],'FontSize',m,'FontName','Avenir','LineWidth',0.8);
grid on 

lgd = legend('500 m','1000 m','1500 m (Panels a-c)','2000 m','2500 m','3000 m','200 kPa');
legend boxoff;
set(lgd,'Position',[0.825 0.335 0.1 0.1],'FontSize',m-1);
text(18,800,'Ice Thickness:','FontSize',m,'FontName','Avenir')

axes(axe5) % parameter space r_m
text(-150,7.5,'j','FontWeight','bold','FontSize',m+1);
hold on;
title('\sigma_{1,bedopen}');

plot(dmax_y(5,:),dmax_open(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax_open(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax_open(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax_open(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax_open(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot([0 5000],[R(5) R(5)]./1e3,':','LineWidth',1.0,'Color',[0 104 55]./250);
plot([0 5000],[R(4) R(4)]./1e3,':','LineWidth',1.0,'Color',[49 163 84]./250);
plot([0 5000],[R(3) R(3)]./1e3,':','LineWidth',1.0,'Color',[120 198 121]./250);
plot([0 5000],[R(2) R(2)]./1e3,':','LineWidth',1.0,'Color',[173 221 142]./250);
plot([0 5000],[R(1) R(1)]./1e3,':','LineWidth',1.0,'Color',[217 240 163]./250);

plot(dmax_y(5,:),dmax_open(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax_open(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax_open(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax_open(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax_open(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot(dmax_y(1,5),dmax_open(5,4),'p','LineWidth',1.2,'Color','k',...
    'MarkerFaceColor','none','MarkerSize',8.5);

xlim([500 3000]); ylim([0 8]);
set(gca,'xtick',[500:500:3000],'ytick',[0:2:8],'tickdir','in','FontSize',m,'FontName','Avenir','LineWidth',0.8)
ylabel('r_{m}(\sigma_{1}>200 kPa)  [ km ]'); 
xlabel('Ice Thickness, H  [ m ]');
text(2100, 7.25, '\mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
grid on 

axes(axe6)
text(220,7.5,'k','FontWeight','bold','FontSize',m+1);
hold on;
title('\sigma_{1,bedslip}');

plot(dmax_y(5,:),dmax_slip(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax_slip(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax_slip(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax_slip(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax_slip(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot([0 5000],[R(5) R(5)]./1e3,':','LineWidth',1.0,'Color',[0 104 55]./250);
plot([0 5000],[R(4) R(4)]./1e3,':','LineWidth',1.0,'Color',[49 163 84]./250);
plot([0 5000],[R(3) R(3)]./1e3,':','LineWidth',1.0,'Color',[120 198 121]./250);
plot([0 5000],[R(2) R(2)]./1e3,':','LineWidth',1.0,'Color',[173 221 142]./250);
plot([0 5000],[R(1) R(1)]./1e3,':','LineWidth',1.0,'Color',[217 240 163]./250);

plot(dmax_y(5,:),dmax_slip(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax_slip(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax_slip(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax_slip(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax_slip(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot(dmax_y(1,5),dmax_slip(5,4),'p','LineWidth',1.2,'Color','k',...
    'MarkerFaceColor','none','MarkerSize',8.5);

xlim([500 3000]); ylim([0 8]);
set(gca,'xtick',[500:500:3000],'ytick',[0:2:8],'tickdir','in','FontSize',m,'FontName','Avenir','LineWidth',0.8)
xlabel('Ice Thickness, H  [ m ]');
text(2100, 7.25, '\mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
grid on

axes(axe7)
text(220,7.5,'l','FontWeight','bold','FontSize',m+1);
hold on;
title('\sigma_{1,bedopen+bedslip}');
plot(dmax_y(5,:),dmax(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot(dmax_y(1,5),dmax(5,4),'p','LineWidth',1.2,'Color','k',...
    'MarkerFaceColor','none','MarkerSize',8.5);

% L1A 2011 and 2012 (L1A center is at [1.5,0.5] in map space)
load NLdmax2011_200kPa_20230621.mat
load NLdmax2012_200kPa_20230621.mat

% dock at eel pond 
sky_blue = [111, 169, 228]./255; % L1B
metal = [87, 115, 131]./255; % L1C
oar = [251, 219, 154]./255; % L1D
handle = [161, 37, 49]./255; % L1A
dark_oar = [164, 114, 63]./255;

plot(980,NLdmax2011.mu15_inland-1.5,...
    'sb','LineWidth',0.9,'MarkerFaceColor',handle,'Color',handle,'MarkerSize',5);
plot(980,NLdmax2012.mu15_inland-1.5,...
    'sr','LineWidth',0.9,'MarkerFaceColor',dark_oar,'Color',dark_oar,'MarkerSize',5);

plot(dmax_y(5,:),dmax(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot([0 5000],[R(5) R(5)]./1e3,':','LineWidth',1.0,'Color',[0 104 55]./250);
plot([0 5000],[R(4) R(4)]./1e3,':','LineWidth',1.0,'Color',[49 163 84]./250);
plot([0 5000],[R(3) R(3)]./1e3,':','LineWidth',1.0,'Color',[120 198 121]./250);
plot([0 5000],[R(2) R(2)]./1e3,':','LineWidth',1.0,'Color',[173 221 142]./250);
plot([0 5000],[R(1) R(1)]./1e3,':','LineWidth',1.0,'Color',[217 240 163]./250);

plot(dmax_y(5,:),dmax(:,5),'-o','LineWidth',1.0,'Color',[0 104 55]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',6);
plot(dmax_y(4,:),dmax(:,4),'-o','LineWidth',1.0,'Color',[49 163 84]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5.5);
plot(dmax_y(3,:),dmax(:,3),'-o','LineWidth',1.0,'Color',[120 198 121]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',5);
plot(dmax_y(2,:),dmax(:,2),'-o','LineWidth',1.0,'Color',[173 221 142]./250,...FBlister
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',4.5);
plot(dmax_y(1,:),dmax(:,1),'-o','LineWidth',1.0,'Color',[217 240 163]./250,...
    'MarkerFaceColor',[0.4 0.4 0.4],'MarkerSize',3.5);

plot(980,NLdmax2011.mu15_inland-1.5,...
    'sb','LineWidth',0.9,'MarkerFaceColor',handle,'Color',handle,'MarkerSize',5);
plot(980,NLdmax2012.mu15_inland-1.5,...
    'sr','LineWidth',0.9,'MarkerFaceColor',dark_oar,'Color',dark_oar,'MarkerSize',5);

plot(dmax_y(1,5),dmax(5,4),'p','LineWidth',1.2,'Color','k',...
    'MarkerFaceColor','none','MarkerSize',8.5);

xlim([500 3000]); ylim([0 8]);
set(gca,'xtick',[500:500:3000],'ytick',[0:2:8],'tickdir','in','FontSize',m,'FontName','Avenir','LineWidth',0.8);
xlabel('Ice Thickness, H  [ m ]');
text(2100, 7.25, '\mu = 1.5 GPa','FontName','Avenir','FontSize',m-1)
grid on

lgd = legend('0.1 km^{3}','0.05 km^{3}','0.01 km^{3}','0.005 km^{3}','0.001 km^{3}',...
    'Cases in Panels a-c','L1A 2011 (V=0.008 km^{3})','L1A 2012 (V=0.008 km^{3})');
legend boxoff;
set(lgd,'Position',[0.83 0.075 0.1 0.1],'FontSize',m-1);
text(3200,6.75,'Lake Volume:','FontSize',m,'FontName','Avenir')

%% print figure
% figurename=sprintf('paperfig8_idealized_15GPa_20230630.png');
% print(gcf,'-dpng','-r600',figurename);
