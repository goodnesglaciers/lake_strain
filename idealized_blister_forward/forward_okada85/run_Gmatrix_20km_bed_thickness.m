%% RUNFILE FOR NIF IDEALIZED BLISTER RANGING OVER ICE THICKNESS
%  LAS November 20 2020 -- run the 2011 NIF with a 20km by 20 km basal patch.
%  LAS 2023 June 30 -- Use with makegeom_20km_bed.m (simplified geometry
%  that does not include computations for a vertical crack)

path(path,'../software/general')
path(path,'../software/functions/dataIO')
path(path,'../software/objects')
path(path,'../software/functions/elmat')
path(path,'../software/functions/time')
path(path,'../software/filter_tools')
path(path,'../software/matools')
path(path,'../software/okada85')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load a .mat file that has tsb and tsb_enu in them which are time series
% objects from readsinex.m that are in GPS and ENU coordinates respectively
% preferably already downsampled in time to make this go quicker.

         %load ts_displ_week_4stations.mat % just four stations
         load ts_displ_week_NLBS.mat
         ts_displ = ts_displ_week_NLBS;
         
         tstart=167.1;
         tend=172.0;
         [a,I1]=min(abs(ts_displ.epochs-tstart));
         [a,I2]=min(abs(ts_displ.epochs-tend));
         
% Select epochs by subsetting 
% May 24 2018: Changed subsetting from '1' to '3' to downsample timeseries object
         ts2=ts_displ(:,I1:3:I2); % change to 1 for MLE runs
         
         figure
         for i=1:length(ts2.epochs)-1
             plot(i,length(ts2(:,i).d),'*')
             hold on;
         end
         xlabel('epoch #')
         ylabel('number of observations')
         
%%%%%%%% SELECT STATIONS %%%%%%%%%%%%%%%%%%%%%%%
tsb_enu=ts2;   %can't do anything here until you get lats+lons in here instead of hardwired.
load apcoords_lle.mat; lons=apcoords_lle(2,:); lats=apcoords_lle(1,:);
Nsites=length(lats);
plotyn=1;
        if(plotyn)
         figure
         plot(lons,lats,'ko','MarkerFaceColor','r','MarkerSize',4)
         hold on;
         for j=1:length(lats)
           text(lons(j),lats(j),tsb_enu.sites(j));
         end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PICK EPOCHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        J=find(tsb_enu.epochs>=2011.3 & tsb_enu.epochs<=2011.6);
%        J=J(1:2:end);       
%        ts_enu=tsb_enu(:,J); ts=tsb(:,J);  %subsetting

ts_enu=tsb_enu;   %already subsetted above

        % uncomment below if you think you have some bad solutions to throw
        % out a priori; fix threshold to suit
        %disp('finding strange epochs')
	    %for i=1:length(ts_enu.epochs);
    	%	[U,S,V]=svd(ts_enu(:,i).dcov);
     	%	mineig(i)=min(diag(S));
	    %end
        %KK=find(mineig>=1e-10);
        %ts_enu=ts_enu(:,KK);
        %disp('done')

%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP FAULT GEOMETRY%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
dosetup=1;
if(dosetup)
  origin=[68.72, -49.53];  %% SET COORDINATE SYSTEM ORIGIN IN MIDDLE OF ARRAY
 
% greens function solution positions on half-space surface
[xx,yy] = meshgrid(-15:0.25:15, -15:0.25:15); % 250-m grid over 30 km by 30 km region
xy_surf = horzcat(reshape(xx,121*121,1), reshape(yy,121*121,1));
Nsurface=length(xy_surf);

depthc = [0.50:0.25:3.00]; % [ km ] DEPTH OF BEDPLANE

% % output G matrices strains  500m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(1));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface500m_strain.mat', 'Gsurface', '-v7.3'); 

% % output G matrices strains  750m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(2));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface750m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  1000m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(3));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface1000m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  1250m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(origin,depthc(4));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface1250m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  1500m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(5));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface1500m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  1750m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(6));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface1750m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  2000m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(7));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface2000m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  2250m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(8));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface2250m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  2500m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(9));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface2500m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  2750m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(10));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface2750m_strain.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  3000m
[patchesB, G2bv, G3bv, G2bv_strain, G3bv_strain]=makegeom_20km_bed(ts_enu,origin,depthc(11));
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesB = patchesB; Gsurface.Nsurface = Nsurface; 
save('Gsurface3000m_strain.mat', 'Gsurface', '-v7.3'); 

end
%%
disp('Done saving output')       
