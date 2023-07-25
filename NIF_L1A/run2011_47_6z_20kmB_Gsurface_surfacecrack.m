%% RUNFILE FOR 2011 NIF
%  March 22 2018  run the 2011 NIF with only 4 stations - LAS
%  November 20 2020: run the 2011 NIF with all stations and a 20km by 20 km basal patch.
%  Saves Green function matrices for surface as: Gsurface500_strain.mat
%  Use with makegeom_caseC_20km_strains.m for 20km basal plane + L1A hydro-fracture

path(path,'software/general')
path(path,'software/functions/dataIO')
path(path,'software/objects')
path(path,'software/functions/elmat')
path(path,'software/functions/time')
path(path,'software/filter_tools')
path(path,'software/matools')
path(path,'software/okada85')

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
         ts2=ts_displ(:,I1:1:I2); % change to 1 for MLE runs
         
         figure
         for i=1:length(ts2.epochs)-1
             plot(i,length(ts2(:,i).d),'*')
             hold on;
         end
         xlabel('epoch #')
         ylabel('number of observations')
         
%%%%%%%% SELECT STATIONS %%%%%%%%%%%%%%%%%%%%%%%
tsb_enu=ts2;   %can't do anything here until you get lats+lons in here instead of hardwired.

load apcoords_lle.mat;
lons=apcoords_lle(2,:);
lats=apcoords_lle(1,:);
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

% output strains  
[~, ~, ~, ~, ~, Rc, Rb, patchesC, patchesB, ...
    ~, ~, G3v, G2bv, G3bv, G3v_strain, G2bv_strain, G3bv_strain]=...
    makegeom_caseC_20km_strains_surfacecrack(ts_enu,origin);

end

% save Green function matrices for surface
[xx,yy] = meshgrid(-20:0.5:20, -20:0.5:20); % 500-m grid
xy_surf = horzcat(reshape(xx,81*81,1), reshape(yy,81*81,1));
Nsurface=length(xy_surf);
Gsurface500.G3v = G3v; % vertical crack opening
Gsurface500.G2bv = G2bv; % bed crack dip slip
Gsurface500.G3bv = G3bv; % bed crack opening
Gsurface500.Nsurface = Nsurface;
Gsurface500.xy_surf = xy_surf;
Gsurface500.G3v_strain = G3v_strain;
Gsurface500.G2bv_strain = G2bv_strain;
Gsurface500.G3bv_strain = G3bv_strain;
save Gsurface500_strain_surfacecrack.mat Gsurface500;
  
% % % % %
disp('Done saving output')       
