%% RUNFILE FOR NIF IDEALIZED BLISTER RANGING OVER ICE THICKNESS
%  LAS November 20 2020 -- run the 2011 NIF with a 20km by 20 km basal patch.
%  LAS 2023 June 30 -- Use with makegeom_20km_bed.m (simplified geometry
%  that does not include computations for a vertical crack)
%  LAS 2023 July 05 -- include E-W vertical HF opening of 0.5 m 

path(path,'../software/general')
path(path,'../software/functions/dataIO')
path(path,'../software/objects')
path(path,'../software/functions/elmat')
path(path,'../software/functions/time')
path(path,'../software/filter_tools')
path(path,'../software/matools')
path(path,'../software/okada85')

%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP FAULT GEOMETRY%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
origin=[68.72, -49.53];  %% SET COORDINATE SYSTEM ORIGIN IN MIDDLE OF ARRAY
 
% Green function positions on half-space surface
[xx,yy] = meshgrid(-15:0.25:15, -15:0.25:15); % 250-m grid over 30 km by 30 km region
xy_surf = horzcat(reshape(xx,121*121,1), reshape(yy,121*121,1));
Nsurface=length(xy_surf);

depthc = [0.50:0.25:3.00]; % [ km ] DEPTH OF BEDPLANE

% output G matrices strains  1000m -- include E-W vertical crack
[patchesC, patchesB, G3v, G2bv, G3bv, G3v_strain, G2bv_strain, G3bv_strain] = makegeom_20km_bed_EWcurve(origin,depthc(3));
Gsurface.G3v = G3v; Gsurface.G3v_strain = G3v_strain;
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesC = patchesC; Gsurface.patchesB = patchesB; 
Gsurface.Nsurface = Nsurface; 
save('Gsurface1000m_strain_EWcurve.mat', 'Gsurface', '-v7.3'); 

% output G matrices strains  1250m -- include E-W vertical crack
[patchesC, patchesB, G3v, G2bv, G3bv, G3v_strain, G2bv_strain, G3bv_strain] = makegeom_20km_bed_EWcurve(origin,depthc(4));
Gsurface.G3v = G3v; Gsurface.G3v_strain = G3v_strain;
Gsurface.G2bv = G2bv; Gsurface.G3bv = G3bv;
Gsurface.G2bv_strain = G2bv_strain; Gsurface.G3bv_strain = G3bv_strain;
Gsurface.xy_surf = xy_surf; Gsurface.xx = xx; Gsurface.yy = yy;
Gsurface.patchesC = patchesC; Gsurface.patchesB = patchesB; 
Gsurface.Nsurface = Nsurface; 
save('Gsurface1250m_strain_EWcurve.mat', 'Gsurface', '-v7.3'); 

%%
disp('Done saving output')       
