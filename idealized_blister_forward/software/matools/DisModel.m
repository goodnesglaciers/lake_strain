% figure; fig = gcf; set(fig,'color',[0,0,0.8])
fig = figure('Color',[0,0,0.8],'Name','DisModelv')

%see if a FaultGeom file exists and open it
fid = fopen('FaultGeom');

if fid == -1 
	nf = 1;
	geom = zeros(1,10);
else
	fclose(fid);
	load FaultGeom -ascii
	geom = FaultGeom;
	[nf,nel] = size(geom);
end

%see if a VolGeom file exists and open it
fid = fopen('VolGeom');

if fid == -1 
	nv = 0;
	vgeom = [];
else
	fclose(fid);
	load VolGeom -ascii
	vgeom = VolGeom;
	[nv,nvel] = size(vgeom);
end

% Specify the Number of Faults and Volume Sources
	txt1   = uicontrol(fig, 'Style', 'text',  'String', '# of Faults', 'Position', [5 175 110 25]);
	nfaults = uicontrol(fig, 'Style', 'edit',  'String', num2str(nf), 'Position', [120 175 50 25]);

	txt1   = uicontrol(fig, 'Style', 'text',  'String', '# of Vol Sources', 'Position', [180 175 110 25]);
	nvols  = uicontrol(fig, 'Style', 'edit',  'String', num2str(nv), 'Position', [295 175 50 25]);


% DATA FILE NAMES
frm1 = uicontrol(fig, 'Style', 'frame',  'Position',[10 5 210 165]);

cb_gpsdat = uicontrol(gcf, 'Style', 'checkbox',  'String','GPS Data Filename', 'Position',[15 35 200 25]);
getname_gps   = uicontrol(fig, 'Style', 'edit',  'String', 'Data.vel',  'Position', [15 10 200 25]);

cb_levdat = uicontrol(fig, 'Style', 'checkbox',  'String','Level Data Filename', 'Position',[15 85 200 25]);
getname_lev   = uicontrol(fig, 'Style', 'edit',  'String', 'Level.dat',  'Position', [15 60 200 25]);

cb_edmdat = uicontrol(fig, 'Style', 'checkbox',  'String','EDM Data File', 'Position',[15 135 200 25]);
getname_edm   = uicontrol(fig, 'Style', 'edit',  'String', 'EDM.dat',  'Position', [15 110 200 25]);


%SETUP GEOMETRY
setp =  uicontrol(fig, 'Style', 'push',  'String', 'Setup Geometry', 'Position', [430 125 125 25], 'CallBack','SetupGeom', 'BackgroundColor','white', 'ForegroundColor','black');

%LOAD DATA
distart =  uicontrol(fig, 'Style', 'push',  'String', 'Load Data', 'Position', [430 100 125 25], 'CallBack','LoadData', 'BackgroundColor','white', 'ForegroundColor','black');

% PLOT FILE NAME: CAN'T CHANGE NOW
frm2 = uicontrol(fig, 'Style', 'frame',  'Position',[295 10 115 80]);
filelabel2 = uicontrol(fig, 'Style', 'text',  'String','Plot Filename', 'Position',[300 50 100 25]);
getname2  = uicontrol(fig, 'Style', 'text',  'String', 'test.map',  'Position', [300 15 100 25]);


%Plot GMT file?
cb_plt = uicontrol(gcf, 'Style', 'checkbox',  'String', 'Plot', ...
	'Position', [295 95 115 25]);

origin = [];
