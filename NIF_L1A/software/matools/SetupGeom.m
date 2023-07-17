nf = eval( get(nfaults,'string') );
nv = eval( get(nvols,'string') );

%test whether to resize GEOM:
        [nf_geom,nel] = size(geom);
	[nv_geom,nvel] = size(vgeom);
	if nf < nf_geom 
		frm1 = uicontrol(fig,'Style','frame','BackgroundColor',[0,0,0.8], ...
			'Position',[5  207  15+75*(nf_geom+max(nv,nv_geom)+2) 204]);
%			'Position',[5  215  15+75*(nf_geom+1) 185]);
			geom = geom(1:nf,:);
		end
		if nv > 0
		frm2 = uicontrol(fig, 'Style', 'frame', 'BackgroundColor',...
		[0,0,0.8],'Position',[25+75*(nf_geom+1)  312  10+75*(nv_geom+1) 99]);
		end
	if nf > nf_geom 
		geom = [geom;zeros(nf-nf_geom,10)];
	end

%test whether to resize VGEOM:
if nv>0  | nv_geom >0      
	if nv < nv_geom 
	frm2 = uicontrol(fig, 'Style', 'frame', 'BackgroundColor',...
		[0,0,0.8],'Position',[25+75*(nf_geom+1)  312  10+75*(nv_geom+1) 99]);
		vgeom = vgeom(1:nv,:);
	end
	if nv > nv_geom 
		vgeom = [vgeom;zeros(nv-nv_geom,4)];
	end
end

% Draw Frame and Boxes For Dislocation Elements
frm1 = uicontrol(fig, 'Style', 'frame',  'Position',[5  207  15+75*(nf+1) 204]);

%multilabels = uicontrol(fig, 'Style', 'edit',  'String','Fault#|Latit1|Long1|Latit2|Long2|Width|Depth|Dip|Strike-slip|Dip-slip|Open', 'Position', [10 216 70 180], 'Max', 2);

txt = uicontrol(fig, 'Style', 'text', 'String', 'Fault#', 'Position', [10 390 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Latit1', 'Position', [10 372 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Long1', 'Position', [10 354 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Latit2', 'Position', [10 336 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Long2', 'Position', [10 318 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Width', 'Position', [10 300 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Depth', 'Position', [10 282 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Dip', 'Position', [10 264 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Strike-slip', 'Position', [10 246 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Dip-slip', 'Position', [10 228 70 18]);
txt = uicontrol(fig, 'Style', 'text', 'String', 'Open', 'Position', [10 210 70 18]);

for i = 1:nf
txt        = uicontrol(fig, 'Style', 'text',  'String', num2str(i), 'Position', [ 10+75*i  390 70 18]);
getlat1(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,4),6), 'Position', [ 10+75*i  372 70 18]);
getlon1(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,5),6), 'Position', [ 10+75*i  354 70 18]);
getlat2(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,6),6), 'Position', [ 10+75*i  336 70 18]);
getlon2(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,7),6), 'Position', [ 10+75*i  318 70 18]);
getwid(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,1)), 'Position', [ 10+75*i  300 70 18]);
getdep(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,2)), 'Position', [ 10+75*i  282 70 18]);
getdip(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,3)), 'Position', [ 10+75*i  264 70 18]);
getsss(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,8)), 'Position', [ 10+75*i  246 70 18]);
getdss(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,9)), 'Position', [ 10+75*i  228 70 18]);
getops(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,10)), 'Position', [ 10+75*i  210 70 18]);
end


% Draw Frame and Boxes For Volume Source Elements

if nv>0
	offset = 30+75*(nf+1);
	frm1 = uicontrol(fig, 'Style', 'frame',  'Position',...
	   [offset-5  312  10+75*(nv+1) 99]);

%	multilabels = uicontrol(fig, 'Style', 'edit', ...  
%	'String','VSource#|Latit|Long|Depth|Vol(km^3)', 'Position', ...
%	[offset 316 70 82], 'Max', 2);

	txt = uicontrol(fig, 'Style', 'text', 'String', 'VSource#', 'Position', [offset 390 70 18]);
	txt = uicontrol(fig, 'Style', 'text', 'String', 'Latit', 'Position', [offset 372 70 18]);
	txt = uicontrol(fig, 'Style', 'text', 'String', 'Long', 'Position', [offset 354 70 18]);
	txt = uicontrol(fig, 'Style', 'text', 'String', 'Depth', 'Position', [offset 336 70 18]);
	txt = uicontrol(fig, 'Style', 'text', 'String', 'Vol(km^3)', 'Position', [offset 318 70 18]);

	for j = 1:nv
	  txt2        = uicontrol(fig,'Style','text', 'String', num2str(j), ...
	   'Position', [offset+75*j   390 70 18]);
	  getlatv(j) = uicontrol(fig,'Style','edit','String',...
	    num2str(vgeom(j,1),6), 'Position', [ offset+75*j  372 70 18]);
	  getlonv(j) = uicontrol(fig, 'Style', 'edit',  'String', ...
	    num2str(vgeom(j,2),6), 'Position', [ offset+75*j  354 70 18]);
	  getdepv(j) = uicontrol(fig, 'Style', 'edit',  'String', ...
	    num2str(vgeom(j,3)), 'Position', [ offset+75*j  336 70 18]);
	  getvol(j) = uicontrol(fig, 'Style', 'edit',  'String', ...
	    num2str(vgeom(j,4)), 'Position', [ offset+75*j  318 70 18]);
	end
end

distart =  uicontrol(fig, 'Style', 'push',  'String', 'Forward Model', 'Position', [430 75 125 25], 'CallBack','ForModel', 'BackgroundColor','white', 'ForegroundColor','red');

distart2 =  uicontrol(fig, 'Style', 'push',  'String', 'Estimate Slip', 'Position', [430 50 125 25], 'CallBack','SlipEst', 'BackgroundColor','white', 'ForegroundColor','red');

distart3 =  uicontrol(fig, 'Style', 'push',  'String', 'Estimate Geometry', 'Position', [430 25 125 25], 'CallBack','GeomEst', 'BackgroundColor','white', 'ForegroundColor','red');

distart4 =  uicontrol(fig, 'Style', 'push',  'String', 'Save Geometry', 'Position', [430 150 125 25], 'CallBack','SaveGeom', 'BackgroundColor','white', 'ForegroundColor','black');
