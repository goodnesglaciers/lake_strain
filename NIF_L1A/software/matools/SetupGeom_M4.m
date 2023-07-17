nf = eval( get(nfaults,'string') );
nv = eval( get(nvols,'string') );

%test whether to resize GEOM:
        [nf_geom,nel] = size(geom);
	[nv_geom,nvel] = size(vgeom);
	if nf < nf_geom 
		frm1 = uicontrol(fig,'Style','frame','BackgroundColor',[0,0,0.8], ...
			'Position',[5  215  15+75*(nf_geom+max(nv,nv_geom)+2) 185]);
%			'Position',[5  215  15+75*(nf_geom+1) 185]);
			geom = geom(1:nf,:);
		end
		if nv > 0
		frm2 = uicontrol(fig, 'Style', 'frame', 'BackgroundColor',...
		[0,0,0.8],'Position',[25+75*(nf_geom+1)  312  10+75*(nv_geom+1) 88]);
		end
	if nf > nf_geom 
		geom = [geom;zeros(nf-nf_geom,10)];
	end

%test whether to resize VGEOM:
if nv>0  | nv_geom >0      
	if nv < nv_geom 
	frm2 = uicontrol(fig, 'Style', 'frame', 'BackgroundColor',...
		[0,0,0.8],'Position',[25+75*(nf_geom+1)  312  10+75*(nv_geom+1) 88]);
		vgeom = vgeom(1:nv,:);
	end
	if nv > nv_geom 
		vgeom = [vgeom;zeros(nv-nv_geom,4)];
	end
end

% Draw Frame and Boxes For Dislocation Elements
frm1 = uicontrol(fig, 'Style', 'frame',  'Position',[5  215  15+75*(nf+1) 185]);

multilabels = uicontrol(fig, 'Style', 'edit',  'String','Fault#|Latit1|Long1|Latit2|Long2|Width|Depth|Dip|Strike-slip|Dip-slip|Open', 'Position', [10 216 70 180], 'Max', 2);

for i = 1:nf
txt        = uicontrol(fig, 'Style', 'text',  'String', num2str(i), 'Position', [ 10+75*i  378 70 16]);
getlat1(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,4),6), 'Position', [ 10+75*i  362 70 16]);
getlon1(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,5),6), 'Position', [ 10+75*i  346 70 16]);
getlat2(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,6),6), 'Position', [ 10+75*i  330 70 16]);
getlon2(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,7),6), 'Position', [ 10+75*i  314 70 16]);
getwid(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,1)), 'Position', [ 10+75*i  298 70 16]);
getdep(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,2)), 'Position', [ 10+75*i  282 70 16]);
getdip(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,3)), 'Position', [ 10+75*i  266 70 16]);
getsss(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,8)), 'Position', [ 10+75*i  250 70 16]);
getdss(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,9)), 'Position', [ 10+75*i  234 70 16]);
getops(i) = uicontrol(fig, 'Style', 'edit',  'String', num2str(geom(i,10)), 'Position', [ 10+75*i  218 70 16]);
end


% Draw Frame and Boxes For Volume Source Elements

if nv>0
	offset = 30+75*(nf+1);
	frm1 = uicontrol(fig, 'Style', 'frame',  'Position',...
	   [offset-5  312  10+75*(nv+1) 88]);
	multilabels = uicontrol(fig, 'Style', 'edit', ...  
	'String','VSource#|Latit|Long|Depth|Vol(km^3)', 'Position', ...
	[offset 316 70 82], 'Max', 2);

	for j = 1:nv
	  txt2        = uicontrol(fig,'Style','text', 'String', num2str(j), ...
	   'Position', [offset+75*j   378 70 16]);
	  getlatv(j) = uicontrol(fig,'Style','edit','String',...
	    num2str(vgeom(j,1),6), 'Position', [ offset+75*j  362 70 16]);
	  getlonv(j) = uicontrol(fig, 'Style', 'edit',  'String', ...
	    num2str(vgeom(j,2),6), 'Position', [ offset+75*j  346 70 16]);
	  getdepv(j) = uicontrol(fig, 'Style', 'edit',  'String', ...
	    num2str(vgeom(j,3)), 'Position', [ offset+75*j  330 70 16]);
	  getvol(j) = uicontrol(fig, 'Style', 'edit',  'String', ...
	    num2str(vgeom(j,4)), 'Position', [ offset+75*j  314 70 16]);
	end
end

distart =  uicontrol(fig, 'Style', 'push',  'String', 'Forward Model', 'Position', [430 75 125 25], 'CallBack','ForModel', 'BackgroundColor','white', 'ForegroundColor','red');

distart2 =  uicontrol(fig, 'Style', 'push',  'String', 'Estimate Slip', 'Position', [430 50 125 25], 'CallBack','SlipEst', 'BackgroundColor','white', 'ForegroundColor','red');

distart3 =  uicontrol(fig, 'Style', 'push',  'String', 'Estimate Geometry', 'Position', [430 25 125 25], 'CallBack','GeomEst', 'BackgroundColor','white', 'ForegroundColor','red');

distart4 =  uicontrol(fig, 'Style', 'push',  'String', 'Save Geometry', 'Position', [430 150 125 25], 'CallBack','SaveGeom', 'BackgroundColor','white', 'ForegroundColor','black');
