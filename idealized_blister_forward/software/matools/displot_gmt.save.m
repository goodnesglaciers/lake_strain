function  displot_gmt(disgeom, latlon1, latlon2) 
%
% Plot surface projection of the dislocation.
%    Input: Dislocation Geometry in a Global E-N-Up coordinate system
%		disgeom(1) = len = fault length in strike direction (km)
%		disgeom(2) = wid = fault width in dip direction (km)
%		disgeom(3) = dep = fault depth (km)
%		disgeom(4) = dip = dip angle (degrees)
%		disgeom(5) = strik =  strike, counter clockwise from E (degrees)
%		disgeom(6) = delE  = East offset of midpoint from origin (km)
%		disgeom(7) = delN  = North offset of midpoint from origin (km)
%		disgeom(8) = ss  =  strike slip motion (m)	
%		disgeom(9) = ds  =  dip slip motion (m)	
%		disgeom(10) = op  =  opening  motion (m)	
%%
len =disgeom(1); wid = disgeom(2); dep = disgeom(3); delE=disgeom(6); delN=disgeom(7);
% corners in local coordinates
i = sqrt(-1);
dipr = disgeom(4)*pi/180;
ll = 0 +i*0;  lr = len + i*0; ul = 0 +i*wid*cos(dipr); ur = len + i*wid*cos(dipr);
%
%transform corners to global coordinates
strkr = disgeom(5)*pi/180;
% radius Rp (Polar) and Rm (Meridianal) in meters
Rp = 6378206.4;
Rm = Rp*cos(latlon1(1)*pi/180);
%%
gll = (latlon1(1) + i*latlon1(2));  
glr = (latlon2(1) + i*latlon2(2)); 
%offset =  ul*1e3*(cos(strkr)/Rm + i*sin(strkr)/Rp)*180/pi;
offset =  ul*1e3*(cos(strkr)/Rm + i*sin(strkr)/Rm)*180/pi;

%factor of 1e3 converts km to meters, then convert radians to degrees
gul = (latlon1(1) + i*latlon1(2)) + offset;
gur = (latlon2(1) + i*latlon2(2)) + offset;
%
%plot second line just in from the upper edge
%offset2 =  0.9*ul*1e3*(cos(strkr)/Rm + i*sin(strkr)/Rp)*180/pi;
offset2 =  0.9*ul*1e3*(cos(strkr)/Rm + i*sin(strkr)/Rm)*180/pi;
gulm = (latlon1(1) + i*latlon1(2)) + offset2;
gurm = (latlon2(1) + i*latlon2(2)) + offset2;

outline = [latlon1(2), latlon1(1); latlon2(2), latlon2(1); imag(gur), real(gur); imag(gul), real(gul);latlon1(2), latlon1(1); imag(gulm), real(gulm); imag(gurm), real(gurm) ];

%save temp.gmtlin outline -ascii
uo = fopen('temp.gmtlin', 'a');
fprintf(uo, '>  \n ');
fprintf(uo, '%12.8f %12.8f\n', outline');
fclose(uo);

