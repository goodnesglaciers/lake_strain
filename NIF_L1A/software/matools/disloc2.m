function [u,tilt] = disloc2(nu, disgeom, x)
%
% An M-File for the forward dislocation problem.
%    Input: Material Constants
%		nu = Poisson's ratio
%           Dislocation Geometry in a Global E-N-Up coordinate system
%		disgeom(1) = len = fault length in strike direction (km)
%		disgeom(2) = wid = fault width in dip direction (km)
%		disgeom(3) = dep = depth  of lower edge of fault (km)
%		disgeom(4) = dip = dip angle (degrees)
%		disgeom(5) = strik =  strike, clockwise from N (degrees)
%		disgeom(6) = delE  = East offset of midpoint from origin (km)
%		disgeom(7) = delN  = North offset of midpoint from origin (km)
%		disgeom(8) = ss  =  strike slip motion (m)	
%		disgeom(9) = ds  =  dip slip motion (m)	
%		disgeom(10) = op  =  opening  motion (m)	
%    	    Vector of receiver coordinates at which displacement are computed
%		x  = vector [East (km), North (km)]
%    Output:  
%		u  = Displacements in the Global coordinate system
%		tilt  = Tilts in the Global coordinate system
%
alp = 1-2*nu;
i = sqrt(-1);

%
% handle dip
  dip = disgeom(4);
  if dip == 90
     cs = 0.0;
     sn = 1.0;
  else
     dipr = dip*pi/180;
     cs = cos(dipr);
     sn = sin(dipr);
  end
%
source_geom = [disgeom(3); disgeom(1); 0; disgeom(2); 0; sn; cs; disgeom(8); disgeom(9); disgeom(10)];
%
%strkr = disgeom(5)*pi/180;
strkr = (90 - disgeom(5))*pi/180;
[n,m] = size(x);
%
delE = disgeom(6); delN = disgeom(7);

for k= 1:n,
%    transform the global receiver coordinates into local coordinates
     zlocal = (  (x(k,1)+i*x(k,2)) - (delE+delN*i) )*exp(-i*strkr) + 0.5*disgeom(1);
     xrec = [real(zlocal); imag(zlocal)]; 
 
%  call okadaio
     ut  = okadaio(alp, xrec, source_geom);

%  transformaiton back to global coordinates
     zglob = (ut(1)+i*ut(2) )*exp(i*strkr);
     u(k,1) = real(zglob); u(k,2) = imag(zglob); u(k,3) = ut(3);

     zt = (ut(8)+i*ut(9) )*exp(i*strkr);
     tilt(k,1) = real(zt); tilt(k,2) = imag(zt); 

end
