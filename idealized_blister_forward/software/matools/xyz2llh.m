function [lat, lon, h] = xyz2llh(x, y, z)
%function [lat, lon, h] = xyz2llh(x, y, z)
%
%Converts XYZ coordinates to latitude, longitude, height
%
%Input coordinates in meters. Output angles in radians, height in meters.
%

a = 6378137.0;
f = 1.0/298.2572235630;
e2 = 2*f - f*f;
p = sqrt(x*x + y*y);
r = sqrt(p*p + z*z);
lon = atan2(y,x);

%     ... First iteration on lat and h
%            - assumes h = 0

lat = atan2(z/p, 0.01);
n = a/sqrt((1.0 - e2*sin(lat)^2));
h = p/cos(lat) - n;

%     ... Iterate until h converges (should be quick since h << n)

oldh = -10^9;
iter = 0;
num = z/p;
while (abs(h - oldh) > 0.0001)
   iter = iter + 1;
   oldh = h;
   den = 1.0 - e2*n/(n+h);
   lat = atan2(num,den);
   n = a/sqrt((1.0 - e2*sin(lat)^2));
   h = p/cos(lat) - n;
end
