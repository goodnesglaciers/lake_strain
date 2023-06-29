function [lat,lon]=clarke(dlat0,dlon0,x,y)
% clarke
%   convert x,y in meters to lat lon in degrees
%   Based on Clarke spheroid, 1866 as in Bowditch
%   emulates ALVIN TRANSLAT utility
%   
% Usage [lat,lon]=clarke(dlat0,dlon0,x,y)
%   input:    lat. lon. origin (decimal degrees)
%             x, y coords in meters from net origin
%   output:   lat,lon of converted position
% Maurice A. Tivey
%   default origing (EPR): 9.133333, -104.33333

% dlat0 = 9.133333333;
% dlon0 = -104.33333333;

%--------work out scaling (metres per minute of lat, lon):
      radlat = dlat0/57.2957795;
      c1 = cos(radlat);
      c2 = cos(2.*radlat);
      c3 = cos(3.*radlat);
      c4 = cos(4.*radlat);
      c5 = cos(5.*radlat);
      c6 = cos(6.*radlat);
      sclat = (111132.09-566.05*c2+1.20*c4-.002*c6);
      sclon = (111415.13*c1-94.55*c3+.012*c5);
%
 lat=dlat0+(y./sclat);
 lon=dlon0+(x./sclon);
