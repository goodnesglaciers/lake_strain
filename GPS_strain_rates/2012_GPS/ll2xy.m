function [x,y]=iclarke(dlat0,dlon0,lat,lon)
% iclarke
%   convert lat lon to x y coords in meters given a
%   net origin in degrees
%   based on clarke spheroid, 1866 as in Bowditch
%   emulates ALVIN TRANSLAT utility
%   
% Usage [x,y]=iclarke(dlat0,dlon0,lat,lon)
%   input:    dlat0,dlon0 net origin (decimal degrees)
%             lat,lon positions to translate in degrees
%   output:   x,y, coordinates within net in meters from origin
% Maurice A. Tivey
% see clarke.m

%--------work out scaling (metres per minute of lat, lon):
      radlat = dlat0/57.2957795;
      c1 = cos(radlat);
      c2 = cos(2.*radlat);
      c3 = cos(3.*radlat);
      c4 = cos(4.*radlat);
      c5 = cos(5.*radlat);
      c6 = cos(6.*radlat);
      sclat = (111132.09-566.05*c2+1.20*c4-.002*c6)
      sclon = (111415.13*c1-94.55*c3+.012*c5)
%
 y=(lat-dlat0).*sclat;
 x=(lon-dlon0).*sclon;
