function [T] = rotmatrix(lat, lon)

% function [T] = rotmatrix(lat, lon)
%
% Calculates the transform matrix from XYZ system 
% to ENU system given the reference latitude and 
% longtitude
%
% 02/24/99 Yosuke Aoki

sphi = sin(lat*pi/180);
cphi = cos(lat*pi/180);
slmb = sin(lon*pi/180);
clmb = cos(lon*pi/180);

T = [-slmb, clmb, 0;  -sphi*clmb , -sphi*slmb, cphi; ... 
cphi*clmb, cphi*slmb, sphi];

