
function [dist,az,delta] = distaz(lta,lna,ltb,lnb);

% DISTAZ computes the distance and azimuth between two locations on the 
% surface of the earth. The azimuth is measured clockwise from the 
% local direction of North at the location (lta,lna) to the great circle arc 
% connecting the two locations.  The location coordinates are assumed to be 
% in decimal degrees.  East longitudes are positive.  To determine
% the back-azimuth from a seismic station to an earthquake epicenter,
% let (lta,lna) be the coordinates of the seismic station and (ltb,lnb) be
% the cooordinates of the epicenter. 
%
% Reference:  Gubbins, "Seismology and Plate Tectonics", Appendix A.
%
% USAGE: [dist,az,delta] = distaz(lta,lna,ltb,lnb);
%                                                            j.a.collins
%-----------------------------------------------------------------------

radius_earth = 6371.0;

%convert to radians
clta = deg2rad(90-lta);   % co-latitude in radians
lna = deg2rad(lna);
cltb = deg2rad(90-ltb); 
lnb = deg2rad(lnb);

% determine distance
delta = acos( cos(clta)*cos(cltb) + sin(clta)*sin(cltb)*cos(lna-lnb) );
dist = delta * radius_earth;
%%dist = delta;
if (abs(delta) < eps)
    az = 0.0;
    disp('Input coordinates are identical!');
    return;
end

% determine azimuth
cos_az = 1/sin(delta) * (sin(clta)*cos(cltb) - ...
                         sin(cltb)*cos(clta)*cos(lnb-lna));
sin_az = 1/sin(delta) * (sin(cltb)*sin(lnb-lna));
az = atan2(sin_az,cos_az);
if (az < 0)
    az = az + 2*pi;
end
az = rad2deg(az);

delta = rad2deg(delta);
return


