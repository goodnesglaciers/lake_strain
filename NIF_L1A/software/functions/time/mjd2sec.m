function sec=mjd2sec(mjd)
%MJD2SEC    sec=mjd2sec(mjd) 
%
%Converts Modified Julian Days to J2K seconds.

sec=(mjd-51544.5)*86400;