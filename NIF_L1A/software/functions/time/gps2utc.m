function utctime=gps2utc(gpstime)
%GPS2UTC   utctime=gps2utc(gpstime)
%
%Coverts gps time (in seconds) to utc time.

utctime=gpstime+leapsecs;
