function taitime=gps2tai(gpstime)
%GPS2TAI    taitime=gps2tai(gpstime)
%
%Convert gps time (in seconds) to International Atomic
%Time.

taitime=gpstime+19;