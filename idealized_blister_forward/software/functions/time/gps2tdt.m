function tdttime=gps2tdt(gpstime)
%GPS2TDT     tdttime=gps2tdt(gpstime)
%
%Converts gps time (seconds) to Terrestrial Dynamic
%time.

tdttime=gps2tai(gpstime)+32.184;