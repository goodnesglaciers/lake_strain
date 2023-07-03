function doy=decyr2doy(decyr)
%decyr2doy    doy=decyr2doy(decyr)
%
%Converts decimal years to days of year.

daysinyear=isleapyear(decyr)+365;
doy=str2num(sprintf('%3.6f\n',(decyr-floor(decyr)).*daysinyear))+1;
