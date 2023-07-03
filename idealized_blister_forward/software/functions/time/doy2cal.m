function cal=doy2cal(doy,yr)
%DOY2CAL   cal=doy2cal(doy,yr)  
%
%Returns year, month, and day given the day of year and an
%optional year.  If the year is omitted, current year is assumed.

%Check input arguments

if nargin < 1 | nargin > 2
    help doy2ymd
    return
end

if nargin ==1
    yr=clock;
    yr=yr(1);
end

%Calculate calendar date

cal=datevec(datenum(yr(:),01,00)+doy(:));
cal=cal(:,1:3);