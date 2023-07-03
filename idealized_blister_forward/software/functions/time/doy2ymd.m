function ymd=doy2ymd(doy,yr)
%DOY2YMD(doy,yr)     cal=doy2ymd(doy,yr)
%
%Returns yymmmdd given doy.  'yr' should be the
%same size as 'doy'.  If 'yr' is omitted, current
%year is used.

%Check input arguments

if nargin < 1 | nargin > 2
    help doy2ymd
    return
end

if nargin ==1
    yr=clock;
    yr=yr(1);
end

%Get calendar date

d=datenum(yr(:),01,00)+doy(:);
ymd=[datestr(d,'yy'),lower(datestr(d,'mmm')),datestr(d,'dd')];

