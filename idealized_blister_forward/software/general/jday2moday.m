function [month,day]=jday2moday(year,jday);

%  [month,day]=jday2moday(year,jday);

%YJULIAN  returns decimal year given the year, month, and date
%
%  year=yjulian(date);  date is a matrix with columns as, date=[yr mo day].
%  Computes the julian day corresponding to the date and converts to a 
%  decimal year by dividing the julian day by 365 days.  Accounts for leap 
%  years which have 366 days.
%
%  See also TIMEBAR

%  David Schaff 9-4-98




mo = 31*ones(12,1);
II=[4 6 9 11];
mo(II)=30*ones(size(II));
mo(2)=28;
days = 365;
LeapYear = ~rem(year,4); % incorrect in general? no OK



numDates=length(year);
month=zeros(numDates,1); day=month;

for i=1:numDates
   leapyr=LeapYear(i);
   if leapyr
      mo(2) = 29;
   else
      mo(2) = 28;
   end
totdays=cumsum(mo);
idx=max(find(totdays < jday(i))) + 1;
if isempty(idx)
   month(i)=1;
   day(i)=jday(i);
else
   month(i)=idx;
   day(i)=jday(i)-totdays(idx-1);
end

end
