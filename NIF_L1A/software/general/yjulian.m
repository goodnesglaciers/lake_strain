function date=yjulian(dates);

%YJULIAN  returns decimal year given the year, month, and date
%
%  year=yjulian(date);  date is a matrix with columns as, date=[yr mo day].
%  Computes the julian day corresponding to the date and converts to a 
%  decimal year by dividing the julian day by 365 days.  Accounts for leap 
%  years which have 366 days.
%
%  See also TIMEBAR

%  David Schaff 9-4-98



y=dates(:,1);
m=dates(:,2);
d=dates(:,3);

mo = 31*ones(12,1);
II=[4 6 9 11];
mo(II)=30*ones(size(II));
mo(2)=28;
days = 365;
LeapYear = ~rem(y,4); % incorrect in general? no OK

numDates=length(y);
date=zeros(numDates,1);

for i=1:numDates
   leapyr=LeapYear(i);
   month=m(i);
   day=d(i);
   if leapyr
      mo(2) = 29;
      days = 366;
   else
      mo(2) = 28;
      days = 365;
   end
   if month == 1
      mday = 0;
   else
      mday = sum(mo(1:month-1));
   end
   tday = mday + day;
   yjday = tday/days;
   %date(i) = y(i) + yjday;
   date(i)=tday;
end
