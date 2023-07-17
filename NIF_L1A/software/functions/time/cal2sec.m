function sec=cal2sec(cal)
%CAL2SEC   sec=cal2sec(cal)  
%
%Returns seconds before 2000/01/01 12:00 given a calendar date.
%Input 'cal' should be: [year month day hour minute second]
%(i.e., the dates are stored in the rows of 'cal'.)


cal(:,1)=yy2yyyy(cal(:,1));
cal(:,7)=0;

sec=(datenum(cal(:,1),cal(:,2),cal(:,3),0,0,0)-730486.5)*86400+cal(:,4)*3600+cal(:,5)*60+cal(:,6);
