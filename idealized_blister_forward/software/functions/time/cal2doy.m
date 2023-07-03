function doy=cal2doy(cal)
%CAL2DOY   doy=cal2doy(cal)  
%
%Returns day of year given a calendar date.
%Input 'cal' should be: [year month day].
%(i.e., the dates are stored in the rows of 'cal'.)


doy=datenum(cal(:,1),cal(:,2),cal(:,3))-datenum(cal(:,1),1,0);
    