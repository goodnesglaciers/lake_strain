function [mjd]=cal2mjd(cal)
%CAL2MJD   mjd=cal2mjd(cal)  
%
%Returns Modified Julian Data given a calendar date.
%Input 'cal' should be: [year month day hour minute second]
%(i.e., the dates are stored in the rows of 'cal'.)


cal(:,1)=yy2yyyy(cal(:,1));
cal(:,7)=0;
mjd=datenum(cal(:,1),cal(:,2),cal(:,3),cal(:,4),cal(:,5),cal(:,6))-678942;
