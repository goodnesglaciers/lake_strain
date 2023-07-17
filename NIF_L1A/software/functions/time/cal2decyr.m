function decyr=cal2decyr(cal)
%CAL2SEC   decyr=cal2decyr(cal)  
%
%Returns decimal years given a calendar date.
%Input 'cal' should be: [year month day hour minute second]
%(i.e., the dates are stored in the rows of 'cal'.)

cal(:,1)=yy2yyyy(cal(:,1));
cal(:,7)=0;

decyr=cal(:,1)+(datenum(cal(:,1),cal(:,2),cal(:,3),cal(:,4),cal(:,5),cal(:,6))-datenum(cal(:,1),1,1))./(datenum(cal(:,1)+1,1,0)-datenum(cal(:,1),1,0));