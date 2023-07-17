function cal=decyr2cal(decyr)
%decyr2cal  cal=decyr2cal(decyr)  
%
%Returns calendar date given decimal years.

cal=datevec(datenum(floor(decyr),1,1)+(datenum(floor(decyr)+1,1,0)-datenum(floor(decyr),1,0)).*(decyr-floor(decyr)));
