function cal=sinex2cal(st)
%SINEX2CAL    cal=sinex2cal(st)
%
%Converts "sinex" time (YY:DOY:SECOND) to calendar time

st(:,1)=yy2yyyy(st(:,1));
cal=datevec(datenum(doy2cal(st(:,2),st(:,1)))+st(:,3)/86400);