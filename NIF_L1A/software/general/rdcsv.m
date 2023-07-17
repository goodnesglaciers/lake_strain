% function a=rdcsv(file)
% file format:
%      station name, reference NEU position (North-lat.,East-lon., elevation-meter)
%      date(yyyy-mm-dd), dN, dE, dU, MJD(Julian day)
%                   ...
%
% dN, dE, dU - unit mm. 
% In test website, displacement is in order E, N, U. Should it be N, E, U?
% 
% ZLIU  -- 01/28/09
%
function D=rdcsv(file)
if ~exist(file)
    fprintf('%s\n',[file, ' does not exist!']);
    return;
end

fid=fopen(file, 'r');
A=textscan(fid,'%s%f%f%f',1,'delimiter',',');
B=textscan(fid,'%s %f %f %f %f','delimiter',',');
fclose(fid);

site=char(A{1});
site=upper(site(:,9:12));
lat=A{2};
lon=A{3};
ele=A{4};

a=char(B{1});
ymd=[str2num(a(:,1:4)) str2num(a(:,6:7)) str2num(a(:,9:10))];
ep=cal2decyr(ymd);
de=B{2};
dn=B{3};
du=B{4};
MJD=B{5};

D=struct('site',{site},'lon',lon,'lat',lat,'ele',ele,'e',de,'n',dn,'u',du,'time',ep);
