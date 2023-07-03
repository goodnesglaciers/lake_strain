function cal=sec2cal(sec)
%SEC2CAL   cal=sec2cal(sec) 
%
%Returns caledar date ([year month day hour minute second])
%given seconds before 2000/01/01 12:00.

cal=datevec(sec/86400+730486.5);

% %Days
% 
%     S=mod(sec,86400);
%     days=(sec-S)/86400;
%     cal=datevec(days+730486.5);
% 
% %Hours
% 
%     Sh=mod(S,3600);
%     hours=(S-Sh)/3600;
%     cal=cal+datevec(hours/24);
%     S=S-hours*3600;
% 
% %Minutes
% 
%     Sm=mod(S,60);
%     minutes=(S-Sm)/60;
%     cal=cal+datevec(minutes/1440);
%     S=S-minutes*60;
% 
% %Seconds
% 
%     cal=cal+datevec(S/86400);
%     
% %Tidy up
% 
%     cal(S==0,6)=0;
