function bit=isleapyear(yr)
%ISLEAPYEAR     isleapyear(yr)
%
%Returns 1 if 'yr' is a leap year, 0 otherwise

yr=floor(yr);
bit=(mod(yr,4)==0 & mod(yr,100)~=0 | mod(yr,400)==0);