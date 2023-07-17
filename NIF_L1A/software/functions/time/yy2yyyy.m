function YR=yy2yyyy(yr,pivot)
%YY2YYYY     YR=yy2yyyy(yr,pivot)
%
%Converts two digit years to four digit years
%given a 'pivot' year.  Pivot defaults to 80
%if omitted.

if nargin<2
    pivot=80;
end

I=yr<pivot;
J=yr>=pivot & yr < 100;

YR=yr;
YR(I)=YR(I)+2000;
YR(J)=YR(J)+1900;