	function y=mycos(x,xdata,partials);
% cosfit y=B+A*cos(x-theta);
%
%  x=[B,A,theta];
%  xdata = azimuths
%  y=predicted dts
%  B=cmt time offset
%  A=dtmax
%  theta=relative location azimuth
   whos
   for j=1:length(xdata)
    y(j)=x(1)+partials(j)*x(2)*cos((x(3)-xdata(j))*(pi/180));
   end
end
