	function y=mycos(x,xdata);
% cosfit y=B+A*cos(x-theta);
%
%  x=[B,A,theta];
%  xdata = azimuths
%  y=predicted dts
%  B=cmt time offset
%  A=dtmax
%  theta=relative location azimuth

   y=x(1)+x(2)*cos((x(3)-xdata)*(pi/180));
end
