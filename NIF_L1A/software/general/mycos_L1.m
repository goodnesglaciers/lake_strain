	function [xopt]=mycos_L1(x,xdata);
% L1 norm minimizations of xdata=B+A*cos(x-theta);
% just a dumb grid search
%
%  x=[B,A,theta];
%  xdata = azimuths in degrees
%  y=predicted dts
%  B=cmt time offset
%  A=dtmax
%  theta=relative location azimuth

ndata=length(xdata);
misfit=zeros(101,201,72);
mis=999999;
i1=0; i2=0; i3=0;
for x1i=-50:1:50   %loop over mean
 i1=i1+1;
 i2=0; i3=0;
 for x2i=-60:.6:60  % loop over peak-peak  %20*3.7km/s =75 km, has to allow more than cutoff
    i2=i2+1;
    i3=0;
    for x3i=0:5:175 % loop over azimuths
     i3=i3+1;
     for j=1:ndata;
      yhat=x1i+x2i*cos((x3i-x(j))*(pi/180));
      misfit(i1,i2,i3)=misfit(i1,i2,i3)+abs(xdata(j)-yhat);
     end
     if(misfit(i1,i2,i3) <=mis)
      mis=misfit(i1,i2,i3);
      xopt=[x1i, x2i, x3i];
     end

    end  %azimuths
 end %peak-peak
end %mean

%keyboard