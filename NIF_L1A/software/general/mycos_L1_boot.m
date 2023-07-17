	function [xopt,stdE]=mycos_L1(x,xdata,dt_stdE);
% L1 norm minimizations of xdata=B+A*cos(x-theta);
% just a dumb grid search
%
%  x=[B,A,theta];
%  xdata = azimuths in degrees
%  y=predicted dts
%  B=cmt time offset
%  A=dtmax
%  theta=relative location azimuth
%  dt_stdE is the standard deviation of the travel time observations which are
%  assumed independent and identical variance.
%  stdE is the standard error of each element of xopt


% best estimate
  xopt=mycos_L1_fine(x,xdata);              %L1 norm often better

% basically want to set up a population with random samples and then
% call the original xopt_L1

BS=100;
N=length(xdata)
zz=ceil(N*rand(BS,N)); % scale?

Bsv=zeros(BS,1); Asv=zeros(BS,1); theta_sv=zeros(BS,1);

for i=1:BS
  %make xi xdatai
  % normally in a bootstrap this is a resample of the actual data, but we don't
  % have very many samples (5-10) and their azimuthal distribution is key to constraining
  % anything, so I don't think any info is gained from reducing the azimuthal distrubiton
  % so we follow shearer 1997 jgr / Billings 1994 bssa.  We assume gaussian picking errors
  % like billings
  xdatai=xdata+randn(1,N)*dt_stdE;
  xi=x;  %no error in azimuth?
  %
  [xopt_i]=mycos_L1_fine(xi,xdatai);
  Bsv(i)=xopt_i(1); Asv(i)=xopt_i(2); theta_sv(i)=xopt_i(3);
end

% now calculate sample standard devaition of each and return

clear stdE
stdE(1)=std(Bsv); stdE(2)=std(Asv); stdE(3)=std(theta_sv);

end
