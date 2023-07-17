function [data1,data2]=bpfilter(data,lcorner,hcorner,dt,order)

%BPFILTER  bandpass filters seismograms (2-pass)
%
%  data=bpfilter(data);  Input data matrix has seismograms arranged as column
%  vectors.  Uses default values of lcorner = 1 Hz, hcorner = 9 Hz, and dt = 
%  0.01 seconds.
%
%  data=bpfilter(data,lcorner,hcorner,dt,order); is the complete function call to
%  specify desired paramters.  To use default values for any argument input
%  the empty matrix.

%  David Schaff 8-8-98	modified from Doug Dodge.

if nargin == 1
   lcorner=1;
   hcorner=9;
   dt=0.01;
   order=2;
elseif nargin < 5
   error('Need either 1 or 5 input arguments.')
else
   if isempty(lcorner), lcorner=1; end
   if isempty(hcorner), hcorner=9; end
   if isempty(dt), dt=0.01; end
   if isempty(order), order=2; end
end


   [b,a]=butter(order,[lcorner*dt*2 hcorner*dt*2]);	
   data2=filtfilt(b,a,data);
   data1=filter(b,a,data);
