function ts=secmomest(x,lat,lon,dep,Lc,Wc,v0,mv0, tauc, M0)
%
%Constructor for the array which keeps track of second moment estimates
%Assumes spherical coordinates
%
%  The way to do this is to make one for each event and then concatenate it to
% the master list and save ie:  newmaster = [masterlist, newevent]; where
% masterlist and newevent are secmomest structures
%
% 
%   x   	The 16? vector of optimal estimates that comes out of minlc.
%   lat		final centroid latitude (deg)
%   lon		final centroid longitude (deg)
%   dep		final centroid depth (km)
%   Lc		estimate of Lc (km)
%   Wc		estimate of Wc (km)
%   v0          estiamte of v0, 3 vector (km/s)
%   mv0		extimate of magnitude v0 (km/s)
%   tauc        estimate of tauc
%   M0		moment in Nm
%
%Set object precendence

    superiorto('double','struct','cell', 'char');

%Do argument checking

    switch nargin

        case 0
            x=[];
            lat=[];
            lon=[];
            dep=[];
        case 10
        otherwise
            error('secmomest requires 10 input arguments.')
    end


    ts=class(struct('x',x, ...
                    'lat',lat, ...
                    'lon',lon, ...
                    'dep',dep, ...
                    'Lc',Lc, ...
                    'Wc',Wc, ...
		    'v0',v0, ...
		    'mv0',mv0, ...
		    'tauc',tauc, ...
		    'M0',M0),'secmomest');