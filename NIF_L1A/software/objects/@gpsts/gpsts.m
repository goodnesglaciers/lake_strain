function ts=gpsts(d,dcov,sites,epochs,siteindex,epochindex,apcoords)
%GPSTS     ts=gpsts(d,dcov,sites,epochs,siteindex,epochindex,apcoords)
%
%Constructor for the gps time series (gpsts) class.  If called
%with no input arguments, an empty object is returned.
%
%The general syntax for a gpsts object is:
%
%  ts(siteindex,epochindex).dtype
%
%Data are returned based on the conjunction of the site index and the epoch
%index.
%
%The siteindex can be: (1) a numerical index into the master list of sites
%for the gpsts object in question (i.e., an index into ts.sites), (2) a
%site name, or (3) a braced ({}), comma separated list of site names.
%
%The epochindex can be: (1) a numerical index into the master list of epochs
%for the gpsts object in question (i.e., an index into ts.epochs) or (2) a
%braced ({}), comma separated pair of times (decimal years) representing the
%desired interval.
%
%The colon (:) and the 'end' keyword function as normal for both indices. If
%no indices are specified, (:,:) is assumed.
%
%Valid entries for 'dtype':
%
%   d              data vector
%   dcov           data covariance
%   sites          list of sites
%   apcoords       3xn matrix of approximate site coordinates
%   epochs         a vector of observation times
%   siteindex      a vector of length(ts(siteindex,epochindex).d) of indices into
%                  ts(siteindex,epochindex).sites.
%   epochindex     a vector of length(ts(siteindex,epochindex).d) of indices into
%                  ts(siteindex,epochindex).epochs.
%   sitemindex     a vector of length(ts(siteindex,epochindex).d) of indices into
%                  ts.sites.
%   epochmindex    a vector of length(ts(siteindex,epochindex).d) of indices into
%                  ts.epochs.
%   ts             a gpsts object containing only the data indexed by siteindex
%                  and epochindex
%
%If no 'dtype' is specified, 'ts' is assumed.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Apr 11, 2001   Peter Cervelli        Reading sinex files now
%                                        occurs in a separate function.
%                                        When called with no arguments
%                                        gpsts returns and empty object.
%   Feb 28, 2001   Peter Cervelli		 Original Code
%
%-------------------------------------------------------------


%Set object precendence

%    superiorto('double','sparse','struct','cell','char','inline');

%Do argument checking

    switch nargin

        case 0
            d=[];
            dcov=[];
            sites=[];
            epochs=[];
            siteindex=[];
            epochindex=[];
            apcoords=[];

        case 7
        otherwise
            error('gpsts requires 7 input arguments.')
    end


    ts=class(struct('d',d, ...
                    'dcov',dcov, ...
                    'sites',{sites}, ...
                    'epochs',epochs, ...
                    'siteindex',siteindex, ...
                    'epochindex',epochindex, ...
                    'apcoords',apcoords),'gpsts');