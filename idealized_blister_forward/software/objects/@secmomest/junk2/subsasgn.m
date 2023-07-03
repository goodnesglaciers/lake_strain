function ts=subsasgn(ts,arg,value)
%
%Subscripted assignment method for gpsts objects.  Mainly for data
%deletion.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Apr 24, 2001   Peter Cervelli		 Original Code
%
%-------------------------------------------------------------

%Check for proper assignment

	if ~isempty(value)
        error('Only deletion (empty set assignment) is defined for gpsts objects.')
	end

%Parse the index 

	[siteindex,epochindex]=deal(arg(1).subs{1:2});
	masterindex=~ParseIndex(siteindex,epochindex,ts.sites,ts.epochs,ts.siteindex,ts.epochindex);

%Perform the deletion

     ts.d=ts.d(masterindex);
     [i,j,s]=find(ts.dcov(masterindex,masterindex));
     ts.dcov=sparse(i,j,s,length(masterindex),length(masterindex));
     oldsites=ts.sites;
     ts.sites=ts.sites(findindex(ts.siteindex(masterindex)));
     ts.epochs=ts.epochs(findindex(ts.epochindex(masterindex)));
     ts.epochindex=reindex(ts.epochindex(masterindex));
     ts.siteindex=reindex(ts.siteindex(masterindex));
     [NULL,b]=intersect(oldsites,setdiff(oldsites,ts.sites));
     ts.apcoords(:,b)=[];

%BEGIN sub-function findindex ------------------------------------------

function I=findindex(J)

    I(J)=1;
    I=logical(I);

%END sub-function findindex --------------------------------------------


%BEGIN sub-function reindex --------------------------------------------

function I=reindex(J);
    
    U=unique(J);
    K(U,1)=(1:length(U))';
    I=K(J);
    
%END sub-function reindex --------------------------------------------
