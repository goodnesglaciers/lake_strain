function output=subsref(ts,arg)
%
%Subscript interpreter for gpsts objects.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Feb 25, 2001   Peter Cervelli		 Original Code
%
%-------------------------------------------------------------


%Parse the index -- if no index is passed, assume all (:,:)

    if strcmp(arg(1).type,'.')
        if size(arg,2)==1 | (size(arg,2)==2 & strcmp(arg(2).type,'()'))
            if strmatch(arg(1).subs,{'d','dcov','sites','epochs','siteindex','epochindex','apcoords'},'exact')
                output=builtin('subsref',ts,arg);
                return
            end
        end
        arg=[struct('type','()','subs',{{':',':'}}) arg(1:end)];
    end

    [siteindex,epochindex]=deal(arg(1).subs{1:2});

    masterindex=find(ParseIndex(siteindex,epochindex,ts.sites,ts.epochs,ts.siteindex,ts.epochindex));

    if isempty(masterindex)
        error('The station and epoch combination you requested is empty.')
    end

%Switch over the request specified by the first word after the '.'
    
    if size(arg,2)>1
        request=arg(2).subs;
    else
        request='gpsts';
    end

    switch lower(request)

        case {'apcoords','apcoord'}
            
            output=ts.apcoords(:,findindex(ts.siteindex(masterindex)));

        case 'd'

            output=ts.d(masterindex);

        case 'dcov'

            output=ts.dcov(masterindex,masterindex);
            if nnz(output)/prod(size(output)) >= 0.6
                output=full(output);
            end
            
        case {'date','cal'}

            output=decyr2cal(ts.epochs(findindex(ts.epochindex(masterindex))));

        case {'epochs','epoch','t'}

            output=ts.epochs(findindex(ts.epochindex(masterindex)));

        case {'sites','sitelist','site','stations','station'}
            
            output=ts.sites(findindex(ts.siteindex(masterindex)),:);

        case 'epochindex'
            
            output=reindex(ts.epochindex(masterindex));

        case 'epochmindex'
            
            output=ts.epochindex(masterindex);

        case 'siteindex'

            output=reindex(ts.siteindex(masterindex));

        case 'sitemindex'
            
            output=ts.siteindex(masterindex);

        case 'ts'

            output=reshape(ts.d(masterindex),3,length(masterindex)/3);

        case {'gpsts','local'}

            [i,j,s]=find(ts.dcov(masterindex,masterindex));

            output=class(struct('d',ts.d(masterindex), ...
                                'dcov',sparse(i,j,s), ...
                                'sites',{ts.sites(findindex(ts.siteindex(masterindex)),:)}, ...
                                'epochs',ts.epochs(findindex(ts.epochindex(masterindex))), ...
                                'siteindex',reindex(ts.siteindex(masterindex)), ...
                                'epochindex',reindex(ts.epochindex(masterindex)), ...
                                'apcoords',ts.apcoords(:,findindex(ts.siteindex(masterindex)))),'gpsts');

        otherwise

            output=[];
    end

    if size(arg,2)==3
        output=subsref(output,arg(3));
    end

%BEGIN sub-function findindex

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

