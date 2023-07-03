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

whos

    if strcmp(arg(1).type,'.')
        if size(arg,2)==1 | (size(arg,2)==2 & strcmp(arg(2).type,'()'))
            if strmatch(arg(1).subs,{'x','lat','lon','dep','Lc','Wc','v0','mv0','tauc','M0'},'exact')
                output=builtin('subsref',ts,arg);
                return
            end
        end
        arg=[struct('type','()','subs',{{':',':'}}) arg(1:end)];
    end
output
%Switch over the request specified by the first word after the '.'
    
    if size(arg,2)>1
        request=arg(2).subs;
    else
        request='gpsts';
    end
%
    switch lower(request)

        case 'lat'

            output=ts.lat(masterindex);

        case 'lon'

            output=ts.lon(masterindex);           

        otherwise

            output=[];
    end


