function varargout=size(ts)
%
%Size function for gpsts (GPS Time Series) objects.

ns=size(ts.sites,1);
ne=length(ts.epochs);

switch nargout

    case {0,1}
        varargout{1}=[ns ne];

    case 2
        varargout{1}=ns;
        varargout{2}=ne;

end

