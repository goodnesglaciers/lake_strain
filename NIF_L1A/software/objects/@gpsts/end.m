function output=end(ts,i,j)
%
%end function for gps ts epochs

switch i
    case 1
        output=size(ts.sites,1);
    case 2
        output=length(ts.epochs);
end