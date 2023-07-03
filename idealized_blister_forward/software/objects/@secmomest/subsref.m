function b=subsref(a,s)
%
%
%
%


switch s.type
case 'lat'
 b=get(s,'lat')
otherwise
  error('WHAT THE FUCK')
end
