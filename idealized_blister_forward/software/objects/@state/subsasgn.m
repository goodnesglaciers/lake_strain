function x=subsasgn(x,arg,value)
%
%Subscript assignment interpreter for state objects.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Apr 14, 2001   Peter Cervelli		 Original Code
%
%-------------------------------------------------------------


i=arg(1).subs{1}+1;
j=arg(1).subs{2}+1;

[m,n]=size(x.index);

if i>m | j>n | x.index(i,j)==0
    x.index(i,j)=length(x.value)+1;
end

if ~isempty(x.swapdir)
    value=swapmanager(i,j,value,x.swapdir);
end

x.value{x.index(i,j)}=value;