function output=subsref(x,arg)
%
%Subscript interpreter for state objects.

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

if i>m | j>n
    index=0;
else
    index=x.index(i,j);
end

if index==0
    error('No such index.')
end

output=x.value{index};

if ischar(output)
    output=swapmanager(i,j,output,x.swapdir);
end

if size(arg,2)==2
    output=subsref(output,arg(2));
end


