function x=state(value,swapdir)
%STATEVEC    x=state(value,swapdir)
%
%Creates a state object.  If no inputs are passed then
%an empty object is returned.  Optional swapdir input is a path
%to a file where the data in the state object will be
%stored (i.e., not in RAM).

%Set object precendence

    superiorto('double','struct','cell','char','inline');

%Switch over input arguments

	if nargin < 2
	    swapdir=[];
	end

%Create state object
    
    x=class(struct('index',sparse(1,1,1), ...
                   'value',{{[]}}, ...
                   'swapdir',swapdir),'state');

    if nargin>0
        x=subsasgn(x,struct('type','()','subs',{{0,0}}),value);
    end
