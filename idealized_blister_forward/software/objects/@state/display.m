function display(x)
%
%Display function for state objects

var=inputname(1);

s=size(x);

if s(1)==-1 & s(2)==-1
    fprintf('\n\tEmpty state object.\n\n') 
    return
end

if isempty(x.swapdir)
    x.swapdir='<N/A>';
    db=0;
else
    w=dir([x.swapdir,'\*.mat']);
    db=sum(cat(1,w.bytes));
end
fprintf('\n\tProperties of state object %s:\n',var) 
fprintf('\n\t\tMaximum i index: %d\n',s(1)) 
fprintf('\t\tMaximum j index: %d\n',s(2))
fprintf('\t\tTotal number of entries: %d\n',nnz(x.index));
fprintf('\t\tSwap Directory: %s\n',x.swapdir) 
w=whos('x');
fprintf('\t\tBytes in RAM: %d\n',w.bytes) 
fprintf('\t\tBytes on disk: %d\n\n',db) 