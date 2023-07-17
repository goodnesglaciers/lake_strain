function display(ts)
%
%Command windows display for gpsts (GPS Time Series) objects.

var=inputname(1);

ne=length(ts);  %number of earthquakes


fprintf('\n%s = \n\nSecond Moments Object\n',var);
fprintf('%s %s \n', num2str(ne), 'Earthquakes');