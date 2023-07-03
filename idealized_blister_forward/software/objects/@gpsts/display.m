function display(ts)
%
%Command windows display for gpsts (GPS Time Series) objects.

var=inputname(1);

ns=size(ts.sites,1);
ne=length(ts.epochs);

if ns > 1
   sitestr='sites';
else
   sitestr='site';
end

if ne > 1
   epochstr='epochs';
else
   epochstr='epoch';
end


fprintf('\n%s = \n\n     GPS Time Series Object\n',var);
fprintf('     [%d %s; %d %s]\n\n',ns,sitestr,ne,epochstr);