function  volplot_gmt(vgeom)

% function to plot volume sources on GMT plots

[nv, nvels] = size(vgeom);

uo = fopen('volums.gmtlin', 'a');
for i=1:nv
	fprintf(uo, '%12.8f %12.8f\n', vgeom(i,2), vgeom(i,1)');
end
fclose(uo);
