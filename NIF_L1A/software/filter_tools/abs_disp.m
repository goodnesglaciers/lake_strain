	function u = abs_disp(nu, dis_geom, xy)

	[nsta,m] = size(xy);

	u = zeros(nsta,3);
	for i= 1:nsta
		u(i,:) = disloc(nu, dis_geom, [xy(i,1), xy(i,2)]); 
	end


	u = u';
