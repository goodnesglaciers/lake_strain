	function u_rel = relative_disp(nu, dis_geom, xy)

	[nsta,m] = size(xy);
	nvec = nsta-1;

	u = zeros(nsta,3);
	for i= 1:nsta
		u(i,:) = disloc(nu, dis_geom, [xy(i,1), xy(i,2)]); 
	end

%  compute  displacements realtive to LAST station
 	  u_rel = zeros(3,nvec);
	  u_rel(1,:) = u(1:nvec,1)' - u(nvec+1,1).*ones(1,nvec);
	  u_rel(2,:) = u(1:nvec,2)' - u(nvec+1,2).*ones(1,nvec);
	  u_rel(3,:) = u(1:nvec,3)' - u(nvec+1,3).*ones(1,nvec);

