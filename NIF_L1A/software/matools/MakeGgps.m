function G = MakeGgps(dis_geom, xy, nu)

%Subroutine to make design matrix for slip

%Input:
%	dis_geom 	- Dislocation Geometry
%	xy		- xy coordinates of statinos
%	nu		- Poisson's ratio

	[nf,nel] = size(dis_geom);
	[nsta,nel2] = size(xy);
	ndata = (nsta-1)*3;  % this works for GPS data

	G = zeros(ndata,3*nf);

for i=1:nf
	if dis_geom(i,8)~=0 
	  u_rel_1 = relative_disp(nu, [dis_geom(i,1:7),1,0,0], xy); 
	  G(:,3*(i-1)+1) = u_rel_1(:);
	end

	if dis_geom(i,9)~=0 
	  u_rel_2 = relative_disp(nu, [dis_geom(i,1:7),0,1,0], xy); 
	  G(:,3*(i-1)+2) = u_rel_2(:);
	end

	if dis_geom(i,10)~=0 
	  u_rel_3 = relative_disp(nu, [dis_geom(i,1:7),0,0,1], xy); 
	  G(:,3*(i-1)+3) = u_rel_3(:);
	end
end

