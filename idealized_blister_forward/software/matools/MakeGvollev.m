function G = MakeGvollev(vol_geom, xy, nu)

%Subroutine to make design matrix for volume sources
% with leveling data

%Input:
%       vol_geom        - Source Geometry
%       xy              - xy coordinates of statinos
%       nu              - Poisson's ratio
 
[nv,nvels]  = size(vol_geom);
[nsta,nel2] = size(xy);
ndata = (nsta-1);  % appropriate for leveling data

G = zeros(ndata,nv);

% note relative_disp computes displacements relative to LAST station
% so put reference station last

	xy = [xy(2:nsta,:);xy(1,:)];


for i=1:nv
	u_rel = vol_rel_disp(nu, [vol_geom(i,1:3),1], xy);
        G(:,i) = u_rel(3,:)';
end
