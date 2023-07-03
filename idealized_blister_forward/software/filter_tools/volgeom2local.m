function vol_geom = volgeom2local(vgeom, origin)
%
%
% 	vol_geom = volgeom2local(vgeom, origin)


	[nv, el] = size(vgeom);

for i=1:nv


        offset = llh2localxy([vgeom(i,1);  vgeom(i,2)], origin);

        vol_geom(i,1) = offset(1);
        vol_geom(i,2) = offset(2);
        vol_geom(i,3) = vgeom(i,3);
        vol_geom(i,4) = vgeom(i,4);



end



