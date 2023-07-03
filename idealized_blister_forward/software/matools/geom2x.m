function x = geom2x(geom, origin)
%	Input dislocation geometry and output vector of unknowns
%		for optimization.
%	Modified so that slip magnitude is not included in x.
%	x -- model vector
%	geom -- input geometry vector
	

%    Input: Dislocation Geometry
%               geom(1) = wid = fault width in dip direction (km)
%               geom(2) = dep = depth of fault bottom (km)
%               geom(3) = dip = dip angle (degrees)
%               geom(4) = latitude of one end (decimal degrees)
%               geom(5) = longitude of one end (decimal degrees)
%               geom(6) = latitude of second end (decimal degrees)
%               geom(7) = longitude of second end (decimal degrees)
%               geom(8) = srike slip (m)
%               geom(9) = dip slip (m)
%               geom(10) = opening slip (m)

[nf,nel] = size(geom);
tempgeom = zeros(nf,7);

for i=1:nf
        xy_1end(i,:)  = llh2localxy([geom(i,4);  geom(i,5)], origin);
       	xy_2end(i,:)  = llh2localxy([geom(i,6);  geom(i,7)], origin);
        len(i) = sqrt(  (xy_1end(i,1) -  xy_2end(i,1) )^2 + ...
			(xy_1end(i,2) - xy_2end(i,2))^2);
       	strik(i) = (180/pi)*atan2((xy_2end(i,2) -  xy_1end(i,2) ), ...
			(xy_2end(i,1) -  xy_1end(i,1) ));
	strik(i) = 90.0 - strik(i);
	midpointlat(i)  = 0.5*(geom(i,4) + geom(i,6));
	midpointlong(i) = 0.5*(geom(i,5) + geom(i,7));

	tempgeom(i,:) = [len(i),geom(i,1),geom(i,2),geom(i,3), strik(i),...
		midpointlat(i),midpointlong(i)];
end

	geomt = tempgeom';
	x = geomt(:);
