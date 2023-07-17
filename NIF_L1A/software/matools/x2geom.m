function [geom, dis_geom] = x2geom(x, origin)
%	Given the model vector x, return the geometry and dis_geom arrays
%	x 	-- model vector
%	geom 	--  geometry array
%	dis_geom --  geometry array
	

nf = length(x)/7;
tempgeom = zeros(nf,7);geom = zeros(nf,7);dis_geom = zeros(nf,7);

% N = radius of earth, e = eccentricity
        N = 6378206.4;
        e = 1./298.257;

for i=1:nf
	first = 7*(i-1)+1;
	tempgeom(i,:) = x(first:first+6)';  
	theta = (90 - tempgeom(i,5))*pi/180;

	EastOffset = 0.5*tempgeom(i,1)*cos(theta);
	NorthOffset = 0.5*tempgeom(i,1)*sin(theta);

	  	 W = sqrt(1 -e^2*sin(tempgeom(i,6)*pi/180));
       		 Rp = (N/W)*cos(tempgeom(i,6)*pi/180);
       		 Rm = N*(1-e^2)/W^3;

	North_Deg = (1.0e3*NorthOffset/Rm)*180/pi;
	East_Deg = (1.0e3*EastOffset/Rp)*180/pi;

	geom(i,1) = tempgeom(i,2); 	%fault width in dip direction (km)
	geom(i,2) = tempgeom(i,3);	%depth of fault bottom (km)
	geom(i,3) = tempgeom(i,4);	%dip angle (degrees)
	geom(i,4) = tempgeom(i,6) - North_Deg;
					%latitude of one end (decimal degrees)
	geom(i,5) = tempgeom(i,7) - East_Deg;
					%longitude of one end (decimal degrees)
	geom(i,6) = tempgeom(i,6) + North_Deg;
					%latitude of 2 end (decimal degrees)
	geom(i,7) = tempgeom(i,7) + East_Deg;
					%longitude of 2 end (decimal degrees)


        xy_1end(i,:)  = llh2localxy([geom(i,4);  geom(i,5)], origin);
       	xy_2end(i,:)  = llh2localxy([geom(i,6);  geom(i,7)], origin);

        delE(i) = 0.5*(  xy_1end(i,1) +  xy_2end(i,1) ) ;
        delN(i) = 0.5*(  xy_1end(i,2) +  xy_2end(i,2) ) ;

        dis_geom(i,:) = [tempgeom(i,1),geom(i,1),geom(i,2),geom(i,3),...
		tempgeom(i,5),delE(i),delN(i)];
end

