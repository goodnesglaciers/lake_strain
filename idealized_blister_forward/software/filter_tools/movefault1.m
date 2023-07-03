function outgeom = movefault1(ingeom)


%outgeom = movefault(ingeom);
% move the fault so that the coordinates of the midpoint refer to the
% fault bottom

% input ingeom
	%ingeom(1) = length
	%ingeom(2) = width
	%ingeom(3) = depth
	%ingeom(4) = dip
	%ingeom(5) = strike
	%ingeom(6) = East
	%ingeom(7) = North


% dip direction 
	dipdir = (180 + ingeom(:,5) + 90)*pi/180;
	offset = abs(ingeom(:,2).*cos(ingeom(:,4)*pi/180));
	EastOff = offset.*sin(dipdir);
	NorthOff = offset.*cos(dipdir);

	outgeom = ingeom;
	outgeom(:,3) = abs(ingeom(:,2).*sin(ingeom(:,4)*pi/180));
	outgeom(:,4) = 180 + ingeom(:,4);
	outgeom(:,6) = ingeom(:,6) + EastOff;
        outgeom(:,7) = ingeom(:,7) + NorthOff;
