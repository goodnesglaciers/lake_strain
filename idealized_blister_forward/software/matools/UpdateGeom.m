% Update Fault Geometry


if nf==0 & nv>0
	disp('You have not defined any faults, this will be trouble')
	disp('as the coordinate origin is defined relative to the mid-point')
	disp('of the first fault')
	disp('This is  a bug that needs to be fixed')
	
	load junk
	%this just terminates the program
end 

for i=1:nf
	geom(i,1) = eval( get(getwid(i),'string') );
	geom(i,2) = eval( get(getdep(i),'string') );
	geom(i,3) = eval( get(getdip(i),'string') );
	geom(i,4) = eval( get(getlat1(i),'string') );
	geom(i,5) = eval( get(getlon1(i),'string') );
	geom(i,6) = eval( get(getlat2(i),'string') );
	geom(i,7) = eval( get(getlon2(i),'string') );
	geom(i,8) = eval( get(getsss(i),'string') );
	geom(i,9) = eval( get(getdss(i),'string') );
	geom(i,10) = eval( get(getops(i),'string') );

	% make midpoint of first fault the origin
	if i==1
	  origin = [0.5*(geom(i,4) + geom(i,6)); 0.5*(geom(i,5) + geom(i,7))];
       	end

	xy_1end(i,:)  = llh2localxy([geom(i,4);  geom(i,5)], origin);
	xy_2end(i,:)  = llh2localxy([geom(i,6);  geom(i,7)], origin);
	len(i) = sqrt(  (xy_1end(i,1) -  xy_2end(i,1) )^2 + ...
		(xy_1end(i,2) - xy_2end(i,2) )^2);
	strik(i) = (180/pi)*atan2((xy_2end(i,2) -  ...
		xy_1end(i,2) ), (xy_2end(i,1) -  xy_1end(i,1) )); 
	strik(i) = 90.0 - strik(i);

	delE(i) = 0.5*(  xy_1end(i,1) +  xy_2end(i,1) ) ;
	delN(i) = 0.5*(  xy_1end(i,2) +  xy_2end(i,2) ) ;
	dis_geom(i,:) =  [len(i),geom(i,1),geom(i,2),geom(i,3), strik(i),...
                delE(i),delN(i),geom(i,8), geom(i,9),geom(i,10)];
end


for i=1:nv
	vgeom(i,1) = eval( get(getlatv(i),'string') );
	vgeom(i,2) = eval( get(getlonv(i),'string') );
	vgeom(i,3) = eval( get(getdepv(i),'string') );
	vgeom(i,4) = eval( get(getvol(i),'string') );

	offset = llh2localxy([vgeom(i,1);  vgeom(i,2)], origin);
	vol_geom(i,1) = offset(1);
	vol_geom(i,2) = offset(2);
	vol_geom(i,3) = vgeom(i,3);
	vol_geom(i,4) = vgeom(i,4);
end
