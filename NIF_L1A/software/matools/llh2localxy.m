	function [xy] = llh2localxy(llh,ll_org)

	[rows, nsta] = size(llh);
	% change from decimal degrees to decimal seconds
	lat = 3600.0 * llh(1,:);
	lon = 3600.0 * llh(2,:);
	
	Lat_Orig = 3600.0 * ll_org(1);
	Diff_long  = 3600.0 * ll_org(2)*ones(size(lon)) - lon;

	xy = zeros(nsta,2);
	for i=1:nsta
         xy(i,:) = polyconic(lat(i), Diff_long(i), Lat_Orig);
	end
 
	% convert units from meter into kilometer and flip x-axis
 
         xy(:,1) = -xy(:,1) / 1000.0;
         xy(:,2) =  xy(:,2) / 1000.0;
 
 
