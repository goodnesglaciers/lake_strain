%SlipEst.m
%  Script to Estimate Slip on Dislocations and Volume changes for volcanic
%  sources for Fixed Fault Geometry.
%  Paul Segall and Roland Burgmann 1995


%Poisson's ratio
        nu = 0.25;

%Update fault geometry and station coordinates in local coords
        UpdateGeom

% CREATE DESIGN MATRIX

% GPS DATA
if length(d_gps) > 0
        xy_gps = llh2localxy(llh_gps, origin);

	if nf > 0
        	G_dis = MakeGgps(dis_geom, xy_gps, nu);
	else
		G_dis = [];
	end
	if nv > 0
		G_vol = MakeGvolgps(vol_geom, xy_gps, nu);
	else
		G_vol = [];
	end


	G = [G_dis, G_vol];
else
	G = [];
end 

% LEVELING DATA

if length(d_lev) > 0
        xy_lev = llh2localxy(llh_lev, origin);

	if nf > 0
		clear G_dis;
        	G_dis = MakeGlev(dis_geom, xy_lev, nu);
	else
		G_dis = [];
	end
	if nv > 0
		clear G_vol;
		G_vol = MakeGvollev(vol_geom, xy_lev, nu);
	else
		G_vol = [];
	end


	G = [G; G_dis, G_vol];

end 


% invert and estimate slip
	Ginv = pinv(W*G);
	slip = Ginv*W*d;
	cov_slip = Ginv*Ginv';
	sig_slip = sqrt(diag(cov_slip));
	dhat = G*slip;
        r = W*(d - dhat);


% Output results
	PrintSlip(slip, sig_slip, nf, nv);
	PrintResids(r, d_gps, d_lev, G);


% Clear Plot Files
	!rm temp.gmtlin
	!rm volums.gmtlin
	!rm obs.gmtvec
	!rm mod.gmtvec

% IF horizontal data exists -- plot  vectors

if length(d_gps) > 0
	%plot observed vectors
        [a,b,az] = ErrlpsGMT( cov3d2cov2d(cov_gps)  );
	nvec = length(d_gps)/3;
        outmatrix1 = [llh_gps(2,1:nvec)', llh_gps(1,1:nvec)', vel(1,:)', ...
		vel(2,:)',a',b',az'];

	% plot predicted vectors
	u_pre = zeros(3,nvec);
	for i = 1:nvec
		u_pre(:,i)  = dhat( 3*(i-1)+1 : 3*(i-1)+3 );
	end
	null = zeros(size(u_pre(1,:)))';
	outmatrix2 = [llh_gps(2,1:nvec)', llh_gps(1,1:nvec)', u_pre(1,:)', ...
		 u_pre(2,:)', null, null, null];

	save obs.gmtvec outmatrix1 -ascii
	save mod.gmtvec outmatrix2 -ascii
end

% Plot Fault Outlines and Volume Sources in GMT plot

	if nf > 0
		for i=1:nf
	  	 displot_gmt( geom(i,:), dis_geom(i,:) ) 
		end
	end
	
	if nv >0
	 	 volplot_gmt( vgeom ); 
	end

% Run Script
	if get(cb_plt, 'value') == 1
		!test.map
	end


% IF Leveling data exists plot observed and predicted

if length(d_lev) > 0
 
	% plot observed level data
        figure
        errorbar (Len,d_lev, sigma_lev, sigma_lev, 'y+')
        hold on
 
	% plot modeled level data
        plot (Len,dhat(length(d_gps)+1:length(d_gps)+length(d_lev)), 'g')
        hold off;
end

disp('Done')

