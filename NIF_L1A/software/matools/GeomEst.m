%GeomEst.m
%  Script to Estimate  Fault Geometry.
%  Currently implemented for GPS data only
%  Paul Segall 1995

%Poisson's ratio
        nu = 0.25;

        UpdateGeom
        xy = llh2localxy(llh, origin);

        d = vel(:);
	Icov = inv(cov);
        W = chol(Icov);

        ndata = length(d);
        nvec = ndata/3;

%save the starting model
       !rm temp.gmtlin
	[nf, nel] = size(geom);
        for i=1:nf
           displot_gmt( geom(i,:), dis_geom(i,:) )
        end
	!mv temp.gmtlin start.gmtlin
%

%Starting Value
	x0 = geom2x(geom, origin);

% Set Options for Optimization program

options(1) = 1;
options(2) = 1.0e-2;
options(3) = 1.0e-2;
%options(5) = 1;                        % use Gauss-newton
%options(9) = 1;
%options(14) = 100;
[x, options] = leastsq('Disl_wr',x0, options, [ ], origin, xy, nu, d, W, Icov, geom);

disp(' number of function evals'), options(10)
disp(' estimated model' ), x

[geom_est, dis_geom_est] = x2geom(x, origin);
dis_geom = [dis_geom_est, geom(:,8:10)];
geom     = [geom_est, dis_geom(:,8:10)];

        G = MakeGgps(dis_geom, xy, nu);
        Ginv = pinv(W*G);
        slip = Ginv*W*d;
        cov_slip = Ginv*Ginv';
        sig_slip = sqrt(diag(cov_slip));
        dhat = G*slip;
 
% Output results
        PrintSlip(slip, sig_slip);


% plot observed vectors
        [a,b,az] = ErrlpsGMT( cov3d2cov2d(cov)  );
        outmatrix1 = [llh(2,1:nvec)', llh(1,1:nvec)', vel(1,:)',...
			 vel(2,:)',a',b',az'];

% plot predicted vectors
        u_pre = zeros(3,nvec);
        for i = 1:nvec
                u_pre(:,i)  = dhat( 3*(i-1)+1 : 3*(i-1)+3 );
        end
        null = zeros(size(u_pre(1,:)))';
        outmatrix2 = [llh(2,1:nvec)', llh(1,1:nvec)', u_pre(1,:)',...
		u_pre(2,:)',null, null, null];

        !rm obs.gmtvec
        !rm mod.gmtvec
        save obs.gmtvec outmatrix1 -ascii
        save mod.gmtvec outmatrix2 -ascii

        for i=1:nf
           displot_gmt( geom(i,:), dis_geom(i,:) )
        end
if get(cb_plt, 'value') == 1
        !test.map
end

disp('Done')
