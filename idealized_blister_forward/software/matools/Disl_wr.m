function r = Disl_wr(x, origin, xy, nu, d, W, Icov, geom_start)
%
% Compute the weighted residual for Dislocation Model
%	x -- model vector
%	d -- data
%	W -- Weight Matrix
%	nu -- poissson's ratio
%	xy coordinates of stations
	
    	nf = length(x)/7;
   	ndata = length(d);
        nvec = ndata/3;
        u_rel = zeros(3,nvec);

	[geom_short, dis_geom_short] = x2geom(x, origin);
	dis_geom = [dis_geom_short, geom_start(:,8:10)];

%For Testing
%        figure
%        plot(xy(:,1), xy(:,2), 'o'), axis('equal'), hold on
%        for i=1:nf
%         displot(dis_geom(i,:))
%        end
%	 title('Geometry in Local Cartesian Coordiantes')

%        G = MakeGgps(dis_geom, xy, nu);

	clear G
	G = [];
for i=1:nf
        if dis_geom(i,8)~=0
          g = relative_disp(nu, [dis_geom(i,1:7),1,0,0], xy);
          G= [G, g(:)];
        end

        if dis_geom(i,9)~=0
          g = relative_disp(nu, [dis_geom(i,1:7),0,1,0], xy);
          G= [G, g(:)];
        end

        if dis_geom(i,10)~=0
          g = relative_disp(nu, [dis_geom(i,1:7),0,0,1], xy);
          G= [G, g(:)];
        end
end
%        dhat = G*pinv(W*G)*W*d;
	dhat = G*inv(G'*Icov*G)*G'*Icov*d;


% compute weighted residual sum of squares
        r = W*(d - dhat);

