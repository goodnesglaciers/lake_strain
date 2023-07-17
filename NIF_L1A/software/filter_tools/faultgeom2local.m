function dis_geom = faultgeom2local(geom, origin)


% dis_geom = faultgeom2local(geom, origin)

[nf, el] = size(geom);

for i=1:nf

        % make midpoint of first fault the origin
%        if i==1
%          origin = [0.5*(geom(i,4) + geom(i,6)); 0.5*(geom(i,5) + geom(i,7))];
%        end

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
