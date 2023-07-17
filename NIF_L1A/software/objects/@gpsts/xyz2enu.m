function tsenu=xyz2enu(tsxyz,origin)
%XYZ2ENU    ts_enu=xyz2enu(ts_xyz,origin)
%
%Transforms from a gpsts in global cartestian (XYZ) coordinates
%to a local coordinate system aligned with the geographic directions at
%the position specified by 'origin'.
%
%Input 'tsxyz' is a gpsts object.  Input 'origin' should
%be a vector of length 2 or 3.  In the former case, the function
%assumes the origin is specified as a longitude, latitude pair (degrees),
%while in the latter case the origin is assumed to a cartesian (XYZ)
%triple.  
%
%Output 'ts_enu' is a gpsts object containing the transformed coordinates. 
%
%-------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Apr 14, 2001   Peter Cervelli		 Original Code
%
%-------------------------------------------------------------

tsenu=tsxyz;

[tsenu.d,tsenu.dcov]=xyz2enu(tsenu.d,tsenu.dcov,origin);

% the following doesn't work anymore 11/2005, can't deal with blank covariance matrix
%[tsenu.apcoords]=xyz2enu(tsenu.apcoords,origin);
[tsenu.apcoords]=xyz2enu_vect(tsenu.apcoords(:),origin);
