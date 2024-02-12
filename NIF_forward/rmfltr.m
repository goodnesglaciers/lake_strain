function xavg = rmfltr ( x, npts, lng, ndec )

% function xavg = rmfltr ( x, npts, lng, ndec )
%
% Running mean low pass filter of data vector 'x', length 'npts', with 
% filter of length 'lng' and decimation factor 'ndec'.  Filtered values 
% start at xavg[1] regardless of the value of 'lng'.  For 'ndec' = 1 the 
% number of averaged points is (npts-(lng-1)).  For ndec>1 averaged data 
% is decimated, returning every 'ndec' point after the first point.  The 
% returned vector is of length 'nleft', the number of points left after 
% averaging and decimation.
%
% AJP, 11-14-94

% check input parameters
if ( npts <= 0 | lng <= 0 | ndec <= 0)
  fprintf ('   rmfltr error: illegal parameter values \n');
  fprintf ('     npts= %d, lng= %d, ndec= %d \n', npts, lng, ndec);
  return;
end
if( lng > npts ) 
  fprintf ('   rmfltr error: lng= %d > npts= %d \n', lng, npts);
  return;
end 
if( ndec > lng ) 
  fprintf('   rmfltr warning: ndec= %d > lng= %d \n', ndec, lng);
end

% compute points left after averaging, average the array
nleft = npts - (lng - 1);
for i = 1:nleft
  x(i) = mean( x(i:(i+lng-1)) );  
end

% decimate if necessary, reset number points left
% use definition of modulo(a,b) = a - ( fix(a/b)*b )
if (ndec > 1)
  j = 1;
  for i = 2:nleft
    mod = (i-1) - ( fix((i-1)/ndec) * ndec );
    if( mod == 0 )
      j = j + 1;
      x(j) = x(i);
    end
  end
  nleft = j;
end

%fprintf ('   rmfltr: lng= %d, ndec= %d, nleft= %d \n', lng,ndec,nleft);
xavg = x(1:nleft);
return
