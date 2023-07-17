function [d2d, cov2d] =  3dto2d(d3d, cov3d)
%  [d2d, cov2d] =  3dto2d(d3d, cov3d)
%
% Input: 
%	d3d 	= 	three dimensional data vector
%	cov3d 	= 	three dimensional covriance matrix
% Output: 
%	d2d	=	two dimensional data vector
%	cov2d 	= 	two dimensional covriance matrix


	n = length(d3d)/3;
	d2d = zeros(1,2*n);
	cov2d = zeros(2*n);



      for i = 1:n
        d2d( (i-1)*2 + 1 ) = d3d( (i-1)*3 + 1);
        d2d( (i-1)*2 + 2 ) = d3d( (i-1)*3 + 2);

	for j = 1:n
      	  for k = 1:2
            cov2d((i-1)*2 + 1, (j-1)*2 + k) = cov3d((i-1)*3 + 1, (j-1)*3 +k);
          end 

      	  for k = 1:2
             cov2d((i-1)*2 + 2, (j-1)*2 + k) = cov3d((i-1)*3 + 2, (j-1)*3 +k);
          end 
        end  
      end  
