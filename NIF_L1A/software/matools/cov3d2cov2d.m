	function [cov2d] = cov3d2cov2d(cov3d)

	nvec = length(diag(cov3d))/3;


	for i=1:nvec 
		for j = 1:nvec
			for k = 1:2
			cov2d((i-1)*2 + 1, (j-1)*2 + k) = cov3d((i-1)*3 + 1, (j-1)*3 +k);
			end 
			for k = 1:2
			cov2d((i-1)*2 + 2, (j-1)*2 + k) = cov3d((i-1)*3 + 2, (j-1)*3 +k);
			end 
		end 
	end  

