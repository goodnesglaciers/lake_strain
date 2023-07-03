	function [rate, rate_sig] = rate_recon(x, sigma, Nbasis)

%%
%% function [rate, rate_sig] = rate_recon(x, sigma, Nbasis)
%%
%%  Reconstruct slip-rate and uncertainty at each epoch from
%%  state vector and covariance
%%

%%  Determine some constants
	[Nepochs, statedim] = size(x);
	Nsites = statedim - 3*Nbasis;

        rate = zeros(Nbasis,Nepochs);
        Tr   = zeros(Nbasis, 3*Nbasis);
        tmpvar = zeros(Nbasis,Nbasis);
        rate_sig = zeros(Nbasis,Nepochs);
        sig = eye(statedim);

for k = 1:Nepochs
        for i = 1:Nbasis
                Tr(i,3*(i-1)+1:3*(i-1)+3) = [1,0,1];
        end

        rate(:,k) = Tr*x(k,1:3*Nbasis)';

	%recunstruct covariance from vectorized form       
	sig(:) = sigma(k,:);     
        tmpvar = Tr*sig(1:3*Nbasis,1:3*Nbasis)*Tr';
        rate_sig(:,k) = sqrt(diag(tmpvar));
end
