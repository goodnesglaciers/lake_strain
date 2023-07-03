	function [slip, slip_sig] = slip_recon(t, x, sigma, Nbasis)

%%
%% function [slip, slip_sig] = slip_recon(x, sigma)
%%
%%  Reconstruct slip and uncertainty at each epoch from
%%  state vector and covariance
%%

%%  Determine some constants
	[Nepochs, statedim] = size(x);
	Nsites = statedim - 3*Nbasis;

        slip = zeros(Nbasis,Nepochs);
        Ts   = zeros(Nbasis, 3*Nbasis);
        tmpvar = zeros(Nbasis,Nbasis);
        slip_sig = zeros(Nbasis,Nepochs);
        sig = eye(statedim);

for k = 1:Nepochs
        for i = 1:Nbasis
                Ts(i,3*(i-1)+1:3*(i-1)+3) = [ t(k)-t(1),1,0];
        end

        slip(:,k) = Ts*x(k,1:3*Nbasis)';

	%recunstruct covariance from vecotrized form       
	sig(:) = sigma(k,:);     
        tmpvar = Ts*sig(1:3*Nbasis,1:3*Nbasis)*Ts';
        slip_sig(:,k) = sqrt(diag(tmpvar));
end
