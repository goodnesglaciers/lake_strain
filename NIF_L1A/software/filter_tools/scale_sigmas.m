function [E_Sig, N_Sig, U_Sig] = scale_sigmas(t, E_Sig, N_Sig, U_Sig, sigma)

%% function [E_Sig, N_Sig, U_Sig] = scale_sigmas(t, E_Sig, N_Sig, U_Sig, sigma)

%% Scale the uncertainties to reflect the fact that 
%% Bernese software underestimates the errors
%% This scaling is based on JGR results

%% Also scale by sigma, the rms of the solution

for i = 1:length(t)
        if t(i) < 91.4603               %scaling based on 
              E_Sig(:,i) = sqrt(6)*E_Sig(:,i); 
              N_Sig(:,i) = sqrt(6)*N_Sig(:,i); 
              U_Sig(:,i) = sqrt(6)*U_Sig(:,i); 
        else
              E_Sig(:,i) = sqrt(10)*E_Sig(:,i); 
              N_Sig(:,i) = sqrt(10)*N_Sig(:,i);  
              U_Sig(:,i) = sqrt(10)*U_Sig(:,i); 
        end
end

	E_Sig = sigma*E_Sig;
	N_Sig = sigma*N_Sig;
	U_Sig = sigma*U_Sig;
