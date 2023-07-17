function plotOPnew3(t, E,N,U, E_Sig, N_Sig, U_Sig, Ep, Np, Up, sites, psites, rel_site)


%% plotOPnew3(t, E,N,U, E_Sig, N_Sig, U_Sig, Ep, Np, Up, sites, psites, rel_site)
%% plot time series of observed and predicted station positions relative to
%% "rel_site".  Only the stations in the vector "psites" are ploted. Plot shows
%% predicted as line, and observed values with one standard deviation error bars
%%
%%   Inputs:
%%	t	=  vector of observation times.
%%	E,N,U   =  Observed values of East North and Up components
%%	E_sig,..=  One standard deviation of observed value in E, N, U
%%	Ep, Np,.=  Predicted values in E, N, U
%%	sites   =  vector of all site names
%% 	psites	=  vector of site names to be plotted
%%	rel_site=  name of station that all stations are plotted relative to

%%  Paul Segall
%%  Modified to plot +/- one sigma error bars

	[NPsites, four] = size(psites);
	[Ncomps, Nepochs] = size(E);


	E_Plot = []; E_pre_Plot = [];N_Plot = []; N_pre_Plot = [];
	Es_Plot = [];  Ns_Plot = [];

%% Pull out only the stations we want to plot 
	for i =1:NPsites
		k = GetIndex(sites, psites(i,:));
		E_Plot = [E_Plot;E(k,:)];
		E_pre_Plot = [E_pre_Plot;Ep(k,:)];
		Es_Plot = [Es_Plot;E_Sig(k,:)];
		N_Plot = [N_Plot;N(k,:)];
		N_pre_Plot = [N_pre_Plot;Np(k,:)];
		Ns_Plot = [Ns_Plot;N_Sig(k,:)];
	end

%%% change t into an array
	tm = [];
	for i = 1:NPsites
		tm = [tm;t];
	end
	

%%Plot East

	[O, P] =  OffsetPlot(E_Plot, E_pre_Plot);
	figure
	errorbar(tm', O', Es_Plot', 'o'), hold on
	plot( t, P, '--')
	title(['East Observed and Predicted Relative to ', rel_site]);
	xlabel('year')
	ylabel('East Displacement (m)')


%% Print the name of the station in the margin of the time series
	delt = 0.10*( t(length(t)) - t(1));
for i = 1:NPsites
	text(t(Nepochs)+delt, P(i, Nepochs), psites(i,:))
end

%%Plot North

	[O, P] =  OffsetPlot(N_Plot, N_pre_Plot);
	figure
	errorbar(tm', O', Ns_Plot', 'o'), hold on
	plot( t, P, '--')
	title(['North Observed and Predicted Relative to ',rel_site])
	xlabel('year')
	ylabel('North Displacement (m)')


%% Print the name of the station in the margin of the time series
	delt = 0.10*( t(length(t)) - t(1));
for i = 1:NPsites
	text(t(Nepochs)+delt, P(i, Nepochs), psites(i,:))
end


return

%%%%%% DONT DO FAULT PARALLEL ETC FOR NOW

%%try rotating into Fault Normal and Parallel

	strike_deg = -45;
	strike_rad = strike_deg*pi/180;
	FN = E_Plot*cos(strike_rad)  - N_Plot*sin(strike_rad);	
	FP = E_Plot*sin(strike_rad)  + N_Plot*cos(strike_rad);

	FN_pre = E_pre_Plot*cos(strike_rad)  - N_pre_Plot*sin(strike_rad);	
	FP_pre = E_pre_Plot*sin(strike_rad)  + N_pre_Plot*cos(strike_rad);

%%Plot Fault Normal

	[O, P] =  OffsetPlot(FN, FN_pre);
	figure
	plot(t, O, 'o'), hold on
	plot( t, P, '--')
	title(['Fault Normal Observed and Predicted Relative to ',rel_site])
	xlabel('year')
	ylabel('Fault Normal Displacement (m)')

for i = 1:NPsites
	text(t(Nepochs)+0.50, P(i, Nepochs), psites(i,:))
end


%%Plot Fault Parallel

	[O, P] =  OffsetPlot(FP, FP_pre);
	figure
	plot(t, O, 'o'), hold on
	plot( t, P, '--')
	title(['Fault Parallel Observed and Predicted Relative to ',rel_site])
	xlabel('year')
	ylabel('Fault Parallel Displacement (m)')

for i = 1:NPsites
	text(t(Nepochs)+0.50, P(i, Nepochs), psites(i,:))
end

