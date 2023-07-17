	function [colorsnew] = setcolors;

	a = get(gca,'ColorOrder');
	colorsnew = a(2:6,:);

	set(gcf,'DefaultAxesColorOrder',colorsnew)

% NOTE: may have to execute this command outside the
% subroutinte to have it apply to all subsequent plots

