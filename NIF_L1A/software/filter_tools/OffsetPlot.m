function [O, P] =  OffsetPlot(Obs, Pred)


        [Nobs,Nepochs] = size(Obs);

%% Offset lines so they don't overlap

% Displace lines by constant so they don't overlap
        O = zeros(size(Obs));
        P = zeros(size(Pred));
 
% first line doesn't change
        O(1,:) = Obs(1,:);
        P(1,:) = Pred(1,:);

for k=2:Nobs

        %% Determine the maximum offset


%        maxobs_offset = max( stripnan(Obs(k-1,:))) - min( stripnan(Obs(k,:)));
%        maxpre_offset = max( Pred(k-1,:)) - min( Pred(k,:));

        maxobs_offset = max( stripnan(O(k-1,:))) - min( stripnan(Obs(k,:)));
        maxpre_offset = max( P(k-1,:)) - min( Pred(k,:));


        off = max(maxobs_offset, maxpre_offset);
	if off>= 0
		offset = 1.2*off;
	else
		offset = 0.8*off;
	end

        O(k,:) = Obs(k,:)  + offset*ones(size(Obs(k,:)));
        P(k,:) = Pred(k,:) + offset*ones(size(Obs(k,:)));
end
