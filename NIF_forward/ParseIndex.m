function masterindex=ParseIndex(siteindex,epochindex,sites,epochs,SITEINDEX,EPOCHINDEX)

    NS=size(sites,1);
    NE=size(epochs,1);

%Determine which stations are indexed

    switch class(siteindex)

        case 'double'   %Explicit index reference

            if any(siteindex>NS)
                error('Index exceeds number of stations.')
            end

        case {'char','cell'}    %Colon, list of stations (in cell array or char array form)

            siteindex=cellstr(siteindex);
            if strcmp(siteindex,':')
                siteindex=(1:NS)';
            else
                [NULL,NULL,siteindex]=intersect(siteindex,sites);
            end
            if isempty(siteindex)
                error('Could not match any of the specified stations to the masterlist.')
            end

    end

%Determine which epochs are indexed

    switch class(epochindex)

        case 'double'   %Explicit index reference
            
            if any(epochindex>NE)
                error('Index exceeds number of epochs.')
            end

        case 'cell'     %Range of times

            epochindex=find(epochs>=epochindex{1} & epochs<=epochindex{2});
            if isempty(epochindex)
                error('Index falls outside epoch range.')
            end

        case 'char'     %The colon

            if epochindex==':'
                epochindex=(1:NE)';
            else
                error('Cannot parse epochs index.')
            end

    end

%Find master logical index

    SI=zeros(NS,1);
    EI=zeros(NE,1);
    SI(siteindex)=1;
    EI(epochindex)=1;
    masterindex=SI(SITEINDEX)& EI(EPOCHINDEX);
