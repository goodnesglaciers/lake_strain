function output=swapmanager(i,j,X,swapdir)
%
%Manages the swap file for state objects

persistent stackindex
persistent stack

if isempty(stackindex)
    stackindex={'','','','','','',''};
end

switch class(X)

    case 'char'  %Reference

        filename=X;
        I=strmatch(filename,char(stackindex),'exact');
	%I
        if isempty(I)
            stackindex{7}=filename;
            stackindex(1)=[];

            A=struct2cell(load([swapdir,'\',filename]));

            stack{7}=A{1};
            stack(1)=[];
            I=6;
        end
        output=stack{I};
                

    case 'double'   %Assignment

        value=X;
        filename=['x_',num2str(i-1),'_',num2str(j-1)];
	%whos value filename
	%filename
        save([swapdir,'\',filename],'value')
        I=strmatch(filename,char(stackindex),'exact');
        if isempty(I)
            I=1;
        end

        stackindex{7}=filename;
        stackindex(I)=[];

        stack{7}=value;
        stack(I)=[];

        output=filename;

end