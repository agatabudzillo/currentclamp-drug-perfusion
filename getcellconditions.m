function [ conditions, numconditions ] = getcellconditions( numfiles )
%getcellconditions asks the user to define cell conditions and filename
%starts for those conditions

numconditions = input('Enter number of cell conditions in folder: ');

conditions = cell(numconditions,2);

if numconditions == 1;
    conditions{1,1} = input('Enter name of condition: ','s');
    conditions{1,2} = 0;
else
    for c = 1:numconditions
        conditions{c,1} = input(['Enter name of condition ' num2str(c) ': '],'s');
        conditions{c,2} = input(['Enter last three digits of the first file for ' '"' conditions{c,1} '"' ': ']);
    end
end


end

