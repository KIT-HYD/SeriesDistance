function [vals,count] = f_ReplaceEqualNeighbours(vals)
% Replaces equal neighbouring values in a series
% Uwe Ehret, 15.Nov.2013

% INPUT
    % vals: (n,1) matrix with values
% OUTPUT
    % vals: (n,1) matrix, same as input, but equal neighbors replaced
    % count: number of equal neighbors found
% METHOD
    % If a series of equal neigbouring values is found, they are successively raised by 1/1000 
    % --> in a sequence of equal values, the last will become the largest

count = 0;
for z=1:length(vals)-1; % loop over all values
    if vals(z) == vals(z+1) % Equal neighboriung values? 
        count = 1;

        while (z+count)<=length(vals) && vals(z+count)==vals(z) % Loop over all neighboring equal values
            vals(z+count)=vals(z+count-1) + 0.0001;   % raise the value of the previous point by 0.0001
            count = count +1;
        end
    end
end

end

