function [hydcase] = f_CalcHydCase(vals)
% Returns the hydrological case for each timestep of a hydrological timeseries
% Uwe Ehret, 15.Nov.2013

% INPUT
    % vals: (n,1) matrix with values
% OUTPUT
    % hydcase: (n,1) matrix with hydrological case: -2=valley -1=drop, 1=rise 2=peak
% METHOD
    % for each point, calcultes the gradient to the previous and the next value
    % drop-rise: valley   drop-drop: drop   rise-rise: rise   rise-drop : peak  

len = length(vals);
hydcase = NaN(len,1);   % initialize result array

% Find hydrological case for first and last timestep (special case due to incomplete neighborhood)
if (vals(2)-vals(1)) < 0
    hydcase(1)=-1;
else
    hydcase(1)=1;
end
if (vals(end)-vals(end-1)) < 0
    hydcase(end)=-1;
else
    hydcase(end)=1;
end

% Find hydrological case for each timestep
for z=2:len-1   % loop over all values except the first and last
    if (vals(z)-vals(z-1))<0 && vals(z+1)-vals(z)>0, hydcase(z)= -2;end;  %drop-rise: valley
    if (vals(z)-vals(z-1))<0 && vals(z+1)-vals(z)<0, hydcase(z)= -1;end;  %drop-drop: drop
    if (vals(z)-vals(z-1))>0 && vals(z+1)-vals(z)>0, hydcase(z)= 1;end;  %rise-rise: rise
    if (vals(z)-vals(z-1))>0 && vals(z+1)-vals(z)<0, hydcase(z)= 2;end;  %rise-drop : peak  
end

end

