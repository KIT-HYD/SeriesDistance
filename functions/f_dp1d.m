function [xy_dp] = f_dp1d(xy,tol,numpoints,nselimit)
% Line simplification with a modified Douglas-Peucker algorithm
% Uwe Ehret, 15.Nov.2013

% INPUT
    % xy:           (n,2) matrix with n x and y data of a line (x in ascending order)
    % tol:          maximum tolerance (in vertical, i.e. y-direction) allowed between the original and simplified series. In units of the y-values
    % numpoints:    optional: states by how many points the line shall be approximated. If numpoints is set, tol is not used
    % nselimit:     optional: states which Nash-Sutcliffe efficiency should be achieved between the original and simplified line. 
    %               If nselimit is set, tol and numpoints are not used

% OUTPUT
    % xy_dp:       (m,2) vector with m x and y data of the simplified line (x in ascending order) m <= n
    %               xy_dp(:,1) --> x-values of the simplified series
    %               xy_dp(:,2) --> y-values of the simplified series
% METHOD
    % this function approximates a line by fewer of its points
    % the simplified line starts with the first and the last point of xy, then the points of maximum distance are subsequently added until
    % the selected cutoff criterion is met (tol, numpoints, or nselimit)

% check the input
if exist('nselimit', 'var')
    calccase = 3;
elseif exist('numpoints', 'var')
        calccase = 2;
else
    calccase = 1;
end
   
% initialize the simplified line by the first and last value
xy_dp = [];
xy_dp(1,1) = xy(1,1);
xy_dp(1,2)= xy(1,2);
xy_dp(2,1) = xy(end,1);
xy_dp(2,2)= xy(end,2);

alldone = 0;
switch calccase
    case 1  % max tol limit
        while alldone == 0  % loop until the stop criterion is fulfilled
            intp_dp=interp1(xy_dp(:,1),xy_dp(:,2),xy(:,1),'linear');    % sample the simplified line at all x-locations of the original line

            % find the maximum distance between simplied and original line
            d = abs(xy(:,2)-intp_dp);   % determine the absolute distance at each x-location of the original line (vertical, i.e. in y-direction)
            index = find(d == max(d));  % find the maximum of the distances ...
            index = index(1);           % ... if there are several maxima of the same value, take only the first

            if d(index)> tol            % if the max distance exceeds the tolerance ...
                x_add = xy(index,1);    % find the x-value of the maximum distance
                y_add = xy(index,2);    % find the y-value of the maximum distance
                if x_add > xy_dp(end,1) % if the maximum value was at the very end of the line (should not, as we take the first and last as initial values)
                    xy_dp = cat (1,xy_dp,[x_add y_add]);   
                else
                    insertpos = find(xy_dp(:,1)>x_add ,1,'first');  % find the position where to insert the new point
                    xy_dp = cat (1,xy_dp(1:insertpos-1,:),[x_add y_add],xy_dp(insertpos:end,:)); % add the new point to the simplified line  
                end
            else                    % if the max distance is smaller than the tolerance
                alldone = 1;        % we're done
            end      
        end     
        
    case 2  % numpoints limit
        while alldone == 0  % loop until the stop criterion is fulfilled
            intp_dp=interp1(xy_dp(:,1),xy_dp(:,2),xy(:,1),'linear');    % sample the simplified line at all x-locations of the original line

            % find the maximum distance between simplied and original line
            d = abs(xy(:,2)-intp_dp);   % determine the absolute distance at each x-location of the original line (vertical, i.e. in y-direction)
            index = find(d == max(d));  % find the maximum of the distances ...
            index = index(1);           % ... if there are several maxima of the same value, take only the first

           if length(xy_dp) < (numpoints*2) % if the number of points is not enough ...
                x_add = xy(index,1);    % find the x-value of the maximum distance
                y_add = xy(index,2);    % find the y-value of the maximum distance
                if x_add > xy_dp(end,1) % if the maximum value was at the very end of the line (should not, as we take the first and last as initial values)
                    xy_dp = cat (1,xy_dp,[x_add y_add]);   
                else
                    insertpos = find(xy_dp(:,1)>x_add ,1,'first');  % find the position where to insert the new point
                    xy_dp = cat (1,xy_dp(1:insertpos-1,:),[x_add y_add],xy_dp(insertpos:end,:)); % add the new point to the simplified line  
                end
            else                    % if the the simplified line consists of numpoints
                alldone = 1;        % we're done
           end
        end
        
    case 3  % nse limit
        while alldone == 0  % loop until the stop criterion is fulfilled
            intp_dp=interp1(xy_dp(:,1),xy_dp(:,2),xy(:,1),'linear');    % sample the simplified line at all x-locations of the original line

            % find the maximum distance between simplied and original line
            d = abs(xy(:,2)-intp_dp);   % determine the absolute distance at each x-location of the original line (vertical, i.e. in y-direction)
            index = find(d == max(d));  % find the maximum of the distances ...
            index = index(1);           % ... if there are several maxima of the same value, take only the first

            % compute the NSE
            A = xy(:,1);
            B = intp_dp;
            C = cat(2,A,B);
           [nse metric_id] = f_nashsutcliffe(xy,C);
           
           if nse < nselimit % the nse is not good enough yet ...
                x_add = xy(index,1);    % find the x-value of the maximum distance
                y_add = xy(index,2);    % find the y-value of the maximum distance
                if x_add > xy_dp(end,1) % if the maximum value was at the very end of the line (should not, as we take the first and last as initial values)
                    xy_dp = cat (1,xy_dp,[x_add y_add]);   
                else
                    insertpos = find(xy_dp(:,1)>x_add ,1,'first');  % find the position where to insert the new point
                    xy_dp = cat (1,xy_dp(1:insertpos-1,:),[x_add y_add],xy_dp(insertpos:end,:)); % add the new point to the simplified line  
                end
            else                    % if the nse is good enough
                alldone = 1;        % we're done
           end
        end           
        
end






