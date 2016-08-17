function [segs] = f_DefineSegments(x,y)
% Defines all segments in an event (starttime, endtime, realtive duration, relative magnitude change (dQ fraction), importance)
% Uwe Ehret, 15.Nov.2013, modified Simon Seibert, March 2014

% INPUT
    % x: (n,1) matrix with time position (x-position) of values
    % y: (n,1) matrix with values
% OUTPUT
    % segs: (1,x) struct, where x is the number of segments found in the entire event
    % Note: A segment always includes its first and last point (start, valley, peak or end) --> peaks and valleys are used twice!

% initialize structure
segs = struct('starttime_local',{},'endtime_local',{},'starttime_global',{},'endtime_global',{}, 'length',{}, 'rel_length',{},'sum_dQ',{}, 'rel_dQ',{},'relevance',{});

% find all segments
segs(1).starttime_global = x(1);
segs(1).starttime_local = 1;

for z = 2 : length(y)-1 % loop over all values except the first and last
    hydcase = f_CalcHydCase(y);  
    if hydcase(z)==2 || hydcase(z)==-2  % peak or valley
        segs(end).endtime_global = x(z);
        segs(end).endtime_local = z;
        segs(end+1).starttime_global = x(z); % (end+1) is on purpose! It creates a new segment template in the structure
        segs(end).starttime_local = z;          % (end) is on puropse! As a new segment has been created in the previous line, the entry is now made in the new segment
    end
end
segs(end).endtime_global = x(end);
segs(end).endtime_local = length(y);

% compute segment properties
for z= 1 : length (segs)
    % length of the segment
    segs(z).length = segs(z).endtime_local - segs(z).starttime_local; 
    % relative length to the entire time series [0,1]
    segs(z).rel_length = (segs(z).endtime_local - segs(z).starttime_local) / (length(y)-1); 
    % sum of segment slopes
    segs(z).sum_dQ = sum(diff(y(segs(z).starttime_local:segs(z).endtime_local)));
    % sum of slopes relative to entire event [0,1]
    segs(z).rel_dQ = (sum(abs(diff(y(segs(z).starttime_local:segs(z).endtime_local)))))/sum(abs(diff(y)));   % discharge changes relative to entire event [0,1]
    % relative importance of the segement
    segs(z).relevance = sqrt(segs(z).rel_length^2 + segs(z).rel_dQ^2); % relevance as the euclidean distance of rel_length and rel_dQ
end

% normalize the segment relevance with overall relevance of the entire event [0,1]
Relevance = sum([segs.relevance]);
for z = 1:length(segs)
    segs(z).relevance= segs(z).relevance/Relevance;
end