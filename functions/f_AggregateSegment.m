function [segs, hydcase] = f_AggregateSegment(segs,hydcase,y,seg2erase)
% Erases a specified or the least relevant segment of an event. The segment is then merged with its neighbors
% Uwe Ehret, 15.Nov.2013, modified: Simon Seibert March 3rd 2014

% INPUT
    % segs: (1,x) struct, where x is the number of segments found in the entire event
    % hydcase: (n,1) matrix with hydrological case: -2=valley -1=drop, 1=rise 2=peak  
    % seg2erase: optional, number of the segment to erase
% OUTPUT
    % segs: struct with reduced number of segments and adjusted segement properties 
    % hydcase: time series with adjusted hydrologlical cases
% METHOD
    % Note: The first or the last segment can ONLY erased if only two segments are left. This assures that obs and sim event both 
    % - still start with the same hydcase (rise or fall)
    % - end with the same hydcase (rise or fall)

% identify the least relevant segment if the segment to erase is not specified
if ~exist('seg2erase', 'var')   % the segment to erase is not specified, 
    if length(segs) > 2         % if more than the first and last segment are left ...
        dum = [segs.relevance];    
        dum = dum(2:end-1);     % exclude the first and last from deletion
        seg2erase =  find (dum == min(dum));    % find the least relevant segment
        seg2erase = seg2erase(1)+1; % consider the first segement only if several segements of identical relevance were found (add 1 to the position, as we deleted the first (and last) segment before  
    elseif length(segs) == 2    % only the first and last of the segments are left
        dum = [segs.relevance];     
        seg2erase =  find (dum == min(dum)); 
        seg2erase = seg2erase(1); % reduce to size 1 in case several equally relevant segments were found  
    else                        % less than 2 segments left
        error('f_AggregateSegment: less than 2 segments left!');
    end
else                            % the segment to erase is specified
    if length(segs) > 2         % if more than the first and last segment are left ...
        if (seg2erase == 1 || seg2erase == length(segs))  % if more than the first and last segment are left ...
            error('f_AggregateSegment: more than 2 segments left: cannot erase the first or last segment!');
        end
    elseif length(segs) == 2    % only the first and last of the segments are left 
        % no problem
    else                        % less than 2 segments left
        error('f_AggregateSegment: less than 2 segments left!');
    end
    
end

%% join and modify segment properties

% delete the first segment and join it with the second
if seg2erase == 1  
    % adjust the hydrological cases (take over the value of the following segment)
    hydcase(segs(seg2erase).starttime_local:segs(seg2erase).endtime_local) = hydcase(segs(seg2erase).starttime_local + 1);   
    
    % join 2 segments: the segment to erase, and the following. the following takes it all
    segs(seg2erase+1).starttime_local = segs(seg2erase).starttime_local;
    segs(seg2erase+1).starttime_global = segs(seg2erase).starttime_global;
    segs(seg2erase+1).length = segs(seg2erase+1).endtime_local - segs(seg2erase).starttime_local;
    segs(seg2erase+1).rel_length = (segs(seg2erase+1).endtime_local - segs(seg2erase).starttime_local) / (length(y)-1);
    segs(seg2erase+1).sum_dQ = sum(diff(y(segs(seg2erase).starttime_local : segs(seg2erase+1).endtime_local)));
    segs(seg2erase+1).rel_dQ = sum(abs(diff(y(segs(seg2erase).starttime_local : segs(seg2erase+1).endtime_local))))/sum(abs(diff(y))); 
    segs(seg2erase+1).relevance = segs(seg2erase+1).relevance + segs(seg2erase).relevance;
    
    % erase the first segment
    segs(seg2erase)=[];    
    
% delete the last segment and join it with the secondlast
elseif seg2erase == length(segs)    
    % adjust the hydrological cases (take over the value of the previous segment)
    hydcase(segs(seg2erase).starttime_local:segs(seg2erase).endtime_local) = hydcase(segs(seg2erase).starttime_local - 1); 
    
    % join 2 segments: the segment to erase, and the previous. the previous takes it all
    segs(seg2erase-1).endtime_local = segs(seg2erase).endtime_local;
    segs(seg2erase-1).endtime_global = segs(seg2erase).endtime_global;
    segs(seg2erase-1).length = segs(seg2erase).endtime_local - segs(seg2erase-1).starttime_local;
    segs(seg2erase-1).rel_length = (segs(seg2erase).endtime_local - segs(seg2erase-1).starttime_local) / (length(y)-1);
    segs(seg2erase-1).sum_dQ = sum(diff(y(segs(seg2erase-1).starttime_local : segs(seg2erase).endtime_local)));
    segs(seg2erase-1).rel_dQ = sum(abs(diff(y(segs(seg2erase-1).starttime_local : segs(seg2erase).endtime_local))))/sum(abs(diff(y))); 
    segs(seg2erase-1).relevance = segs(seg2erase-1).relevance + segs(seg2erase).relevance;
   
    % erase the last segment
    segs(seg2erase)=[]; 

% 3 or more segments are left, join all three of them and update segement properties
else  
    % adjust the hydrological cases (take over the value of the previous segment)
    hydcase(segs(seg2erase).starttime_local:segs(seg2erase).endtime_local) = hydcase(segs(seg2erase).starttime_local - 1); 

    % join 3 segments: the segment to erase, the previous and the following. the previous takes it all
    segs(seg2erase-1).endtime_local = segs(seg2erase+1).endtime_local;
    segs(seg2erase-1).endtime_global = segs(seg2erase+1).endtime_global;  
    segs(seg2erase-1).length = segs(seg2erase+1).endtime_local - segs(seg2erase-1).starttime_local;
    segs(seg2erase-1).rel_length = (segs(seg2erase+1).endtime_local - segs(seg2erase-1).starttime_local) / (length(y)-1);
    segs(seg2erase-1).sum_dQ = sum(diff(y(segs(seg2erase-1).starttime_local : segs(seg2erase+1).endtime_local)));
    segs(seg2erase-1).rel_dQ = sum(abs(diff(y(segs(seg2erase-1).starttime_local : segs(seg2erase+1).endtime_local))))/sum(abs(diff(y))); 
    segs(seg2erase-1).relevance =  segs(seg2erase-1).relevance + segs(seg2erase).relevance + segs(seg2erase+1).relevance; 
    
    % erase the segment to erase and the following
    % NOTE: here the order of deletion is important! Work from the end backwards, otherwise enumeration gets screwed up!
    segs(seg2erase+1)=[];
    segs(seg2erase)=[];   
end

end






