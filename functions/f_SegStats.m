function[segment_stats] = f_SegStats(segs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    peaks = length(find([segs.sum_dQ]>0));
    troughs = length(find([segs.sum_dQ]<0));
    rise_dur = sum([segs(find([segs.sum_dQ]>0)).length]);
    fall_dur = sum([segs(find([segs.sum_dQ]<0)).length]);
    tot_dur = segs(end).endtime_global-segs(1).starttime_global;
    % summarize in a single array
    segment_stats = [peaks troughs rise_dur fall_dur tot_dur];
end

