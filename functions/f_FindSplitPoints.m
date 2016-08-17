function [ timeseries_splits ] = f_FindSplitPoints( obs, sim, split_frequency )
% method
% - overall, we want to split the time series at times where both obs and sim are in low flow
% - at each candidate split point, sample obs-sim pairs within a searchrange of
%   'perc' percent of the split frequency (forward and backward in time)
% - from each sample set, determine the ranks of the obs and sim values after sorting by size
% - the time step that has min(rank(obs) + rank(sim)) is the best split canditate
% - compare the obs and sim values of that time step with the cdf of the entire time series: 
%   In order to make a good split, the values should be also globally of low rank, expressed by paramater 'max_quantile' (probability of unexceedance)
%   If the obs and sim values are BOTH below the Pu limit of their time series, a split is done. If not, not. UNLESS ...
%   there are more than 'max_num_segs' segments in the obs or sim timeseries since the last split. This will make the iterative SD
%   algorithm very slow. Therefore, in this case, keep the candidate split point nevertheless

if(~exist('split_frequency', 'var') | isnan(split_frequency))
    warning('split_frequency IS EITHER MISSING OR NaN')
    return
end

max_quantile = 0.50;    % candidate split points are only kept if both the obs and sim value are low enough to be below this probability of unexceedance of the entire time series
max_num_segs = 15;      % a candidate split point is kept if since the last split point, more segments are contained in the obs or sim series
perc = 15;              % the region around a split point used for searching the optimal split point (percent of the split frequency)

searchrange = round((perc/100) * split_frequency);

% determine t
obs_max_global = quantile(obs,max_quantile);
sim_max_global = quantile(sim,max_quantile);

%% build the list of splits

% add the start of the time series (mandatory)
timeseries_splits(1) = 1;

for i = split_frequency : split_frequency : length(obs)-split_frequency
    obs_sample = obs(i-searchrange : i+searchrange);
    sim_sample = sim(i-searchrange : i+searchrange);
    [~,~,obs_ranks] = unique(obs_sample);
    [~,~,sim_ranks] = unique(sim_sample);
    best_split_time = find(((obs_ranks + sim_ranks) == min(obs_ranks + sim_ranks)),1);
    best_split_time = i - searchrange + best_split_time - 1;
    if (obs(best_split_time) <= obs_max_global) && (sim(best_split_time) <= sim_max_global)
        % if at the candidate split point, both obs and sim are smalle enough, add the point in time to the list
        timeseries_splits = [timeseries_splits; best_split_time];
    else
        % if there are many segments in obs or sim since the last split point, keep the split point anyways
        times = (timeseries_splits(end) : 1: best_split_time)';
        vals_obs = obs(times);
        vals_sim = sim(times);
        segs_obs = f_DefineSegments(times,vals_obs);
        segs_sim = f_DefineSegments(times,vals_sim);
        if (size(segs_obs,2) >= max_num_segs) || (size(segs_sim,2) >= max_num_segs)
            timeseries_splits = [timeseries_splits; best_split_time];
        end
    end
end

% add the end of the time series (mandatory)
timeseries_splits = [timeseries_splits; length(obs)];

end

