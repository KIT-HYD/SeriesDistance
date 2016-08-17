function [obs x_obs sim x_sim] = f_TrimSeries(obs, obs_eventindex, sim, sim_eventindex)
% Trims if necessary the obs and sim series until
    % they both start with either a 'rise' or a 'fall'
    % they both end with either a 'rise' or a 'fall'

% INPUT
    % obs: (n,1) matrix with observed discharge. n = number of time steps
    % hydcase_obs: (n,1) matrix with observed hydrological case: -2=valley -1=drop, 1=rise 2=peak      
    % sim: (n,1) matrix with simulated discharge. n = number of time steps
    % hydcase_sim: (n,1) matrix with simulated hydrological case: -2=valley -1=drop, 1=rise 2=peak 
    
% OUTPUT
    % start_obs: new start of the observed series
    % end_obs: new end of the observed series   
    % start_sim: new start of the simulated series
    % end_sim: new end of the simulated series   
    
% Determine the hydrological case for each timestep in the original time series
    hydcase_obs = f_CalcHydCase(obs);  % -2=valley -1=drop 0=no feature 1=rise 2=peak
    hydcase_sim = f_CalcHydCase(sim);

% find the best starting points
    pos_obs_rise = find(hydcase_obs == 1,1,'first');
    pos_sim_rise = find(hydcase_sim == 1,1,'first'); 
    sum_pos_rise = pos_obs_rise + pos_sim_rise;

    pos_obs_fall = find(hydcase_obs == -1,1,'first');
    pos_sim_fall = find(hydcase_sim == -1,1,'first'); 
    sum_pos_fall = pos_obs_fall + pos_sim_fall;

    % choose the starting point pair where the least trimming is required
    if sum_pos_rise <= sum_pos_fall
        start_obs = pos_obs_rise;
        start_sim = pos_sim_rise;
    else
        start_obs = pos_obs_fall;
        start_sim = pos_sim_fall; 
    end
  
% find the best end points
    pos_obs_rise = find(hydcase_obs == 1,1,'last');
    pos_sim_rise = find(hydcase_sim == 1,1,'last'); 
    sum_pos_rise = pos_obs_rise + pos_sim_rise;

    pos_obs_fall = find(hydcase_obs == -1,1,'last');
    pos_sim_fall = find(hydcase_sim == -1,1,'last'); 
    sum_pos_fall = pos_obs_fall + pos_sim_fall;
    
    % choose the end point pair where the least trimming is required
    if sum_pos_rise > sum_pos_fall
        end_obs = pos_obs_rise;
        end_sim = pos_sim_rise;
    else
        end_obs = pos_obs_fall;
        end_sim = pos_sim_fall; 
    end
    
% memorize temporal offsets of trimmed time series
    offset_obs = obs_eventindex(1) + start_obs-1; % offset of trimmed obs series
    offset_sim = sim_eventindex(1) + start_sim-1; % offset of trimmed sim series

% trim the time series
    obs = obs(start_obs:end_obs);                           % obs trimmed
    x_obs = (offset_obs:(offset_obs + size(obs,1) - 1))';      % the global time position of obs
    %hydcase_obs_orig = hydcase_obs_orig(start_obs:end_obs); % hydrological cases of obs

    sim = sim(start_sim:end_sim);                           % sim trimmed
    x_sim = (offset_sim:(offset_sim + size(sim,1) - 1))';      % the global time position of sim
    %hydcase_sim_orig = hydcase_sim_orig(start_sim:end_sim); % hydrological cases of sim

% update flow condition index/ hydrological cases
    hydcase_obs = f_CalcHydCase(obs); 
    hydcase_sim = f_CalcHydCase(sim);

end