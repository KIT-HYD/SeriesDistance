function [e_q_1d, cons] = f_SD_1dNoEventError(obs, sim, obs_events, sim_events, obs_sim_events_mapped,error_model)
% Determines the 1-d error distribution between obs and sim for all time steps that are neither part of an obs nor a sim event

bad_times = [];  % list of all times that shall NOT be used to determine the NoEventError, because they are part of either an obs or sim event
cons = struct('x_match_obs_global',{},'y_match_obs',{},'x_match_sim_global',{},'y_match_sim',{} );

% loop over all items in the matching events list
    for i = 1 : size(obs_sim_events_mapped,1) 

        start_obs = obs_sim_events_mapped(i,1);
        start_sim = obs_sim_events_mapped(i,2);               
        end_obs = obs_events( find(obs_events(:,1)==start_obs),2);
        end_sim = sim_events( find(sim_events(:,1)==start_sim),2);

        bad_start = min(start_obs,start_sim);   % find the earlier start of the related obs and sim event
        bad_end = max(end_obs, end_sim);        % find the later end of the related obs and sim event

        bad_times = [bad_times bad_start:bad_end];  % exclude these time steps from the 1d error calculation of low-flows

    end

% the list of times that shall be used for the NoEventError are ALL times minus the times that shall NOT be used
    good_times = 1:size(obs,1);
    good_times(bad_times) = [];

% determine the errors 
    if strcmp(error_model,'standard') % compute the simple difference 
        e_q_1d = obs(good_times) - sim(good_times);  % > 0 means obs is larger than sim
    elseif strcmp(error_model,'relative') % compute a scaled difference
        e_q_1d = (obs(good_times) - sim(good_times)) ./ ((obs(good_times) + sim(good_times)) .* 0.5); % > 0 means obs is larger than sim
    else
        error('distance function not properly specified');
    end

% save matching points(for plotting)
    cons(1).x_match_obs_global =  good_times;
    cons(1).y_match_obs = obs(good_times)';
    cons(1).x_match_sim_global = good_times;
    cons(1).y_match_sim = sim(good_times)'; 

end

