function [obs, sim, obs_events, sim_events, obs_sim_events_mapped] = f_inputsubset(to, obs, sim, obs_events, sim_events, obs_sim_events_mapped)
    obs = obs(1:to);
    sim = sim(1:to);
    if isempty(obs_events) 
        sim_events = [];
        obs_events = [];
        obs_sim_events_mapped = [];
    else
        index = find(obs_events(:,1)>1 & obs_events(:,2)<to);
        sim_events = sim_events(index,:);
        obs_events = obs_events(index,:);
        obs_sim_events_mapped = obs_sim_events_mapped(index,:);
    end
end
