function [contingency_table] = f_ComputeContingencyTable(obs_events, sim_events, obs_sim_events_mapped )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
contingency_table = struct('hits',NaN,'misses',NaN,'false_alarms',NaN);

contingency_table.hits = size(obs_sim_events_mapped,1);
contingency_table.misses = size(obs_events,1) - contingency_table.hits;
contingency_table.false_alarms = size(sim_events,1) - contingency_table.hits;

end

