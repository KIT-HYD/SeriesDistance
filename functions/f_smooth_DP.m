function [obs sim] = f_smooth_DP(obs, sim, nse_smooth_limit)

% calcualte statistics to compare original and smoothed time series
    % number of extremes
    local_mins = length(find(diff(sign(diff(obs)))==+2)+1); % number of positive sign changes in the temporal derivative of obs (local maxima) 
    local_maxs = length(find(diff(sign(diff(obs)))==-2)+1); % number of negative sign changes in the temporal derivative of obs (local minima) 
    obs_tot_extremes = local_mins + local_maxs;
    
    local_mins = length(find(diff(sign(diff(sim)))==+2)+1); % same for sim
    local_maxs = length(find(diff(sign(diff(sim)))==-2)+1); % same for sim
    sim_tot_extremes = local_mins + local_maxs;

    SumAbsSIM = sum(abs(diff(sim)));
    SumAbsOBS = sum(abs(diff(obs)));
    
    disp(strcat('orginal obs: var: ', num2str(var(obs)), ', # extremes: ', num2str(obs_tot_extremes), ', diff(obs)=', num2str(SumAbsOBS), fprintf('\n'), ...
        '; orginal sim: var: ', num2str(var(sim)), ', # extremes: ', num2str(sim_tot_extremes), ...
        ', diff(sim)=', num2str(SumAbsSIM)))
    
    obs_old = obs;              % keep the original time series
    xes = (1:1:size(obs,1))';   % create x data (points in time)
    xy = cat(2, xes, obs_old);  % prepare input for f_dp1d
    % added: % simplify obs up to a point number criterion
    xy_dp = f_dp1d(xy,-999,sim_tot_extremes);  
        % former version:
        % simplify obs up to a nse-agreement level specified by 'nse_limit_obs'
        %xy_dp = f_dp1d(xy,-999,-999,nse_smooth_limit);  
    % sample the simplified line at the original x-locations (points in time)  
    obs = interp1(xy_dp(:,1),xy_dp(:,2),xes,'linear'); 
    
%     % repeat for sim 
%     sim_old = sim;              % keep the original time series
%     xes = (1:1:size(sim,1))';   % create x data (points in time)
%     xy = cat(2, xes, sim_old);  % prepare input for f_dp1d
%     % simplify 'sim' such that it has the same number of landmarks (or peaks and troughs) as 'obs'
%     num_points = size(xy_dp,1); % the number of points of the simplified 'sim' time series   
%     xy_dp = f_dp1d(xy,-999,num_points);  
%     % sample the simplified line at the original x-locations (points in time)  
%     sim = interp1(xy_dp(:,1),xy_dp(:,2),xes,'linear');  
 
     % number of extremes smoothed time series
    local_mins = length(find(diff(sign(diff(obs)))==+2)+1); % number of positive sign changes (local maxima) 
    local_maxs = length(find(diff(sign(diff(obs)))==-2)+1); % number of negative sign changes (local minima) 
    obs_tot_extremes = local_mins + local_maxs;
    
    local_mins = length(find(diff(sign(diff(sim)))==+2)+1); % number of positive sign changes (local maxima) 
    local_maxs = length(find(diff(sign(diff(sim)))==-2)+1); % number of negative sign changes (local minima) 
    sim_tot_extremes = local_mins + local_maxs;
    
    SumAbsSIM = sum(abs(diff(sim)));
    SumAbsOBS = sum(abs(diff(obs)));

    % plot statistics on the command window
    disp(strcat('smoothed obs: var: ', num2str(var(obs)), ', # extremes: ', num2str(obs_tot_extremes), ', diff(obs)=', num2str(SumAbsOBS), fprintf('\n'), ...
        '; smoothed sim: var: ', num2str(var(sim)), ', # extremes: ', num2str(sim_tot_extremes), ...
        ', diff(sim)=', num2str(SumAbsSIM)))
   
end

