function [segs_obs_opt_all, segs_sim_opt_all, cons_all, e_sd_t_all, e_sd_q_all] = ...
    f_CoarseGraining_Continuous(obs, sim, timeseries_splits,weight_nfc,weight_rds,weight_sdt,weight_sdv,error_model)
        
%% initialize arrays
    cons_all = struct('x_match_obs_global',{},'y_match_obs',{},'x_match_sim_global',{},'y_match_sim',{} );   % connectors between matching points in 'obs' and 'sim'
    segs_obs_opt_all = [];  % coarse-grained segments of entire 'obs' time series
    segs_sim_opt_all = [];  % coarse-grained segments of entire 'sim' time series
    e_sd_rise_all = [];     % SD errors rising segments (final output)
    e_sd_fall_all = [];     % SD errors falling segments (final output)
    e_sd_t_all = [];        % SD timing errors (final output)
    e_sd_q_all = [];        % SD magnitude errors (final output)
    obs_org = obs;          % backup original time series
    sim_org = sim;          % backup original time series

%% pre-processing 
for i = 1 : size(timeseries_splits,1)-1 % loop over all items in timeseries_splits list
    % display progress information
    txt = ['time series split ' num2str(i) ' of ' num2str(size(timeseries_splits,1)-1)];
    disp (txt)
    
    % create subset/ split the time series
    obs_split = obs_org(timeseries_splits(i):timeseries_splits(i+1));       % extract the obs values within the time series subset
    sim_split = sim_org(timeseries_splits(i):timeseries_splits(i+1));      % extract the obs values within the time series subset

    % Trim sim and obs to ensure that both start and end with either rise or fall (if necessary)
    [obs, x_obs, sim, x_sim] = f_TrimSeries(obs_split, (timeseries_splits(i):timeseries_splits(i+1)), sim_split, (timeseries_splits(i):timeseries_splits(i+1)));
        
    % Determine the hydrological case for each timestep in the original time series
    hydcase_obs_orig = f_CalcHydCase(obs); 
    hydcase_sim_orig = f_CalcHydCase(sim);
    
    hydcase_obs = hydcase_obs_orig;     % will change with increasing segment merging
    hydcase_sim = hydcase_sim_orig;     % will change with increasing segment merging

    % define segments in the two time series
    segs_obs = f_DefineSegments (x_obs,obs);
    segs_sim = f_DefineSegments (x_sim,sim);
     
    % check for differences in the number of segments in obs and sim
    seg_diff = length(segs_obs) - length(segs_sim);
    
    % error checking: events must have either both even or both odd # of segments
    if mod (seg_diff,2) ~= 0
        error('f_SeriesDistance: events must have either both even or both odd # of segments!');
    end
    
    % equalize the # of segments starting with the least relevant segment in the time series which has more segments
    while seg_diff ~= 0
        if seg_diff > 0     % more obs than sim segments
           [segs_obs, hydcase_obs] = f_AggregateSegment(segs_obs,hydcase_obs,obs);      % erase the least relevant segment
        else                % more sim than obs segments
           [segs_sim, hydcase_sim] = f_AggregateSegment(segs_sim,hydcase_sim,sim);       % erase the least relevant segment
        end    
        seg_diff = length(segs_obs) - length(segs_sim);                             % number of segments still unequal?
    end

% cleanup
clear seg_diff

%% iterative reduction of segments and calculation of the selected statistics of aggreement

% determine number of reduction steps
    num_red = floor((length(segs_obs)/2))-1 ;  % reduce to 2 (if starting with an even number of segments)
                                               % reduce to 3 (if starting with an odd number of segments)
% initialize arrays
    percfalsecase = NaN(num_red+1,1);   % number of wrong hydcase assignments
    mafdist_t = NaN(num_red+1,1);       % Timing Error of SD
    mafdist_v = NaN(num_red+1,1);       % Magnitude Error of SD 
    segment_data = cell(num_red+1,1);   % (num_red,3) cell array which contains the best obs and sim segments (col 1 and 2) found for each reduction step and joint SD proerties (col 3)
    connector_data = cell(num_red+1,1); % cell array which contains the connectors of the different coarse graining steps
    e_sd_rise = cell(num_red+1,2);      % timing & magnitude errors of each coarse-graining step (rising)
    e_sd_fall = cell(num_red+1,2);      % timing & magnitude errors of each coarse-graining step (falling)
    % NOTE: first entries contain the initial state before coarse graining

% Apply SD and calculate all three statistics for the initial conditions (no reduction, only equalized # of segments)
    % apply SD 
    [fdist_q, fdist_t, ~, e_q_rise, e_t_rise, ~, e_q_fall, e_t_fall, ~, cons] = ...
            f_SD(obs,segs_obs,sim,segs_sim,error_model);
   
    % store segments and connectors for initial conditions
    segment_data{1,1} = segs_obs;
    segment_data{1,2} = segs_sim;
    connector_data{1,1} = cons;
    connector_data{1,2} = 0; % field stores coarse graining step; 0=inital conditions
    e_sd_rise{1,1} = e_t_rise; % for timing errors of rising segments
    e_sd_rise{1,2} = e_q_rise; % for magnitude errors of rising segments
    e_sd_fall{1,1} = e_t_fall; % for timing errors of falling segments
    e_sd_fall{1,2} = e_q_fall; % for magnitude errors of falling segments
        
    % calculate objective function inputs for initial conditions
    percfalsecase(1) = length(find(hydcase_obs_orig~=hydcase_obs))/length(obs) ...
        + length(find(hydcase_sim_orig~=hydcase_sim))/length(sim);
    mafdist_t(1)=mean(abs(fdist_t));
    mafdist_v(1)=mean(abs(fdist_q)); 

    % apply coarse-graining to all time series splits (big for-loop): Jointly reduce obs/sim segments, one by one, until only one obs and one sim segment are left
    for z = 1: num_red  % reduce until only 2 or 3 segments are left (2: when started with even # of segments, 3: when started with odd # of segments)
        % some error checking
        if hydcase_obs(1)~= hydcase_sim(1) || length(segs_obs) ~= length(segs_sim)
            error('error in big loop');
        end

        % initialize temporary variables (temporary for one reduction step)
        tmp_percfalsecase = NaN(length(segs_obs));  % (m,m) matrix with percentage of false hydcases (obs + sim) for all possible reduction combinations
        tmp_rel_del_seg = NaN(length(segs_obs));    % (m,m) matrix with relevance of deleted segments (obs + sim)
        tmp_mafdist_t = NaN(length(segs_obs));      % (m,m) matrix with timing error
        tmp_mafdist_v = NaN(length(segs_obs));      % (m,m) matrix with value error

        % loop over all possible segment reduction combinations
        for z_obs = 2:length(segs_obs)-1  % loop over all observed segments, except the first and last

            % for each new loop, start with the best selection of the previous z-round
            tmp_segs_obs = segs_obs;        
            tmp_hydcase_obs = hydcase_obs;

            tmp_rel_obs = tmp_segs_obs(z_obs).relevance;    % save the relevance of the segment before it is deleted
            [tmp_segs_obs, tmp_hydcase_obs] = f_AggregateSegment(tmp_segs_obs,tmp_hydcase_obs,obs,z_obs); % erase the specified segment

            for z_sim = 2:length(segs_sim)-1 % loop over all simulated segments, except the first and last
                % for each new loop, start with the best selection of the previous z-round
                tmp_segs_sim = segs_sim;       
                tmp_hydcase_sim = hydcase_sim;

                tmp_rel_sim = tmp_segs_sim(z_sim).relevance; % save the relevance of the segment before it is deleted
                [tmp_segs_sim, tmp_hydcase_sim] = f_AggregateSegment(tmp_segs_sim,tmp_hydcase_sim,sim,z_sim); % erase the specified segment

                % compute the objective function
                tmp_percfalsecase(z_obs,z_sim)= (length(find(hydcase_obs_orig~=tmp_hydcase_obs))/length(obs)) ...
                                              + (length(find(hydcase_sim_orig~=tmp_hydcase_sim))/length(sim));
                tmp_rel_del_seg(z_obs,z_sim) = tmp_rel_obs + tmp_rel_sim;

                [fdist_q, fdist_t,~,~,~,~,~,~,~,cons,~,~] = f_SD(obs,tmp_segs_obs,sim,tmp_segs_sim,error_model, 'true');  % UE: used to be 'standard'
                tmp_mafdist_t(z_obs,z_sim)=mean(abs(fdist_t));
                tmp_mafdist_v(z_obs,z_sim)=mean(abs(fdist_q));
            end
        end   

        %% find the best erase-combination for the given reduction step
        % norm and weight the criteria. NOTE: For all criteria: the smaller = the better 0=best, 1=worst
        norm_tmp_percfalsecase = weight_nfc * f_normalize(tmp_percfalsecase);
        norm_tmp_rel_del_seg = weight_rds * f_normalize(tmp_rel_del_seg);
        norm_tmp_mafdist_t = weight_sdt * f_normalize(tmp_mafdist_t);
        norm_tmp_mafdist_v = weight_sdv * f_normalize(tmp_mafdist_v);

        % join the criteria to calculate the objective function (euclidean distance)   
        tmp_opt = NaN(length(segs_obs));        
        for zzz =1: length(tmp_opt)
            for yyy = 1: length(tmp_opt)
                tmp_opt(yyy,zzz) = sqrt(norm_tmp_percfalsecase(yyy,zzz)^2 + norm_tmp_rel_del_seg(yyy,zzz)^2 + norm_tmp_mafdist_t(yyy,zzz)^2 + norm_tmp_mafdist_v(yyy,zzz)^2 );
            end
        end

        % find the minimum (=best) value
        [pos_obs,pos_sim] = find (tmp_opt == min(min(tmp_opt)));
        pos_obs = pos_obs(1);   % reduce to size 1 in case several equally small values were found
        pos_sim = pos_sim(1);   % reduce to size 1 in case several equally small values were found  

        % execute the change on the real events
        [segs_obs, hydcase_obs] = f_AggregateSegment(segs_obs,hydcase_obs,obs,pos_obs); % erase the specified segment
        [segs_sim, hydcase_sim] = f_AggregateSegment(segs_sim,hydcase_sim,sim,pos_sim); % erase the specified segment

        %% Calculate Series Distance on the optimized level of aggregated segments return SD errors
        [fdist_q, fdist_t, ~, e_q_rise, e_t_rise, ~, e_q_fall, e_t_fall, ~, cons] = ...
            f_SD(obs,segs_obs,sim,segs_sim,error_model); 
        
        % add segment data, connectors and time/ magnitude errors of the
        % best solution for this time series split to that of the entire time series
        segment_data{z+1,1} = segs_obs;
        segment_data{z+1,2} = segs_sim;
        connector_data{z+1,1} = cons;
        connector_data{z+1,2} = z; % add coarse graining step
        e_sd_rise{z+1,1} = e_t_rise; % for timing errors of rising segments
        e_sd_rise{z+1,2} = e_q_rise; % for magnitude errors of rising segments
        e_sd_fall{z+1,1} = e_t_fall; % for timing errors of falling segments
        e_sd_fall{z+1,2} = e_q_fall; % for magnitude errors of falling segments 
    
        % compute objective function inputs
        percfalsecase(z+1) = (length(find(hydcase_obs_orig~=hydcase_obs))/length(obs)) ...
                           + (length(find(hydcase_sim_orig~=hydcase_sim))/length(sim));
        mafdist_t(z+1) = mean(abs(fdist_t));
        mafdist_v(z+1) = mean(abs(fdist_q));

        % progress info
        txt = ['reduction step ' num2str(z) ' of ' num2str(num_red)];
        disp (txt)
    end

%% Calculate objective function and find the optimal coarse graining step
    ObFuncVal = sqrt(weight_nfc*f_normalize(percfalsecase).^2 + ...
                     weight_sdt*f_normalize(mafdist_t).^2 + ...
                     weight_sdv*f_normalize(mafdist_v).^2);

    opt_step = find(ObFuncVal == min(ObFuncVal));
    opt_step = opt_step(1); % take the first if multiple identical values exist
    if length(ObFuncVal)>1  % display best coarse graining step
        if opt_step==1 
            disp('selected step # intial conditions');     
        else        
            disp(['selected step #' num2str(opt_step)]);     
        end
    end

%% select and return coarse-grained segments and connectors for optimal level of generalization
    segs_obs_opt = cell2mat(segment_data(opt_step,1));  % segment data obs
    segs_sim_opt = cell2mat(segment_data(opt_step,2));  % segment data sim
    cons = cell2mat(connector_data(opt_step,1));        % connectors
    e_sd_rise_opt = cell2mat(e_sd_rise(opt_step,:));        % timing & magnitude errors of rising segments
    e_sd_fall_opt = cell2mat(e_sd_fall(opt_step,:));        % timing & magnitude errors of falling segments
    
%% add segment data and connectors of the splitted subset to that of the entire time series
    if isempty(cons_all)
        cons_all(1).x_match_obs_global =  cons(1).x_match_obs_global;
        cons_all(1).y_match_obs = cons(1).y_match_obs;
        cons_all(1).x_match_sim_global = cons(1).x_match_sim_global;
        cons_all(1).y_match_sim = cons(1).y_match_sim; 
    else
        cons_all(1).x_match_obs_global = [cons_all(1).x_match_obs_global cons(1).x_match_obs_global];
        cons_all(1).y_match_obs = [cons_all(1).y_match_obs cons(1).y_match_obs];
        cons_all(1).x_match_sim_global = [cons_all(1).x_match_sim_global cons(1).x_match_sim_global];
        cons_all(1).y_match_sim = [cons_all(1).y_match_sim  cons(1).y_match_sim]; 
    end
    % segment data
    segs_obs_opt_all = [segs_obs_opt_all segs_obs_opt];
    segs_sim_opt_all = [segs_sim_opt_all segs_sim_opt]; 
    % time-magnitude errors for all rising/ falling sections
    e_sd_rise_all = [e_sd_rise_all; e_sd_rise_opt];      % timing & magnitude errors of entire time series
    e_sd_fall_all = [e_sd_fall_all; e_sd_fall_opt];      % timing & magnitude errors of entire time series
    % time and magnitude errors for entire hydrographs
    e_sd_t_all = [e_sd_t_all; [e_sd_rise_opt(:,1); e_sd_fall_opt(:,1)]];      % timing & magnitude errors of entire time series
    e_sd_q_all = [e_sd_q_all; [e_sd_rise_opt(:,2); e_sd_fall_opt(:,2)]];      % timing & magnitude errors of entire time series
end

end

