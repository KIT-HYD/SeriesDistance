function [segs_obs_opt, segs_sim_opt, cons, connector_data, ObFuncVal, opt_step, CoarseGrain_segs, seg_raw_statistics] = ...
           f_CoarseGraining_Event(obs, obs_eventindex, sim, sim_eventindex, weight_nfc, weight_rds, ...
           weight_sdt, weight_sdv, error_model, plot_intermedSteps)

% open issues %
    % revise/ shorten pre-processing
    % re-structure and comments in the big-loop - finde Prozedere zum Teil noch schwer nachvollziehbar.
    
% INPUT
    % obs: (n,1)        matrix with observed discharge
    % obs_eventindex    index of observed events
    % sim: (n,1)        matrix with simulated discharge
    % sim_eventindex    index of simulated events
    % weight_ndf: (1)   weighting factor for the number of false hydrological cases, used in the objective function.
    %                   The higher the more relevant. Recommended value: 1
    % weight_rds: (1)   weighting factor for the relevance of the deleted segments, used in the objective function.
    %                   The higher the more relevant. Recommended value: 1 
    % weight_sdt: (1)   weighting factor for the series distance time error, used in the objective function.
    %                   The higher the more relevant. Recommended value: 5
    % weight_sdv: (1)   weighting factor for the series distance value (or magnitude) error, used in the objective function.
    %                   The higher the more relevant. Recommended value: 0         
    % error_model: (1)  sets the way the magnitude distance among obs and sim is computed
    %                   but only in the last step, (section %%'compute and return the final, opt_stepimized series distances'
    %                   with the opt_stepimized set of segments. During opt_stepimization, as simple distance (obs - sim) is used
    %                   if 'true': dist_v = (obs - sim) / ((obs + sim)*0.5)
    %                   if 'false': dist_v = (obs - sim)
    %                   Recommended: 'true'    
    % plot_intermedSte  Plots intermediate coarse graining steps
        
%% pre-processing
    obs = obs(obs_eventindex);
    sim = sim(sim_eventindex);

    % Trim sim and obs to ensure that both start and end with either rise or fall (if necessary)
    [obs x_obs sim x_sim] = f_TrimSeries(obs, obs_eventindex, sim, sim_eventindex);
        
    % Determine the hydrological case for each timestep in the original time series
    hydcase_obs_orig = f_CalcHydCase(obs); 
    hydcase_sim_orig = f_CalcHydCase(sim);
    
    hydcase_obs = hydcase_obs_orig;     % will change with increasing segment merging
    hydcase_sim = hydcase_sim_orig;     % will change with increasing segment merging

    % define segments in the two time series
    segs_obs = f_DefineSegments (x_obs,obs);
    segs_sim = f_DefineSegments (x_sim,sim);
    
    % determine statistics of initial segment definition
    seg_raw_statistics = [f_SegStats(segs_obs) f_SegStats(segs_sim)];
    
    % check for differences in the number of segments of obs and sim
    seg_diff = length(segs_obs) - length(segs_sim); 
      
        % error checking: events must have either both even or both odd # of segments
        if mod (seg_diff,2) ~= 0
            error('f_SeriesDistance: events must have either both even or both odd # of segments!');
        end
    
    % equalize the # of segments starting with the least relevant segment in the event which has more segments
    while seg_diff ~= 0     % only required if the number of segments differs 
        if seg_diff > 0     % more obs than sim segments
           [segs_obs, hydcase_obs] = f_AggregateSegment(segs_obs,hydcase_obs,obs);  % erase the least relevant segment
        else                % more sim than obs segments
           [segs_sim, hydcase_sim] = f_AggregateSegment(segs_sim,hydcase_sim,sim);  % erase the least relevant segment
        end    
        seg_diff = length(segs_obs) - length(segs_sim);                         % number of segments still unequal?
    end

    % cleanup
    clear seg_diff x_sim x_obs

%% Reduction of segments (coarse graining) and calculation of statistics of aggreement

% determine number of reduction steps
    num_red = floor((length(segs_obs)/2))-1 ;  % reduce to 2 (if starting with an even number of segments)
                                                   % reduce to 3 (if starting with an odd number of segments)    
% initialize arrays for the objective function
    percfalsecase = NaN(num_red+1,1);       % number of wrong hydcase assignments
    mafdist_t = NaN(num_red+1,1);           % Mean Absolute Time Error of SD [h]
    mafdist_v = NaN(num_red+1,1);           % Mean Absolute Value Error of SD [m3/s]
    segment_data = cell(num_red+1,1);       % (num_red,3) cell array which contains the best obs and sim segments (col 1 and 2) found for each reduction step and joint SD proerties (col 3)
    connector_data = cell(num_red+1,1);     % cell array which contains the connectors of the different coarse graining steps
    % NOTE: first entries contain the initial state before coarse graining

% Apply SD and calculate all three statistics for initial conditions (no reduction, only equalized # of segments)
    % apply SD 
    [fdist_q, fdist_t,~,~,~,~,~,~,~,cons] = f_SD(obs,segs_obs,sim,segs_sim,error_model,'false');
   
    % plot initial conditions
    if plot_intermedSteps==true;
        f_PlotCoarseGrainIntSteps(obs,segs_obs,sim,segs_sim,cons,'initial conditions') 
    end
    
    % store segments and connectors for initial conditions
    segment_data{1,1} = segs_obs;
    segment_data{1,2} = segs_sim;
    connector_data{1,1} = cons;
    connector_data{1,2} = 0; % field stores coarse graining step; 0=inital conditions
    
  % calculate objective function inputs for initial conditions
    percfalsecase(1) = length(find(hydcase_obs_orig~=hydcase_obs))/length(obs)... % initial number of wrong hydcase assignments
        + length(find(hydcase_sim_orig~=hydcase_sim))/length(sim);
    mafdist_t(1)=mean(abs(fdist_t));     % mean SD time error
    mafdist_v(1)=mean(abs(fdist_q));     % mean SD value error    
    
    % store segment combinations for initial conditions
    CoarseGrain_segs = [];
    current_segs = [ones(size([segs_obs.starttime_global],2),1)'; segs_obs.starttime_global; segs_obs.endtime_global; segs_sim.starttime_global; segs_sim.endtime_global]';
    CoarseGrain_segs = [CoarseGrain_segs; current_segs]; % rbind
    
%% Iterative coarse-graining: Jointly aggregate segments in obs and sim, one by one, until the event is represented by two obs and two sim segments
for z = 1: num_red  % reduce until only 2 or 3 segments are left (2: when started with even # of segments, 3: when started with odd # of segments)
    % error checking
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
           
            % compute percentage of false hyd. cases and relative importance
            tmp_percfalsecase(z_obs,z_sim)= (length(find(hydcase_obs_orig~=tmp_hydcase_obs))/length(obs)) ...
                                          + (length(find(hydcase_sim_orig~=tmp_hydcase_sim))/length(sim));
            tmp_rel_del_seg(z_obs,z_sim) = tmp_rel_obs + tmp_rel_sim;
            
            % Apply SD
            [fdist_q, fdist_t,~,~,~,~,~,~,~,cons] = f_SD(obs,tmp_segs_obs,sim,tmp_segs_sim,error_model, 'true');  % UE: used to be 'standard'
            tmp_mafdist_t(z_obs,z_sim)=mean(abs(fdist_t));
            tmp_mafdist_v(z_obs,z_sim)=mean(abs(fdist_q));
        end
    end   
        
    %% find the best erase-combination for the given step uing an objective function
    % norm and weight the criteria. NOTE: For all criteria: the smaller = the better 0=best, 1=worst
    norm_tmp_percfalsecase = f_normalize(tmp_percfalsecase);
    norm_tmp_rel_del_seg = f_normalize(tmp_rel_del_seg);
    norm_tmp_mafdist_t = f_normalize(tmp_mafdist_t);
    norm_tmp_mafdist_v = f_normalize(tmp_mafdist_v);

    % join the criteria to calculate the objective function (euclidean distance)   
    tmp_opt_step = NaN(length(segs_obs));        
    for zzz =1: length(tmp_opt_step)
        for yyy = 1: length(tmp_opt_step)
            tmp_opt_step(yyy,zzz) = sqrt(weight_nfc*norm_tmp_percfalsecase(yyy,zzz)^2 + ...
                                         weight_rds*norm_tmp_rel_del_seg(yyy,zzz)^2 + ...
                                         weight_sdt*norm_tmp_mafdist_t(yyy,zzz)^2 + ...
                                         weight_sdv*norm_tmp_mafdist_v(yyy,zzz)^2);
        end
    end

    % find the minimum (=best) value
    [pos_obs,pos_sim] = find (tmp_opt_step == min(min(tmp_opt_step)));
    pos_obs = pos_obs(1);   % reduce to size 1 in case several equally small values were found
    pos_sim = pos_sim(1);   % reduce to size 1 in case several equally small values were found  

    % execute the change on the real events
    [segs_obs, hydcase_obs] = f_AggregateSegment(segs_obs,hydcase_obs,obs,pos_obs); % erase the specified segment
    [segs_sim, hydcase_sim] = f_AggregateSegment(segs_sim,hydcase_sim,sim,pos_sim); % erase the specified segment

    % calculate Series Distance
    [fdist_q, fdist_t,~,~,~,~,~,~,~,cons] = f_SD(obs,segs_obs,sim,segs_sim,error_model); % UE: used to be 'standard' 

    % store all segment data (needed for final SD calculation and plotting of the overall winner)
    segment_data{z+1,1} = segs_obs;
    segment_data{z+1,2} = segs_sim;
    connector_data{z+1,1} = cons;
    connector_data{z+1,2} = z; % add coarse graining step
    
    % compute objective function for initial conditions (before coarse graining)
    percfalsecase(z+1) = (length(find(hydcase_obs_orig~=hydcase_obs))/length(obs)) ...
                       + (length(find(hydcase_sim_orig~=hydcase_sim))/length(sim));
    mafdist_t(z+1) = mean(abs(fdist_t));
    mafdist_v(z+1) = mean(abs(fdist_q));

    % Store segment combinations and corresponding coarse graining steps
    rowindex = (z+1)*ones(size([segs_obs.starttime_global],2),1); % reduction step index
    current_segs = [rowindex'; segs_obs.starttime_global; segs_obs.endtime_global; segs_sim.starttime_global; segs_sim.endtime_global]';
    CoarseGrain_segs = [CoarseGrain_segs; current_segs]; % rbind
    
    % plot intermediate coarse graining steps
    if plot_intermedSteps==true;
        f_PlotCoarseGrainIntSteps(obs,segs_obs,sim,segs_sim,cons,strcat('coarse graining step:', {' '}, num2str(z))) 
    end
    
    % display progress info
    disp(['coarse graining step ' num2str(z) ' of ' num2str(num_red)]);

end % end of huge loop

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
    segs_obs_opt = cell2mat(segment_data(opt_step,1));
    segs_sim_opt = cell2mat(segment_data(opt_step,2));
    cons = cell2mat(connector_data(opt_step,1));

end
