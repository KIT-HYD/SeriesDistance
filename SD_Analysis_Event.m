%  ====================================================================== 
%                 Series Distance Analysis (Event Mode) 
%  ======================================================================
%  last modifaction 16.08.2016

% This programm (development release) determines SeriesDistance (SD) error distributions 
% between 'events' in observed (obs) and simulated (sim) discharge time series. Using 
% sim and obs, indices of event start and end points and information on valid pairs 
% of events, the program determines seperate error distributions for: 
%  - rising limbs within events (2d error distribution of time- and magnitude errors)
%  - falling limbs within events (2d error distribution of time- and magnitude errors)
%  - no-event times (low-flow, 1d error distribution of magnitude errors)

% Further information:
% Seibert, S. P., Ehret, U., and Zehe, E., 2016: Disentangling timing and amplitude errors in streamflow simulations, Hydrol. Earth Syst. Sci.
% Ehret, U., Zehe, E., 2011. Series distance - An intuitive metric to quantify hydrograph similarity in terms of occurrence, amplitude and timing of hydrological events. Hydrol. Earth Syst. Sci. 15, 877–896. doi:10.5194/hess-15-877-2011

% Dependencies
% All required functions are stored in folder \functions

% Input (ascii format)
%  - obs: [n,1] matrix with equidistant and NaN-free time series of observed discharge data.
%  - sim: [n,1] matrix with equidistant and NaN-free time series of simulated discharge data. 
%  - obs_events: [m,2] matrix with start and end times of events in 'obs'. [m,1]=start, [m,2]=end.
%  - sim_events: [m,2] matrix with start and end times of events in 'sim'. [m,1]=start, [m,2]=end. 
%  - obs_sim_pairing: [m,2] matrix with events in 'obs_events' and 'sim_events' that belong together.  

% Parameters are explained and specified in the parameter block

% Output (all stored within a single binary file .mat)
%  - obs:               observed discharge (smoothed)
%  - sim:               simulated discharge
%  - segs_obs_opt_all:  struct that contains the coarse-grained segments of obs
%  - segs_sim_opt_all:  struct that contains all coarse-grained segments of sim
%  - connectors:        struct with SD connectors for periods of low flow and events
%  - e_sd_q_rise:       SD magnitude errors for rising limbs
%  - e_sd_t_rise:       SD timing errors for rising limbs
%  - e_sd_q_fall:       SD magnitude errors for falling limbs
%  - e_sd_t_fall:       SD timing errors for falling limbs
%  - e_sd_lowFlow:      1D errors for low-flow
%  - error_model:       'relative' or 'standard'
%  - seg_raw_statistics: statistics of the original time series e.g. number of peaks, throughs, duration of rising/ fallings segments
%  - seg_opt_statistics: statistics of the coarse-grained time series e.g. number of peaks, throughs, duration of rising/ fallings segments

%% Start from scratch
    clc;
    clear all;
    close all;

%% specify paths and parameters, read inputs
    cd('ADD_YOUR_PATH_HERE')  % set working directory
    % cd('d:/02_Projekte/03_in_Arbeit_intern/1012_SeriesDistance/01_versions/SD_GitHub_public_v1.0/')  % set working directory
    addpath(genpath(pwd)); % add working directory to searchpath

% read input
    obs = dlmread('data/HOST_timeseries.csv', ';', [1 2 10000 2]);
    sim = dlmread('data/HOST_timeseries.csv', ';', [1 3 10000 3]);
    obs_events = dlmread('data/HOST_obs_events.csv', ';');
    sim_events = dlmread('data/HOST_sim_events.csv', ';');
    obs_sim_pairing = dlmread('data/HOST_event_pairing.csv',';'); % 1 row = 1 event pair. [n,1]=start time of obs event, [p,2]= start time of sim event. (optional)

% output filename
    outfile = './results/output.mat';
   
% smoothing options
    smooth_flag = true;         % smooth both obs and sim (default=true)
    nse_smooth_limit = 0.99;    % specifies degree of smoothing according to NSE criterion (default=0.99)

% specification of the magnitude error model
    error_model = 'relative'; % 'relative' or 'standard'; (standard: (sim-obs), relative: [(sim-obs)/((sim + obs)/2)] (default='relative')

% parametrization of the objective function 
    weight_nfc = 1/7;   % weights number of re-assigned hydrological cases (default= 1/7)    
    weight_rds = 1/7;   % weights the importance of the re-assigned segments (default=1/7) 
    weight_sdt = 5/7;   % weights the SD timing error component (default=5/7)
    weight_sdv = 0;     % weights the SD magnitude error component (default=0)
   
% set plot flags 
    pf_input = true;                    % plots smoothed and original input time series ('obs' and 'sim')
    pf_CoarseGrainSteps = false;        % plots intermediate coarse graining steps and progression of objective function (NOTE: 'true' leads to MANY plots)
    pf_segs_cons_indivEvents = false;   % plots each individual event with optimized segments and connectors in an own figure
    pf_objective_functions = false;     % plots the SD objective function for each event in an own figure
    pf_segs_cons_entireSeries = true;   % plots the entire time series with optimized segments and connectors in an individual figure
    pf_errorDistributions = true;       % SeriesDistance error distributions for rise, fall and low-flow 

%% Data manipulations and pre-processing
% allows to select a (1:to) subset of the provided time series, the remainder is ignored.
%     to = 2000;
%     [obs, sim, obs_events, sim_events, obs_sim_pairing] = f_inputsubset(to, obs, sim, obs_events, sim_events, obs_sim_pairing);
%     clear to

% smooth if required (default=true)
    if smooth_flag == true
        obs_org = obs;
        sim_org = sim;
        [obs, sim] = f_smooth_DP(obs, sim, nse_smooth_limit);    
    end
    
% replace identical neigboring values to avoid problems with assignment of unique peaks and valleys
    obs = f_ReplaceEqualNeighbours(obs);
    sim = f_ReplaceEqualNeighbours(sim);

% plot input data
    if pf_input == true
        f_PlotInput(obs_org, obs, obs_events, sim_org, sim, sim_events, obs_sim_pairing, []);
    end
    
% cleanup
    clear smooth_flag nse_smooth_limit pf_input events_mapped_by_user pf_input
    
%% initialize output arrays for both, coarse graining and SD 
    e_sd_t_rise = [];       % error distribution for events, rise, time component
    e_sd_q_rise = [];       % error distribution for events, rise, magnitude component
    e_sd_t_fall = [];       % error distribution for events, fall, time component
    e_sd_q_fall = [];       % error distribution for events, fall, magnitude component
    segs_obs_opt_all = [];  % coarse-grained segments of 'obs'
    segs_sim_opt_all = [];  % coarse-grained segments of 'sim'
    seg_raw_statistics = [];% segment statistics
    seg_opt_statistics = [];% segment statistics
    seg_event_id = [];      % stores event ID for each segs_opt object
    connectors = struct('x_match_obs_global',{},'y_match_obs',{},'x_match_sim_global',{},'y_match_sim',{} );   % connectors between matching points in 'obs' and 'sim'
    
 %% Apply coarse-graining and SD method to each event 
    for ii = 1:size(obs_sim_pairing,1) % loop over all items in the matching events list            
        % print progress in command window
            disp(['event ' num2str(ii) ' of ' num2str(size(obs_sim_pairing,1))]); 
        
        % get event start and corresponding end points
            obs_eventindex = obs_sim_pairing(ii,1) : obs_events(find(obs_events(:,1)==obs_sim_pairing(ii,1)),2); % period of observed event
            sim_eventindex = obs_sim_pairing(ii,2) : sim_events(find(sim_events(:,1)==obs_sim_pairing(ii,2)),2); % period of simulated event
               
        % apply coarse graining: determines optimal level of aggregation of each event (key output: segs_xxx_opt)
            [segs_obs_opt, segs_sim_opt, cons, connector_data, ObFuncVal, opt_step, CoarseGrain_segs, seg_raw_stats] = ...
                f_CoarseGraining_Event(obs,obs_eventindex,sim,sim_eventindex,weight_nfc,weight_rds,weight_sdt,weight_sdv,error_model,pf_CoarseGrainSteps);
        
        % summarize segment statistics (before and after coarse graining)
            seg_raw_statistics = [seg_raw_statistics; seg_raw_stats]; % cols: obs_peaks, obs_troughs, obs_rise_dur, obs_fall_dur, obs_tot_dur, then analogously for sim
            seg_opt_statistics = [seg_opt_statistics; [ii f_SegStats(segs_obs_opt) ... % cols: event #, obs_peaks, obs_troughs, obs_rise_dur, obs_fall_dur, obs_tot_dur,
                (length(ObFuncVal)>1) opt_step f_SegStats(segs_sim_opt)]];             %       coarse grained (0/1), coarse grained step, sim_troughs, sim_rise_dur, sim_fall_dur, sim_tot_dur
        
        % store optimized segment combination for all events 
            segs_obs_opt_all = [segs_obs_opt_all segs_obs_opt];
            segs_sim_opt_all = [segs_sim_opt_all segs_sim_opt];
        
        % memorize corresponding event ID (only required for plotting)
            seg_event_id = [seg_event_id; ii*ones(length(segs_obs_opt),1)];
            
        % plot each individual event with optimized segments and connectors in an own figure
            if pf_segs_cons_indivEvents == true
               f_PlotConnectedSeries(obs,segs_obs_opt,sim,segs_sim_opt,cons) 
            end

        % plot the objective function values
            if pf_objective_functions == true 
                f_plot_ObjectiveFunction_CoarsGrainStps(ObFuncVal, opt_step, strcat(['event # ' num2str(ii)])) 
            end
        
        % cleanup
            clear cons seg_raw_stats ObFuncVal opt_step obs_eventindex sim_eventindex       
        
        %% compute SD results for the optimized level of generalization
        obs_fromto = segs_obs_opt(1).starttime_global:segs_obs_opt(length(segs_obs_opt)).endtime_global; % get trimmed event index
        sim_fromto = segs_sim_opt(1).starttime_global:segs_sim_opt(length(segs_sim_opt)).endtime_global; % get trimmed event index
        [~,~,~,e_q_rise, e_t_rise,~,e_q_fall, e_t_fall,~,cons, ~, ~] = f_SD(obs(obs_fromto),segs_obs_opt,sim(sim_fromto),segs_sim_opt,error_model,'true'); 
    
        % build overall SD error distributions by successively adding errors from individual events
            e_sd_t_rise = [e_sd_t_rise; e_t_rise]; % for timing errors of rising segments
            e_sd_q_rise = [e_sd_q_rise; e_q_rise]; % for magnitude errors of rising segments
            e_sd_t_fall = [e_sd_t_fall; e_t_fall]; % for timing errors of falling segments
            e_sd_q_fall = [e_sd_q_fall; e_q_fall]; % for magnitude errors of falling segments                

        % store individual connectors of each single event
            if isempty(connectors)
                connectors(1).x_match_obs_global = cons(1).x_match_obs_global;
                connectors(1).y_match_obs = cons(1).y_match_obs;
                connectors(1).x_match_sim_global = cons(1).x_match_sim_global;
                connectors(1).y_match_sim = cons(1).y_match_sim; 
            else
                connectors(1).x_match_obs_global = [connectors(1).x_match_obs_global cons(1).x_match_obs_global];
                connectors(1).y_match_obs = [connectors(1).y_match_obs cons(1).y_match_obs];
                connectors(1).x_match_sim_global = [connectors(1).x_match_sim_global cons(1).x_match_sim_global];
                connectors(1).y_match_sim = [connectors(1).y_match_sim  cons(1).y_match_sim]; 
            end
            
    end % here ends the big loop over all 'ii' items in the matching event list
    
% attach event number to segment objects (for plotting, helps identifying events)
    for kk = 1:length(segs_obs_opt_all)
        segs_obs_opt_all(kk).eventID = seg_event_id(kk);
        segs_sim_opt_all(kk).eventID = seg_event_id(kk);
    end

% cleanup    
    clear ii cons segs_obs_opt segs_sim_opt weight_nfc weight_rds weight_sdt weight_sdv obs_fromto ...
        sim_fromto pf_segs_cons_indivEvents e_t_rise e_q_rise e_t_fall e_q_fall ...
        pf_CoarseGrainSteps seg_event_id kk
 
%% determine the SeriesDistance distribution for no-event periods
[e_sd_lowFlow, cons1D] = f_SD_1dNoEventError (obs, sim, obs_events, sim_events, obs_sim_pairing, error_model);

% add the connectors for no-event cases to those of the events
    connectors(1).x_match_obs_global = [connectors(1).x_match_obs_global cons1D(1).x_match_obs_global];
    connectors(1).y_match_obs = [connectors(1).y_match_obs cons1D(1).y_match_obs];
    connectors(1).x_match_sim_global = [connectors(1).x_match_sim_global cons1D(1).x_match_sim_global];
    connectors(1).y_match_sim = [connectors(1).y_match_sim  cons1D(1).y_match_sim]; 

% determine the contingency table        
    contingency_table = f_ComputeContingencyTable(obs_events, sim_events, obs_sim_pairing);

% plot the entire time series with all coarse-grained segments and connectors
    if pf_segs_cons_entireSeries == true
        f_PlotConnectedSeries (obs,segs_obs_opt_all,sim,segs_sim_opt_all,connectors) 
    end

% cleanup
    clear pf_segs_cons_entireSeries cons1D

%% plot error distributions
    if pf_errorDistributions == true       
        % plot 2d error distributions
        f_PlotSDErrors(e_sd_q_rise, e_sd_t_rise, e_sd_q_fall, e_sd_t_fall,[],[],[],[]); 
        % plot 1d error distribution
        f_Plot1dErrors(e_sd_lowFlow,[],cellstr('SD low-flow'));            
    end        
    
% cleanup
    clear pf_errorDistributions 
    
%% save all inputs, outputs and important parameters into a single .mat file
    save(outfile, 'obs', 'obs_org', 'obs_events', 'sim', 'sim_org', 'sim_events', 'obs_sim_pairing', ...
        'segs_obs_opt_all', 'segs_sim_opt_all', 'seg_raw_statistics', 'seg_opt_statistics', 'connectors', ...
        'e_sd_q_rise', 'e_sd_t_rise', 'e_sd_q_fall', 'e_sd_t_fall', 'e_sd_lowFlow',...
        'error_model', 'contingency_table');   
         
        

        
        
    
    
    





