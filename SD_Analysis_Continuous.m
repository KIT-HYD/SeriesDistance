%  ====================================================================== 
%                 Series Distance Analysis (Continuous Mode) 
%  ======================================================================
%  last modifaction 16.08.2016

% This development release of the SeriesDistance (SD) method is applicable to continuous 
% observed (obs) and simulated (sim) discharge time series. It outputs a 2d-error 
% distribution for the entire time series. Contrary to the event based method the 
% continuous version does not differentiate periods of low-flow from periods of events.  
% Nevertheless, the entire hydrograph is classified into rising and falling limbs. The 
% interpretation of these classes is however only meaningful if the periods of low-flow 
% do not cover a significant portion of the hydrograph. Note: For plotting,
% the error distributions for rise and fall are combined.

% Further information:
% Seibert, S. P., Ehret, U., and Zehe, E., 2016: Disentangling timing and amplitude errors in streamflow simulations, Hydrol. Earth Syst. Sci.
% Ehret, U., Zehe, E., 2011. Series distance - An intuitive metric to quantify hydrograph similarity in terms of occurrence, amplitude and timing of hydrological events. Hydrol. Earth Syst. Sci. 15, 877–896. doi:10.5194/hess-15-877-2011

% Dependencies
% All required functions are stored in folder \functions

% Input (ascii format)
%  - obs: [n,1] matrix with equidistant and NaN-free time series of observed discharge data.
%  - sim: [n,1] matrix with equidistant and NaN-free time series of simulated discharge data. 
%  - timeseries_splits: [n,1] optional matrix with points in time where the 'obs' and 'sim' time series will be split to increase computational speed.
%     It must include the first and the last timestep of 'obs' and 'sim'. E.g. if 'obs' and 'sim' are [1100,1],
%     then 'timeseries_splits' could be (1, 250, 800, 1100)

% Parameters are explained and specified in the parameter block

% Outputs (all stored within a single binary file .mat)
%  - obs:              observed discharge (smoothed)
%  - sim:              simulated discharge
%  - segs_obs_opt_all: optimized segments in obs
%  - segs_sim_opt_all: optimized segments in sim
%  - connectors:       SD connectors
%  - e_sd_q_all:       SD magnitude errors for entire time series
%  - e_sd_t_all:       SD timing errors for entire time series
%  - parameters:       all parameters, i.e. error_model, objective function weights, smoothing parameters,  
%                       are included in the outputfile

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
 
% output filenames
    outfile = './results/output.mat';

% smoothing options
    smooth_flag = true;         % smooth both obs and sim (default=true)
    nse_smooth_limit = 0.99;    % specifies degree of smoothing according to NSE criterion (default=0.99)

% specification of the magnitude error model
    error_model = 'relative'; % 'relative' or 'standard'; (standard: (sim-obs), relative: [(sim-obs)/((sim + obs)/2)] (default='relative')

% options for time series splitting:
    timeseries_split_by_user = false;   % 'true': time series splits provided by user in ascii file. 'false': splits will be placed by the program (default = 'false')
    split_frequency = 250;              % only required if timeseries_split_by_user='false': this is the default distance between 2 splits (default=500)

% parametrization of the objective function 
    weight_nfc = 1/7;   % weights number of re-assigned hydrological cases (default= 1)    
    weight_rds = 1/7;   % weights the importance of the re-assigned segments (default=1) 
    weight_sdt = 5/7;   % weights the SD timing error component (default=5)
    weight_sdv = 0;     % weights the SD magnitude error component (default=0)
   
% set plot flags 
    pf_input = true;                   % plots input time series ('obs' and 'sim')
    pf_segs_cons_entireTS = true;       % plots obs, sim, colour-coded pairs of matching segments, SeriesDistance connectors for the entire time series
    pf_errordist = true;                % plots SeriesDistance error distributions 

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
     
% replace identical neigbouring values to avoid problems with assignment of unique peaks and valleys
    obs = f_ReplaceEqualNeighbours(obs);
    sim = f_ReplaceEqualNeighbours(sim);
    
% Define time series split points to improve coarse-graining performance
    if timeseries_split_by_user == false
        timeseries_splits = f_FindSplitPoints(obs, sim, split_frequency);   % find split points if they are not provided by the user       
    else
        timeseries_splits = dlmread('data/HOST_ts_splits.csv', ';');        % read splits defined by user.
    end

% plot input data
    if pf_input == true
        f_PlotInput([], obs, [], [], sim, [], [], timeseries_splits);          % show time series splits
    end

% cleanup
    clear smooth_flag nse_smooth_limit pf_input timeseries_split_by_user split_frequency

%% Apply coarse-graining and the SD method to the entire time series 
    % note: contrary to the event based method both, the coarse-graining and the SD calculation 
    % take place in the same function here due to the splitting of the time series. To this end 
    % the splitting is solved in a simplistic way and does not support separating the 
    % coarse-graining and SD calculation as in the event-based version.
    
    % apply coarse graining and SD calculation: determines optimal level of segment aggregation for entire time series and applies SD to it
    [segs_obs_opt_all, segs_sim_opt_all, connectors, e_sd_t_all, e_sd_q_all] = f_CoarseGraining_SD_Continuous(obs, sim, timeseries_splits, ...
        weight_nfc,weight_rds,weight_sdt,weight_sdv,error_model);        
        
    % plot time series with optimized segments and connectors in an own figure
    if pf_segs_cons_entireTS == true
       f_PlotConnectedSeries(obs, segs_obs_opt_all, sim,segs_sim_opt_all, connectors) 
    end
    
    % plot 2d-error distributions
    if pf_errordist == true
        f_PlotSDErrors_OnePanel(e_sd_t_all, e_sd_q_all); 
    end

%% save output
    save (outfile, 'obs', 'sim', 'segs_obs_opt_all', 'segs_sim_opt_all', 'connectors', 'e_sd_t_all', 'e_sd_q_all', ...
                   'weight_nfc', 'weight_rds', 'weight_sdt', 'weight_sdv', 'error_model');   

         
        

        
        
    
    
    





