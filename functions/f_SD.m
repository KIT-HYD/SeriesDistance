function [e_q, e_t, e_ysim, e_q_rise, e_t_rise, e_ysim_rise, e_q_fall, e_t_fall, e_ysim_fall, cons, e_rise_MD, e_fall_MD] = f_SD(y_obs,segs_obs,y_sim,segs_sim,error_model, printflag)
% Calculates the distance vectors in time and value between two matching events (obs/sim)

% Modification history
% - 2015/08/03: Uwe Ehret, Simon Seibert: First version
% - 2015/09/16: Uwe Ehret: Included possibility for seperate error distributions for discharge ranges

% OUTPUT
    % e_q           % magnitude errors for both rise and fall, e_q=[e_q_rise; e_q_fall]
    % e_t           % timing errors for both rise and fall,    e_t=[e_t_rise; e_t_fall]
    % e_ysim        % 
    % e_q_rise      % magnitude errors in rising limbs
    % e_t_rise      % time errors in rising limbs
    % e_ysim_rise
    % e_q_fall      % magnitude errors in falling limbs
    % e_t_fall      % time errors in falling limbs
    % e_ysim_fall
    % cons          % SD connectors
    % e_rise_MD     % 1D magnitude errors of corresponding rising limb sections
    % e_fall_MD     % 1D magnitude errors of corresponding falling limbs sections

% OLD! OUTDATED!
% INPUT
    % x_org: (n,1) matrix with time position (x-position) of observed values
    % y_obs: (n,1) matrix with observed values
    % segs_obs: (1,x) struct with observed segments, where x is the number of segments (same as in segs_sim!!
    % x_org: (m,1) matrix with time position (x-position) of simulated values (n can be ~= m)
    % y_sim: (m,1) matrix with simulated values
    % segs_sim: (1,x) struct with simulated segments, where x is the number of segments (same as in segs_obs!!) 
    % vert_dist_fun: (1) sets the way the magnitude distance among obs and sim is computed
    %               if 'relative': dist_v = (obs - sim) / ((obs + sim)*0.5)
    %               if 'standard': dist_v = (obs - sim)
    %               Default: 'standard'
    % savepointsflag: (1) if true, the connector points of series distance will be stored in the global variables *_match_*
    
    % METHOD
    % Note: 
    % - The # of obs and sim segments needs to be equal
    % - The order of the segment types of obs and sim has to be equal: either both start with a 'rise' or a 'fall'
    % the total number of connectors for the event equals mean(length(obs_event),length(sim_event))
    % the number of connectors per segment is determined by the mean importance of the segment (mean of obs and sim relevance)

    if ~exist('printflag','var')
        printflag = 'false';
    end
    
cons = struct('x_match_obs_global',{},'y_match_obs',{},'x_match_sim_global',{},'y_match_sim',{} );

%% specify connectors

    % determine the total number of connectors (average of total length of obs and sim event times the user-specified percentage)
    totnumcons = (length(y_obs) + length(y_sim)) * 0.5;

    % determine the share of connectors for each segment
    num_segs = length(segs_obs);  % # of segments (could also be length (segs_sim)
    segs_cons=NaN(1,num_segs);      % variable for the number of connectors assigned to each segment
    sum_rels=sum([segs_obs.relevance])+sum([segs_sim.relevance]);   % the overall sum of relevance (as relevance is already normalized, should be 1 + 1 = 2

    % the share of connectors for each segment is proportional to its relative relevance
    segs_cons = round(([segs_obs.relevance] + [segs_sim.relevance]) * totnumcons / sum_rels);

%% initialize output variables   

    e_q_rise = [];      % magnitude errors in rising limbs
    e_t_rise = [];      % time errors in rising limbs
    e_q_fall = [];      % magnitude errors in falling limbs
    e_t_fall = [];      % time errors in falling limbs
    e_ysim_rise = [];   % simulated discharge for rising limbs, correspoding to each error (needed to subdivide errors in discharge classes)
    e_ysim_fall = [];   % simulated discharge for falling limbs, correspoding to each error (needed to subdivide errors in discharge classes)
    e_rise_MD = [];     % 1D magnitude errors in corresponding rising limb sections 
    e_fall_MD = [];     % 1D magnitude errors in corresponding falling limb sections 

    % loop over all segments
    for z = 1 : num_segs

        % determine the GLOBAL x-location (time) of the connectors in the current segment  
        num = segs_cons(z);      
        con_x_obs_global_seg = linspace(segs_obs(z).starttime_global,segs_obs(z).endtime_global,num);
        con_x_sim_global_seg = linspace(segs_sim(z).starttime_global,segs_sim(z).endtime_global,num);  

        % determine the LOCAL x-location (time) of the connectors in the current segment     
        con_x_obs_local_seg = linspace(segs_obs(z).starttime_local,segs_obs(z).endtime_local,num);
        con_x_sim_local_seg = linspace(segs_sim(z).starttime_local,segs_sim(z).endtime_local,num);      

        % determine the local x-locations of the segment
        x_obs_local_seg = segs_obs(z).starttime_local : segs_obs(z).endtime_local;
        x_sim_local_seg = segs_sim(z).starttime_local : segs_sim(z).endtime_local;    

        xobs = segs_obs(z).starttime_local:segs_obs(z).endtime_local;
        xsim = segs_sim(z).starttime_local:segs_sim(z).endtime_local;
        xint = intersect(xobs, xsim);

        % show vertically compared segments 
        if false
            plot(xint , y_obs(xint), 'k','LineWidth',2)
            hold on;
            plot(xint , y_sim(xint), 'k','LineWidth',2)
            plot(1:length(y_obs), y_obs, 'r')
            plot(1:length(y_sim), y_sim, 'b')
        end

        % determine the local y-values of the segment (based on the local x-locations)
        y_obs_seg = y_obs(x_obs_local_seg);
        y_sim_seg = y_sim(x_sim_local_seg);  

        % determine the y-values of the connectors in the current segment (based on the local x-locations) with linear interpolation
        con_y_obs_seg = interp1(x_obs_local_seg,y_obs_seg,con_x_obs_local_seg,'linear'); 
        con_y_sim_seg = interp1(x_sim_local_seg,y_sim_seg,con_x_sim_local_seg,'linear');     

        % find out whether the current segment is 'rise' or 'fall'
        if segs_obs(z).sum_dQ > 0   % rise
            % calculate the length of the connectors (distance between connector points on obs and sim) in the current segment    
                % time (x) distances 
                e_t_rise_seg = (con_x_obs_global_seg - con_x_sim_global_seg);       % > 0 means obs is later than sim

                % magnitude distances
                if strcmp(error_model,'standard') % compute the simple difference 
                    e_q_rise_seg = (con_y_obs_seg - con_y_sim_seg); % > 0 means obs is larger than sim
                elseif strcmp(error_model,'relative') % compute a scaled difference
                    e_q_rise_seg = (con_y_obs_seg - con_y_sim_seg) ./ ((con_y_obs_seg + con_y_sim_seg) .* 0.5); % > 0 means obs is larger than sim
                else
                    error('distance function not properly specified');
                end

                % add the errors of the segment to the overall errors of the event
                e_q_rise = [e_q_rise; e_q_rise_seg'];
                e_t_rise = [e_t_rise; e_t_rise_seg'];  
                e_ysim_rise = [e_ysim_rise; con_y_sim_seg'];

                % add vertical 1D error to output array
                e_rise_MD = [e_rise_MD; (y_obs(xint)-y_sim(xint)) ./ ((y_obs(xint)+y_sim(xint)) .* 0.5)];
        else                        % fall
            % calculate the length of the connectors (distance between connector points on obs and sim) in the current segment    
                % time (x) distances
                e_t_fall_seg = (con_x_obs_global_seg - con_x_sim_global_seg);   % > 0 means obs is later than sim

                % magnitude distances
                if strcmp(error_model,'standard') % compute the simple difference
                    e_q_fall_seg = (con_y_obs_seg - con_y_sim_seg);   % > 0 means obs is larger than sim
                elseif strcmp(error_model,'relative') % compute a scaled difference
                    e_q_fall_seg = (con_y_obs_seg - con_y_sim_seg) ./ ((con_y_obs_seg + con_y_sim_seg) .* 0.5);   % > 0 means obs is larger than sim
                else
                    error('distance function not properly specified');
                end      

                % add the errors of the segment to the overall errors of the event
                e_q_fall = [e_q_fall; e_q_fall_seg'];
                e_t_fall = [e_t_fall; e_t_fall_seg']; 
                e_ysim_fall = [e_ysim_fall; con_y_sim_seg'];

                % add vertical 1D error to output array
                e_fall_MD = [e_fall_MD; (y_obs(xint)-y_sim(xint)) ./ ((y_obs(xint)+y_sim(xint)) .* 0.5)];

        end 

        % add matching points(x,y) of the segment to overall matching points(for plotting)
        if isempty(cons)
            cons(1).x_match_obs_global =  con_x_obs_global_seg;
            cons(1).y_match_obs = con_y_obs_seg;
            cons(1).x_match_sim_global = con_x_sim_global_seg;
            cons(1).y_match_sim = con_y_sim_seg; 
        else
            cons(1).x_match_obs_global = [cons(1).x_match_obs_global con_x_obs_global_seg];
            cons(1).y_match_obs = [cons(1).y_match_obs con_y_obs_seg];
            cons(1).x_match_sim_global = [cons(1).x_match_sim_global con_x_sim_global_seg];
            cons(1).y_match_sim = [cons(1).y_match_sim  con_y_sim_seg]; 
        end

    end

    % combine all case-specific error distributions to one for magnitude and one for time
    e_q = [e_q_rise; e_q_fall];
    e_t = [e_t_rise; e_t_fall];
    e_ysim = [e_ysim_rise; e_ysim_fall];

