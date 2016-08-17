function [ ] = f_PlotInput(obs_org, obs, obs_events, sim_org, sim, sim_events, obs_sim_events_mapped, timeseries_splits)
% Plots time series, events, events mappings and time series splits

% initialize
legendtext = cell(1,1);
nextlegendentry = 1;
h = [];
figure;
hold on;

% if specified, plot time series of obs
if ~isempty(obs)
    xes = 1:1:length(obs);
    h = [h, plot(xes,obs,'-b')];
    legendtext (nextlegendentry,1) = cellstr('obs'); nextlegendentry = nextlegendentry +1;     
end

% if specified, plot time series of obs_org
if ~isempty(obs_org)
    xes = 1:1:length(obs);
    h = [h, plot(xes,obs_org,'--b')];
    legendtext (nextlegendentry,1) = cellstr('obs_org'); nextlegendentry = nextlegendentry +1;     
end

% if specified, plot time series of sim
if ~isempty(sim)
    xes = 1:1:length(sim);
    h = [h, plot(xes,sim,'-r')];
    legendtext (nextlegendentry,1) = cellstr('sim'); nextlegendentry = nextlegendentry +1;     
end

% if specified, plot time series of sim_org
if ~isempty(sim_org)
    xes = 1:1:length(sim);
    h = [h, plot(xes,sim_org,'--r')];
    legendtext (nextlegendentry,1) = cellstr('sim_org'); nextlegendentry = nextlegendentry +1;     
end

% if specified, plot obs events
if ~isempty(obs_events)
    for i = 1 : size(obs_events,1)
        xes_obs = (obs_events(i,1):1:obs_events(i,2));
        if i == 1 
           h = [h, plot(xes_obs,obs(xes_obs),'-b','LineWidth',4)]; 
           legendtext (nextlegendentry,1) = cellstr('obs events'); nextlegendentry = nextlegendentry +1;  
        else
            plot(xes_obs,obs(xes_obs),'-b','LineWidth',4); 
        end
    end   
end

% if specified, plot sim events
if ~isempty(sim_events)
    for i = 1 : size(sim_events,1)
        xes_sim = (sim_events(i,1):1:sim_events(i,2));
        if i == 1
            h = [h, plot(xes_sim,sim(xes_sim),'-r','LineWidth',4)];  
            legendtext (nextlegendentry,1) = cellstr('sim events'); nextlegendentry = nextlegendentry +1;  
        else
            plot(xes_sim,sim(xes_sim),'-r','LineWidth',4);  
        end         
    end  
end

% if specified, plot event connections
if ~isempty(obs_sim_events_mapped)
    for i = 1 : size(obs_sim_events_mapped,1)
                   
        % event starts
        xes_connect = [obs_sim_events_mapped(i,1) obs_sim_events_mapped(i,2)];
        yes_connect = [obs(xes_connect(1)) sim(xes_connect(2))];
        plot(xes_connect,yes_connect,'-k','LineWidth',1);

        % event ends
        x_obs_end = obs_events(find(obs_events(:,1) == obs_sim_events_mapped(i,1)) ,2);
        x_sim_end = sim_events(find(sim_events(:,1) == obs_sim_events_mapped(i,2)) ,2);
        xes_connect = ([x_obs_end x_sim_end]);
        yes_connect = [obs(xes_connect(1)) sim(xes_connect(2))];
        plot(xes_connect,yes_connect,'-k','LineWidth',1);
        
    end
    h = [h, plot([1,1],[1,1],'-k','LineWidth',1)];
    legendtext (nextlegendentry,1) = cellstr('connected events'); nextlegendentry = nextlegendentry +1;         
end

% if specified, plot time series splits
if ~isempty(timeseries_splits)
    for i = 1 : size(timeseries_splits,1)
         xes = [timeseries_splits(i) timeseries_splits(i)];
         yes = [min([obs; sim]) max([obs; sim])];
         line(xes, yes, 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
    end
    h = [h, plot([1,1],[1,1],'-k','LineWidth',1)];
    legendtext (nextlegendentry,1) = cellstr('timesplits'); nextlegendentry = nextlegendentry +1; 
end

% plot legend    
legend(h,legendtext);
hold off;    
    
end

