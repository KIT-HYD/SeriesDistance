function [ ] = f_PlotConnectedSeries(obs,segs_obs,sim,segs_sim,connectors)
% plots two related time series (obs and sim) and their series-distance connectors
% 19.Nov.2013 Uwe Ehret

if ~exist('showEventIndex', 'var'); % adds event numbers to the plot
    showEventIndex = false;
end
show_connectors = true;

%distance vectors between matching points obs/sim
if show_connectors == true 
    u = connectors.x_match_sim_global - connectors.x_match_obs_global;
    v = connectors.y_match_sim - connectors.y_match_obs;
end

ff = figure();
hold on;

% Plot the timeseries
plot(1:1:size(obs,1),obs,'Color', [0.2 0.2 0.2]);
plot(1:1:size(sim,1),sim,'--','Color', [0.2 0.2 0.2]);

% Plot Feature Distance lines
%distance vectors between matching points obs/sim
if show_connectors == true 
   quiver(connectors.x_match_obs_global,connectors.y_match_obs,u,v,0,'k','ShowArrowHead','off','Color', [0.6 0.6 0.6]);   
end

% Plot the connected segments in unnique color
num_segs = length(segs_obs);
cmap=[1 0 0;0 0.4 1;0 1 0;0.8 0 1;1 0.8 0.5]; % colormap

cmap_count = 1;
for z = 1: num_segs
    xes_global = (segs_obs(z).starttime_global : segs_obs(z).endtime_global);
    yes_global = obs(xes_global);
    plot(xes_global,yes_global,'-','Color',cmap(cmap_count,:),'LineWidth',2);
    
    xes_global = (segs_sim(z).starttime_global : segs_sim(z).endtime_global);
    yes_global = sim(xes_global);
    plot(xes_global,yes_global,'--','Color',cmap(cmap_count,:),'LineWidth',2); 
    
    cmap_count = cmap_count + 1;
    if cmap_count > length(cmap)
        cmap_count = 1;
    end
end

% Formatting
    if show_connectors == false
        legend('observation','simulation','Location','NorthEast');
        legend('boxoff');
    else
        legend('observation','simulation','connectors','Location','NorthEast');
        legend('boxoff');
    end    
ax1 = gca;
set(ax1,'Fontsize',11,'FontWeight','bold')    

% add event index to plot
if showEventIndex==true
    event_ind = [1 find(diff([segs_obs(:).eventID])==1)+1]; % row index
    for kk = 1:1:length(event_ind)
        text([segs_obs(event_ind(kk)).starttime_global], 0, strcat('# ', num2str(kk)))
    end
end
box on;
hold off


end

