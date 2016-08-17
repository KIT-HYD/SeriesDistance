function [ ] = f_PlotCoarseGrainIntSteps(obs,segs_obs,sim,segs_sim,connectors,titlestring)
% plots two related time series (obs and sim) and their series-distance connectors
% 19.Nov.2013 Uwe Ehret

show_connectors = true;

%distance vectors between matching points obs/sim
if show_connectors == true 
    u = connectors.x_match_sim_global - connectors.x_match_obs_global;
    v = connectors.y_match_sim - connectors.y_match_obs;
end

ff = figure();
hold on;

% Plot the timeseries
plot(segs_obs(1).starttime_global:segs_obs(end).endtime_global,obs,'-b');
plot(segs_sim(1).starttime_global:segs_sim(end).endtime_global,sim,'--r');

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
    plot(xes_global,obs(segs_obs(z).starttime_local:segs_obs(z).endtime_local),'-','Color',cmap(cmap_count,:),'LineWidth',2);
    
    xes_global = (segs_sim(z).starttime_global : segs_sim(z).endtime_global);
    plot(xes_global,sim(segs_sim(z).starttime_local:segs_sim(z).endtime_local),'--','Color',cmap(cmap_count,:),'LineWidth',2); 
    
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
title(titlestring)
box on; 
hold off


end

