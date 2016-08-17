function [] = f_plot_ObjectiveFunction_CoarsGrainStps(CoarseGrainStps, opt_step, titlestring, savefigure)
   
    if ~exist('savefigure', 'var');
    savefigure = false;
    end
    
    % plot the results of the objective function for all aggregation steps
    hfig = figure('units','centimeters','position',[1,1,8,8]);
    hold on;

    x=(1:length(CoarseGrainStps))';  % create x-values
    rr = plot (x,CoarseGrainStps, 'ko');
    set(rr, 'MarkerFaceColor', [0.85 0.85 0.85], 'MarkerSize', 6)
    ylim([0, max(1,max(CoarseGrainStps)*1.1)    ])
    set(gca,'fontsize',10)
    %lbl = num2cell((0:length(CoarseGrainStps)-1));
    box on;      
    
    % highlight the best
    oo = plot (opt_step,CoarseGrainStps(opt_step),'ko');
    set(oo, 'MarkerFaceColor', 'r', 'MarkerSize', 6)
    
    xlabel({'coarse graining step (1=initial cond.)'});
    ylabel({'objective function value (-)'});

    title(titlestring)
    hold off;    
    
    % save output file
    if savefigure == true
        saveas(oo,horzcat('./results/results_of_aggregation.emf'));
    end
end