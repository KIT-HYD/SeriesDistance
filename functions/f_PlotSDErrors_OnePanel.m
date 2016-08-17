function [ ] = f_PlotSDErrors_OnePanel(t_errors, q_errors)
% Plots 2-d error distributions
    
    ff=figure('units','centimeters','position',[1,1,8,8]);
    set(gcf, 'PaperSize', [8 8]);
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    hold on;    
    % magnitude errors
        % plot the error distribution
        scatter1 =scatter(t_errors, q_errors);
        set(scatter1,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75],'Marker','.');
        
        % set x- and ylim, but only if there is at least one non-NaN value in the plot
        if any(~isnan(t_errors));
            ylim([nanmax(nanmax(q_errors), abs(nanmin(q_errors)))*-1.1, nanmax(nanmax(q_errors), abs(nanmin(q_errors)))*1.1]);
            if(all(q_errors==0))
                xlim([-0.1,0.1])
            else
                xlim([nanmax(nanmax(t_errors), abs(nanmin(t_errors)))*-1.1, nanmax(nanmax(t_errors), abs(nanmin(t_errors)))*1.1]);    
            end
        end
        
        % add vertical and horizontal lines
        line(xlim, [0 0], 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
        line([0 0],ylim, 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
        
        % plot the mean (center) of the error distribution
        plot(nanmean(t_errors),nanmean(q_errors),'MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',10,'Marker','+','LineStyle','none');

        % add labels
        xlabel('timing error (>0: obs later sim)');
        ylabel('magnitude error (>0: obs larger sim)');
        %titlestr = strcat('2d error distribution for rising segments, Q range:', {' '}, num2str(error_lvls(i)), {' '},'<= Qsim <', {' '},num2str(error_lvls(i+1))); 
        title('entire time series');

        % format the figure
        axis square;
        box on;
        hold off;
end
