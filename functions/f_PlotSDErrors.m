function [ ] = f_PlotSDErrors(e_sd_q_rise, e_sd_t_rise, e_sd_q_fall, e_sd_t_fall, e_sd_q_rise_subset, e_sd_t_rise_subset, e_sd_q_fall_subset, e_sd_t_fall_subset)
% Plots 2-d error distributions
    
    ff=figure('units','centimeters','position',[1,1,16,8]);
    set(gcf, 'PaperSize', [16 8]);
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    
    hold on;     
    % magnitude errors
        subplot(1,2,1);
        hold on;
        % plot the error distribution
        scatter1 =scatter(e_sd_t_rise, e_sd_q_rise);
        set(scatter1,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75],'Marker','.');
        
        % set x- and ylim, but only if there is at least one non-NaN value in the plot
        if any(~isnan(e_sd_q_rise));
            ylim([nanmax(nanmax(e_sd_q_rise), abs(nanmin(e_sd_q_rise)))*-1.1, nanmax(nanmax(e_sd_q_rise), abs(nanmin(e_sd_q_rise)))*1.1]);
            if(all(e_sd_t_rise==0))
                xlim([-0.1,0.1])
            else
                xlim([nanmax(nanmax(e_sd_t_rise), abs(nanmin(e_sd_t_rise)))*-1.1, nanmax(nanmax(e_sd_t_rise), abs(nanmin(e_sd_t_rise)))*1.1]);    
            end
        end
        
        % add vertical and horizontal lines
        line(xlim, [0 0], 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
        line([0 0],ylim, 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
        
        % if specified, plot subset error distribution
        if ~isempty(e_sd_q_rise_subset) && ~isempty(e_sd_t_rise_subset)
            
            scatter2 =scatter(e_sd_t_rise_subset,e_sd_q_rise_subset);
            set(scatter2,'MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'Marker','.');
        end

        % plot the mean (center) of the error distribution
        plot(nanmean(e_sd_t_rise),nanmean(e_sd_q_rise),'MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',10,'Marker','+','LineStyle','none');

        % add labels
        xlabel('timing error (>0: obs later sim)');
        ylabel('magnitude error (>0: obs larger sim)');
        %titlestr = strcat('2d error distribution for rising segments, Q range:', {' '}, num2str(error(i)), {' '},'<= Qsim <', {' '},num2str(error(i+1))); 
        title('rising segments');

        % format the figure
        axis square;
        box on;
        hold off;
    
%%%% add second plot on timing errors
        subplot(1,2,2);
        hold on;

        % plot the error distribution
        scatter1 =scatter(e_sd_t_fall,e_sd_q_fall);
        set(scatter1,'MarkerFaceColor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75],'Marker','.');
        
        % set x- and ylim, but only if there is at least one non-NaN value in the plot
        if any(~isnan(e_sd_q_fall));
            ylim([nanmax(nanmax(e_sd_q_fall), abs(nanmin(e_sd_q_fall)))*-1.1, nanmax(nanmax(e_sd_q_fall), abs(nanmin(e_sd_q_fall)))*1.1]);
            if(all(e_sd_t_rise==0))
                xlim([-0.1,0.1])
            else
                xlim([nanmax(nanmax(e_sd_t_rise), abs(nanmin(e_sd_t_rise)))*-1.1, nanmax(nanmax(e_sd_t_rise), abs(nanmin(e_sd_t_rise)))*1.1]);    
            end
        end
        
        % add vertical and horizontal lines
        line(xlim, [0 0], 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
        line([0 0],ylim, 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
        
        % if specified, plot subset error distribution
        if ~isempty(e_sd_q_fall_subset) && ~isempty(e_sd_t_fall_subset)
            scatter2 = scatter(e_sd_t_fall_subset,e_sd_q_fall_subset);
            set(scatter2,'MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3],'Marker','.');            
        end
        
        % plot the mean (center) of the error distribution
        plot(nanmean(e_sd_t_fall),nanmean(e_sd_q_fall),'MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',10,'Marker','+','LineStyle','none');


        % add labels
        xlabel('timing error (>0: obs later sim)');
        ylabel('magnitude error (>0: obs larger sim)');  
        %titlestr = strcat('2d error distribution for falling segments, Q range:', {' '}, num2str(error(i)), {' '},'<= Qsim <', {' '},num2str(error(i+1))); 
        title('falling segments');

        % format the figure
        axis square;
        box on;
        hold off;    
end

