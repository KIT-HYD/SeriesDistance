function [ ] = f_Plot1dErrors(e_sd_q_1d, e_sd_q_1d_subset, titlestring )
% Plots 1-d error distributions

% Modification history
% - 2015/08/01:  Uwe Ehret and Simon Seibert: First version
% - 2015/09/16: Uwe Ehret: Included possibility for seperate error distributions for discharge ranges

    ff=figure('units','centimeters','position',[1,1,8,8]);
    set(gcf, 'PaperSize', [8 8]);
    pos = get(gca,'Position');
    ti = get(gca,'TightInset');
    hold on; 
 
    % plot the error distribution
    cdf_e_sd_q_1d = linspace(0,1,length(e_sd_q_1d));
    e_sd_q_1d_sorted = sort(e_sd_q_1d);
    plot(e_sd_q_1d_sorted, cdf_e_sd_q_1d, 'MarkerFaceColor', [0.75 0.75 0.75],....
        'MarkerEdgeColor',[0.75 0.75 0.75],'Marker','.');

    % plot the mean (center) of the error distribution
    e_mean = mean(e_sd_q_1d);
    pu_e_mean = cdf_e_sd_q_1d(find((e_sd_q_1d_sorted >= e_mean), 1));
    plot(e_mean, pu_e_mean, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 10, 'Marker', '+', 'LineStyle', 'none');

    % if specified, plot subset error distribution
%     if ~isempty(e_sd_q_1d_subset) 
%         find the pu-values for all points in the subset
%         indx = [];
%         num_errors_in_subset = sum(~isnan(e_sd_q_1d_subset)); 
%         pos = find((e_sd_q_1d_sorted == e_sd_q_1d_subset),1);
%         pos = find(ismember(e_sd_q_1d_sorted,e_sd_q_1d_subset)==1);
%         indx = [indx pos];
%         plot (e_sd_q_1d_subset,cdf_e_sd_q_1d(indx),'.r');
%     end

    % set x-lim, but only if there is at least one non-NaN value in the plot
    if any(~isnan(e_sd_q_1d_sorted))
        xlim([max(max(e_sd_q_1d_sorted), abs(min(e_sd_q_1d_sorted)))*-1.1, max(max(e_sd_q_1d_sorted), abs(min(e_sd_q_1d_sorted)))*1.1]);
    end

    % Add vertical and horizontal lines and format axes
    line([0 0], [0 1], 'LineStyle', '-.', 'LineWidth',1, 'Color',[0.65 0.65 0.65]);
    xlabel('magnitude error (>0: obs larger sim)');
    ylabel('cum. prob. of unexceedance');
    title(titlestring);    
    axis square;
    box on;
    hold off;    
    
end

