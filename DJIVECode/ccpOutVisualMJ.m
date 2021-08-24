function ccpOutVisualMJ(angleHats, phiBars, dataname, iprint, figdir, figname, angleTrues)
% ccpOutVisualMJ   Visualize the optimization results of a specific
% direction.
%   Detailed explanation goes here
%
% Inputs:
%   angleHats - nb x 1 cell array of projected angles on each data block
%   phiBars - vector of perturbation angle for each data matrix
%   dataname - a nb x 1 cell array of string of data matrix's name
%   iprint - indicator of whether saving the results.
%   figdir - a string of saving path
%   figname - astring of figure name
%
%
%   Copyright (c)  Meilei Jiang 2018
    if ~exist('figname', 'var') || isempty(figname)
            figname = strcat('opt_progress');
    end
        
    nb = length(phiBars);
    T = length(angleHats{1});
    idx = 1:T;
    fig = figure;
    set(fig,'Position', [1, 1, 1500, 500])
    set(gca, 'FontSize',30)   
    
    for ib = 1:nb
        subplot(1, nb, ib)
        plot(idx, angleHats{ib}, ...
            'DisplayName', strcat('$\hat{\theta}_', num2str(ib), '$'),...
            'LineWidth',2)
        hold on
        %{
        plot(idx, angleTrues{ib}, 'r--',...
            'DisplayName', strcat('$\theta_', num2str(ib), '$'),...
            'LineWidth',2)
        %}
        hline = refline([0, phiBars(ib)]);
        hline.Color = 'g';
        hline.LineStyle = '-.';
        hline.DisplayName = strcat('$\phi_', num2str(ib),'$');
        hline.LineWidth = 2;
        xlim([1, T]);
        ylim([0, inf])
        %Solve Legend Issue Later
        %{
        legend(strcat('$\hat{\theta}_', num2str(ib), '$'), strcat('$\phi_', num2str(ib),'$')) 
        lgd1 = legend;
        lgd1.Interpreter = 'Latex';
        lgd1.FontSize = 30;
        lgd1.FontWeight = 'bold';
        %}
        title({strcat('{\bf\fontsize{30} ',dataname{ib}, '}'); ...
            strcat('{\it\fontsize{14} ',figname, '}') },...
            'FontWeight','Normal')
    end
    
    % save the figure
    if exist('iprint', 'var') && iprint == 1
        if ~exist('figdir', 'var') || ~exist(figdir, 'dir')
            disp('No valid figure directory found! Will save to current folder.')
            figdir = '';
        end
        
        savestr = strcat(figdir, figname);        
        try
            % orient landscape
            print(fig, '-dpng', savestr)
            disp('Save figure successfully!')
        catch
            disp('Fail to save figure!')
        end
    end   
    
    close
end

