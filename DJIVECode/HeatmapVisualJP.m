function HeatmapVisualJP(mdata, dataname, fig_pos, iprint, figdir, figname, limup, limd)
% HeatmapVisualMJ   Summary of this function goes here
%   Detailed explanation goes here
%
% Inputs:
%   X - n x d data matrix
%   dataname - nb x 1 cell array of string of data matrix's name
%
%   Copyright (c)  Meilei Jiang, Jack Prothero 2019

    if nargin == 1
        dataname = 'Data'; % set default value if no input
    end

    % set the color scheme 
    b=0.99:-0.01:0;
    bb=0:0.01:0.99;
    blpha=[ones(length(b),1)*1 b' b'];bblpha=[bb' bb' ones(length(bb),1)*1 ];
    map=[bblpha;[1 1 1]; blpha;];

    %vax = axisSM(mdata, (mdata)') ;
    if nargin <= 5
        limup= max(prctile(mdata(:), 95), prctile(-mdata(:), 95));
        limd = min(prctile(mdata(:), 5), prctile(-mdata(:), 5));
    end

    %
    fig = figure;
    if exist('fig_pos', 'var')
        set(fig,'Position', fig_pos)
    end
    %}
    imagesc(mdata);
    ax = get(gca);
    colormap(map);
    caxis([limd, limup]);
    %caxis([min(vax), max(vax)]);
    set(gca, 'XTick', []) ; % remove X-axis
    c=colorbar('southoutside', 'FontSize', 10); %add colorbar
    title(dataname, 'FontSize', 18, 'Interpreter', 'Latex',...
        'FontWeight', 'bold') ;  % Add title for subplot
    % save the figure
    if exist('iprint', 'var') && iprint == 1
        if ~exist('figdir', 'var') || ~exist(figdir, 'dir')
            disp('No valid figure directory found! Will save to current folder.')
            figdir = '';
        end
        
        if ~exist('figname', 'var') || isempty(figname)
            figname = strcat(['Data_Heatmap']);
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
end

