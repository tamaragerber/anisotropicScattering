function plotProfileFig(linenum,mode)

    % Load data
    % ---------------------------------------------------------------------
    if mode == 0
        load(['output/curveFit_',num2str(linenum),'_20trace_average.mat']);
    elseif mode == 1
        load(['output/curveFit_',num2str(linenum),'_20trace_average_highres.mat']);
    end

    % Create figure
    % ---------------------------------------------------------------------
    figure();
    set(gcf, 'color', 'none','windowstate','maximized');

    % Plot the data
    % ---------------------------------------------------------------------
    if mode ==0 
        plotData(Sav)
        export_fig(['figures/',num2str(linenum)], '-nocrop', '-transparent', '-png');
    elseif mode == 1
        plotdepthProfiles(Sav)
        export_fig(['figures/',num2str(linenum),'_highres'], '-nocrop', '-transparent', '-png');
    end
    
    % save csv files 
    % ---------------------------------------------------------------------
    if mode == 0
        saveCSVfile(Sav,linenum)
    end
end



%% functions


function plotData(Sav)
    
    % Parameters for panel arrangement
    % ---------------------------------------------------------------------
    N = length(Sav.lat);
    ii = 1;
    n = floor(N/2);                                                         % Number of subplots
    gap = [0.02 0.02];                                                      % [horizontal gap, vertical gap]
    marg_h = [0.1 0.05];                                                    % [bottom margin, top margin]
    marg_w = [0.1 0.05];                                                    % [left margin, right margin], increased left margin

    % Create tight subplots
    % ---------------------------------------------------------------------
    ha = tight_subplot(3, n, gap, marg_h, marg_w);
    
    for i = 1:N
        if mod(i, 2) ~= 0
            continue;
        end

        % Plot dPHH
        % -----------------------------------------------------------------
        axes(ha(ii));
        plotSingleData(Sav, Sav.dPHH{1, i}, Sav.time, 'dPHH', i, Sav.curvefit180, Sav.curvefit90, Sav.dp,N);
        xticks([90, 180, 270, 360]);
        xticklabels('')

        % Plot dPHV
        % -----------------------------------------------------------------
        axes(ha(n + ii));
        plotSingleData(Sav, Sav.dPHV{1, i}, Sav.time, 'dPHV', i, Sav.dp,N);
        xticks([90, 180, 270, 360]);
        xticklabels('')
           
        % Plot CHHVV
        % -----------------------------------------------------------------
        axes(ha(2*n + ii));
        plotSingleData(Sav, Sav.CHHVV{1, i}, Sav.time, 'CHHVV', i,Sav.dp, N);
        xticks([90, 180, 270, 360]);

        ii = ii + 1;
    end
end

function plotSingleData(Sav, M, time, dataType, i, curvefit180, curvefit90, dp,N)
    
    imagesc(1:360, time(1:11214), M(1:11214, :));
    hold on;
    axis ij;
    colormap((brewermap([], "RdBu")));
    hold on;

    if nargin > 7
        % Plot fitted curves
        % -----------------------------------------------------------------
        for dpi = 1:length(dp)
            plot(1:360, -curvefit180{dpi, i} * 0.5 + time(dp(dpi)), 'k-', 'linewidth', 2);
            hold on;
            plot(1:360, -curvefit90{dpi, i} * 0.5 + time(dp(dpi)), 'k:', 'linewidth', 2);
            hold on;
        end
    end

    formatAxes(i);
    clim([-4 4]);
end

function formatAxes(i)
    if i > 2
        yticklabels('');
        yticks([11.71, 17.63, 23.55, 29.45]);
    else
        set(gca, 'fontsize', 14);
        yticks([11.71, 17.63, 23.55, 29.45]);
        yticklabels([1000, 1500, 2000, 2500]);
        set(gca, 'TickLabelInterpreter', 'latex');
    end
    box on;
end

function plotdepthProfiles(Sav)
    
    % Parameters for panel arrangement
    % ---------------------------------------------------------------------
    N = length(Sav.lat);
    ii = 1;
    n = N;                                                                  % Number of subplots
    gap = [0.04 0.02];                                                      % [horizontal gap, vertical gap]
    marg_h = [0.05 0.05];                                                   % [bottom margin, top margin]
    marg_w = [0.05 0.02];                                                   % [left margin, right margin], increased left margin

    % Create tight subplots
    % ---------------------------------------------------------------------
    ha = tight_subplot(2, n, gap, marg_h, marg_w);
  
  
    for i = 1:N
        % Plot dPHH
        axes(ha(ii));
        plot(Sav.maxAmp90{1,i},Sav.depth(Sav.dp),'b','linewidth',2);hold on
        plot(Sav.maxAmp180{1,i},Sav.depth(Sav.dp),'r','linewidth',2);
        axis ij
        yticks([1000, 1500, 2000, 2500])
        if i>1
          yticklabels('')
        end
        box on
        maxXValue = max([Sav.maxAmp90{1,i},Sav.maxAmp180{1,i}]);
        if maxXValue <= 5
            upperXLimit = 5;
        elseif maxXValue <= 10
            upperXLimit = 10;
        elseif maxXValue <= 15
            upperXLimit = 15;
        else
            upperXLimit = 20;
        end

        % Set the x-axis limits
        % -----------------------------------------------------------------
        xlim([0 upperXLimit]);
        ylim([Sav.depth(Sav.dp(1)) Sav.depth(Sav.dp(end))])
        xticks([0,upperXLimit/2, upperXLimit])
        grid on
        set(gca, 'fontsize',12)
       
        % plot r2
        
        axes(ha(n+ii))
        plot(Sav.rsquare90{1,i},Sav.depth(Sav.dp),'b','linewidth',2);hold on
        plot(Sav.rsquare180{1,i},Sav.depth(Sav.dp),'r','linewidth',2);hold on
        axis ij
      
        if ii>1
          yticklabels('')
        end
        xticks([0, 0.5, 1])
        yticks([1000, 1500, 2000, 2500])
        axis([0 1 Sav.depth(Sav.dp(1)) Sav.depth(Sav.dp(end))]);
        box on
        grid on
        set(gca, 'fontsize',12)

        ii = ii + 1;
    end
end