function simulation7Plots(cosdelphi,cosdelphi_grd,rmse_MEMD,rmse_EWT,rmse_IRCNN,rmse_MVMD,rmse_MSWD,t,varargin)

figure
placeholder = [1,6,11,16,21,7,12,17,22,13,18,23,19,24,25];
h = tight_subplot(5,5,[0.01 0.02],[.05 .001],[.05 .1]);
placeholder_del = [2,3,4,5,8,9,10,14,15,20];
for i = 1:length(placeholder_del)
    delete(h(placeholder_del(i)))
end
for i = 1:length(placeholder)
    for j = 1:length(cosdelphi)
        cosdelphi1(j,:) = cosdelphi{j}(:,i);
    end
    [cosdelphi_mean,cosdelphi_low,cosdelphi_up] = confidence_interval(cosdelphi1);
    axes(h(placeholder(i)));
    plot(t,cosdelphi_mean,"Color",'b','LineWidth',1)
    plot(t,cosdelphi_grd(:,i),"Color",'r','LineWidth',1)
    hold on
    plot(t,cosdelphi_low,"Color",'b')
    plot(t,cosdelphi_up,"Color",'b')
    
    patch([t fliplr(t)],[cosdelphi_low fliplr(cosdelphi_up)],'b', 'FaceAlpha',0.5, 'EdgeColor','none')
    xlim([0,max(t)])
    ylim([-1.5,1.5])
    hold off
    if placeholder(i) ~= 1 && placeholder(i) ~= 6 && placeholder(i) ~= 11 && placeholder(i) ~= 16 && placeholder(i) ~= 21
        set(gca, 'YTick', [])
    end
    if placeholder(i) ~= 21 && placeholder(i) ~= 22 && placeholder(i) ~= 23 && placeholder(i) ~= 24 && placeholder(i) ~= 25
        set(gca, 'XTick', [])
    end
end
if ~isempty(varargin)
    exportgraphics(gcf,[varargin{1},'mean_states.png'],'Resolution',600)
end

if ~isempty(rmse_MEMD) && ~isempty(rmse_EWT) && ~isempty(rmse_IRCNN) && ~isempty(rmse_MVMD) && ~isempty(rmse_MSWD)
    %% MVMD-MSWD
    cmap = [0.2422,0.1504,0.6603;...
        0.2794,0.2653,0.9094;...
        0.2621,0.4088,0.9946;...
        0.1755,0.5554,0.9473;...
        0.1009,0.6750,0.8653;...
        0.1023,0.7510,0.7088;...
        0.3176,0.7994,0.4975;...
        0.6701,0.7796,0.2239;...
        0.9523,0.7284,0.2237;...
        0.9694,0.8591,0.1665;...
        0.9769,0.9839,0.0805];
    figure
    h = tight_subplot(5, 5, [0 0],[0 0],[0 0.1]);
    placeholder_MVMD = [1,6,11,16,21,7,12,17,22,13,18,23,19,24,25];
    placeholder_MSWD = [1,2,3,4,5,7,8,9,10,13,14,15,19,20,25]; 
    for i = 1:length(placeholder_MVMD)
        if rmse_MVMD(i) < 0.0454
            color = cmap(1,:);
        elseif rmse_MVMD(i) >= 0.0454 && rmse_MVMD(i) < 0.0908
            color = cmap(2,:);
        elseif rmse_MVMD(i) >= 0.0908 && rmse_MVMD(i) < 0.1362
            color = cmap(3,:);
        elseif rmse_MVMD(i) >= 0.1362 && rmse_MVMD(i) < 0.1816
            color = cmap(4,:);
        elseif rmse_MVMD(i) >= 0.1816 && rmse_MVMD(i) < 0.227
            color = cmap(5,:);
        elseif rmse_MVMD(i) >= 0.227 && rmse_MVMD(i) < 0.2724
            color = cmap(6,:);
        elseif rmse_MVMD(i) >= 0.2724 && rmse_MVMD(i) < 0.3178
            color = cmap(7,:);
        elseif rmse_MVMD(i) >= 0.3178 && rmse_MVMD(i) < 0.3632
            color = cmap(8,:);
        elseif rmse_MVMD(i) >= 0.3632 && rmse_MVMD(i) < 0.4086
            color = cmap(9,:);
        elseif rmse_MVMD(i) >= 0.4086 && rmse_MVMD(i) < 0.454
            color = cmap(10,:);
        elseif rmse_MVMD(i) >= 0.454 && rmse_MVMD(i) < 0.5
            color = cmap(11,:);
        end
        axes(h(placeholder_MVMD(i)));
        if placeholder_MVMD(i) == 1 || placeholder_MVMD(i) == 7 ||...
                placeholder_MVMD(i) == 13 || placeholder_MVMD(i) == 19 || placeholder_MVMD(i) == 25
            x = [0,1,0];
            y = [0,0,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    for i = 1:length(placeholder_MSWD)
        if rmse_MSWD(i) < 0.0454
            color = cmap(1,:);
        elseif rmse_MSWD(i) >= 0.0454 && rmse_MSWD(i) < 0.0908
            color = cmap(2,:);
        elseif rmse_MSWD(i) >= 0.0908 && rmse_MSWD(i) < 0.1362
            color = cmap(3,:);
        elseif rmse_MSWD(i) >= 0.1362 && rmse_MSWD(i) < 0.1816
            color = cmap(4,:);
        elseif rmse_MSWD(i) >= 0.1816 && rmse_MSWD(i) < 0.227
            color = cmap(5,:);
        elseif rmse_MSWD(i) >= 0.227 && rmse_MSWD(i) < 0.2724
            color = cmap(6,:);
        elseif rmse_MSWD(i) >= 0.2724 && rmse_MSWD(i) < 0.3178
            color = cmap(7,:);
        elseif rmse_MSWD(i) >= 0.3178 && rmse_MSWD(i) < 0.3632
            color = cmap(8,:);
        elseif rmse_MSWD(i) >= 0.3632 && rmse_MSWD(i) < 0.4086
            color = cmap(9,:);
        elseif rmse_MSWD(i) >= 0.4086 && rmse_MSWD(i) < 0.454
            color = cmap(10,:);
        elseif rmse_MSWD(i) >= 0.454 && rmse_MSWD(i) < 0.5
            color = cmap(11,:);
        end
        axes(h(placeholder_MSWD(i)));
        xticklabels([])
        yticklabels([])
        if placeholder_MSWD(i) == 1 || placeholder_MSWD(i) == 7 ||...
                placeholder_MSWD(i) == 13 || placeholder_MSWD(i) == 19 || placeholder_MSWD(i) == 25
            hold on
            x = 0:0.01:1;
            y = 1-x;
            plot(x,y,'LineWidth',3.5,'Color','k')
            x = [1,1,0];
            y = [0,1,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    c=colorbar(h(25),'Position',[0.91 0.03 0.022 0.94]); caxis([0 0.5])
    c.FontSize=12;
    if ~isempty(varargin)
        exportgraphics(gcf,[varargin{2}, '_MSWD-MVMD.png'],'Resolution',600)
    end

    %% MEMD-MSWD
    cmap = [0.2422,0.1504,0.6603;...
        0.2794,0.2653,0.9094;...
        0.2621,0.4088,0.9946;...
        0.1755,0.5554,0.9473;...
        0.1009,0.6750,0.8653;...
        0.1023,0.7510,0.7088;...
        0.3176,0.7994,0.4975;...
        0.6701,0.7796,0.2239;...
        0.9523,0.7284,0.2237;...
        0.9694,0.8591,0.1665;...
        0.9769,0.9839,0.0805];
    figure
    h = tight_subplot(5, 5, [0 0],[0 0],[0 0.1]);
    placeholder_MEMD = [1,6,11,16,21,7,12,17,22,13,18,23,19,24,25];
    placeholder_MSWD = [1,2,3,4,5,7,8,9,10,13,14,15,19,20,25]; 
    for i = 1:length(placeholder_MEMD)
        if rmse_MEMD(i) < 0.0636
            color = cmap(1,:);
        elseif rmse_MEMD(i) >= 0.0636 && rmse_MEMD(i) < 0.1272
            color = cmap(2,:);
        elseif rmse_MEMD(i) >= 0.1272 && rmse_MEMD(i) < 0.1908
            color = cmap(3,:);
        elseif rmse_MEMD(i) >= 0.1908 && rmse_MEMD(i) < 0.2544
            color = cmap(4,:);
        elseif rmse_MEMD(i) >= 0.2544 && rmse_MEMD(i) < 0.3180
            color = cmap(5,:);
        elseif rmse_MEMD(i) >= 0.3180 && rmse_MEMD(i) < 0.3816
            color = cmap(6,:);
        elseif rmse_MEMD(i) >= 0.3816 && rmse_MEMD(i) < 0.4452
            color = cmap(7,:);
        elseif rmse_MEMD(i) >= 0.4452 && rmse_MEMD(i) < 0.5088
            color = cmap(8,:);
        elseif rmse_MEMD(i) >= 0.5088 && rmse_MEMD(i) < 0.5724
            color = cmap(9,:);
        elseif rmse_MEMD(i) >= 0.5724 && rmse_MEMD(i) < 0.6360
            color = cmap(10,:);
        elseif rmse_MEMD(i) >= 0.6360 && rmse_MEMD(i) < 0.7
            color = cmap(11,:);
        end
        axes(h(placeholder_MEMD(i)));
        if placeholder_MEMD(i) == 1 || placeholder_MEMD(i) == 7 ||...
                placeholder_MEMD(i) == 13 || placeholder_MEMD(i) == 19 || placeholder_MEMD(i) == 25
            x = [0,1,0];
            y = [0,0,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    for i = 1:length(placeholder_MSWD)
        if rmse_MSWD(i) < 0.0636
            color = cmap(1,:);
        elseif rmse_MSWD(i) >= 0.0636 && rmse_MSWD(i) < 0.1272
            color = cmap(2,:);
        elseif rmse_MSWD(i) >= 0.1272 && rmse_MSWD(i) < 0.1908
            color = cmap(3,:);
        elseif rmse_MSWD(i) >= 0.1908 && rmse_MSWD(i) < 0.2544
            color = cmap(4,:);
        elseif rmse_MSWD(i) >= 0.2544 && rmse_MSWD(i) < 0.3180
            color = cmap(5,:);
        elseif rmse_MSWD(i) >= 0.3180 && rmse_MSWD(i) < 0.3816
            color = cmap(6,:);
        elseif rmse_MSWD(i) >= 0.3816 && rmse_MSWD(i) < 0.4452
            color = cmap(7,:);
        elseif rmse_MSWD(i) >= 0.4452 && rmse_MSWD(i) < 0.5088
            color = cmap(8,:);
        elseif rmse_MSWD(i) >= 0.5088 && rmse_MSWD(i) < 0.5724
            color = cmap(9,:);
        elseif rmse_MSWD(i) >= 0.5724 && rmse_MSWD(i) < 0.6360
            color = cmap(10,:);
        elseif rmse_MSWD(i) >= 0.6360 && rmse_MSWD(i) < 0.7
            color = cmap(11,:);
        end
        axes(h(placeholder_MSWD(i)));
        xticklabels([])
        yticklabels([])
        if placeholder_MSWD(i) == 1 || placeholder_MSWD(i) == 7 ||...
                placeholder_MSWD(i) == 13 || placeholder_MSWD(i) == 19 || placeholder_MSWD(i) == 25
            hold on
            x = 0:0.01:1;
            y = 1-x;
            plot(x,y,'LineWidth',3.5,'Color','k')
            x = [1,1,0];
            y = [0,1,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    c=colorbar(h(25),'Position',[0.91 0.03 0.022 0.94]); caxis([0 0.5])
    c.FontSize=12;
    if ~isempty(varargin)
        exportgraphics(gcf,[varargin{2},'_MSWD-MEMD.png'],'Resolution',600)
    end
    

    %% EWT-MSWD
    cmap = [0.2422,0.1504,0.6603;...
        0.2794,0.2653,0.9094;...
        0.2621,0.4088,0.9946;...
        0.1755,0.5554,0.9473;...
        0.1009,0.6750,0.8653;...
        0.1023,0.7510,0.7088;...
        0.3176,0.7994,0.4975;...
        0.6701,0.7796,0.2239;...
        0.9523,0.7284,0.2237;...
        0.9694,0.8591,0.1665;...
        0.9769,0.9839,0.0805];
    figure
    h = tight_subplot(5, 5, [0 0],[0 0],[0 0.1]);
    placeholder_EWT = [1,6,11,16,21,7,12,17,22,13,18,23,19,24,25];
    placeholder_MSWD = [1,2,3,4,5,7,8,9,10,13,14,15,19,20,25]; 
    for i = 1:length(placeholder_EWT)
        if rmse_EWT(i) < 0.0636
            color = cmap(1,:);
        elseif rmse_EWT(i) >= 0.0636 && rmse_EWT(i) < 0.1272
            color = cmap(2,:);
        elseif rmse_EWT(i) >= 0.1272 && rmse_EWT(i) < 0.1908
            color = cmap(3,:);
        elseif rmse_EWT(i) >= 0.1908 && rmse_EWT(i) < 0.2544
            color = cmap(4,:);
        elseif rmse_EWT(i) >= 0.2544 && rmse_EWT(i) < 0.3180
            color = cmap(5,:);
        elseif rmse_EWT(i) >= 0.3180 && rmse_EWT(i) < 0.3816
            color = cmap(6,:);
        elseif rmse_EWT(i) >= 0.3816 && rmse_EWT(i) < 0.4452
            color = cmap(7,:);
        elseif rmse_EWT(i) >= 0.4452 && rmse_EWT(i) < 0.5088
            color = cmap(8,:);
        elseif rmse_EWT(i) >= 0.5088 && rmse_EWT(i) < 0.5724
            color = cmap(9,:);
        elseif rmse_EWT(i) >= 0.5724 && rmse_EWT(i) < 0.6360
            color = cmap(10,:);
        elseif rmse_EWT(i) >= 0.6360 && rmse_EWT(i) < 0.7
            color = cmap(11,:);
        end
        axes(h(placeholder_EWT(i)));
        if placeholder_EWT(i) == 1 || placeholder_EWT(i) == 7 ||...
                placeholder_EWT(i) == 13 || placeholder_EWT(i) == 19 || placeholder_EWT(i) == 25
            x = [0,1,0];
            y = [0,0,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    for i = 1:length(placeholder_MSWD)
        if rmse_MSWD(i) < 0.0636
            color = cmap(1,:);
        elseif rmse_MSWD(i) >= 0.0636 && rmse_MSWD(i) < 0.1272
            color = cmap(2,:);
        elseif rmse_MSWD(i) >= 0.1272 && rmse_MSWD(i) < 0.1908
            color = cmap(3,:);
        elseif rmse_MSWD(i) >= 0.1908 && rmse_MSWD(i) < 0.2544
            color = cmap(4,:);
        elseif rmse_MSWD(i) >= 0.2544 && rmse_MSWD(i) < 0.3180
            color = cmap(5,:);
        elseif rmse_MSWD(i) >= 0.3180 && rmse_MSWD(i) < 0.3816
            color = cmap(6,:);
        elseif rmse_MSWD(i) >= 0.3816 && rmse_MSWD(i) < 0.4452
            color = cmap(7,:);
        elseif rmse_MSWD(i) >= 0.4452 && rmse_MSWD(i) < 0.5088
            color = cmap(8,:);
        elseif rmse_MSWD(i) >= 0.5088 && rmse_MSWD(i) < 0.5724
            color = cmap(9,:);
        elseif rmse_MSWD(i) >= 0.5724 && rmse_MSWD(i) < 0.6360
            color = cmap(10,:);
        elseif rmse_MSWD(i) >= 0.6360 && rmse_MSWD(i) < 0.7
            color = cmap(11,:);
        end
        axes(h(placeholder_MSWD(i)));
        xticklabels([])
        yticklabels([])
        if placeholder_MSWD(i) == 1 || placeholder_MSWD(i) == 7 ||...
                placeholder_MSWD(i) == 13 || placeholder_MSWD(i) == 19 || placeholder_MSWD(i) == 25
            hold on
            x = 0:0.01:1;
            y = 1-x;
            plot(x,y,'LineWidth',3.5,'Color','k')
            x = [1,1,0];
            y = [0,1,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    c=colorbar(h(25),'Position',[0.91 0.03 0.022 0.94]); caxis([0 0.5])
    c.FontSize=12;
    if ~isempty(varargin)
        exportgraphics(gcf,[varargin{2},'_MSWD-EWT.png'],'Resolution',600)
    end

    %% IRCNN-MSWD
    cmap = [0.2422,0.1504,0.6603;...
        0.2794,0.2653,0.9094;...
        0.2621,0.4088,0.9946;...
        0.1755,0.5554,0.9473;...
        0.1009,0.6750,0.8653;...
        0.1023,0.7510,0.7088;...
        0.3176,0.7994,0.4975;...
        0.6701,0.7796,0.2239;...
        0.9523,0.7284,0.2237;...
        0.9694,0.8591,0.1665;...
        0.9769,0.9839,0.0805];
    figure
    h = tight_subplot(5, 5, [0 0],[0 0],[0 0.1]);
    placeholder_IRCNN = [1,6,11,16,21,7,12,17,22,13,18,23,19,24,25];
    placeholder_MSWD = [1,2,3,4,5,7,8,9,10,13,14,15,19,20,25]; 
    for i = 1:length(placeholder_IRCNN)
        if rmse_IRCNN(i) < 0.1018
            color = cmap(1,:);
        elseif rmse_IRCNN(i) >= 0.1018 && rmse_IRCNN(i) < 0.2036
            color = cmap(2,:);
        elseif rmse_IRCNN(i) >= 0.2036 && rmse_IRCNN(i) < 0.3054
            color = cmap(3,:);
        elseif rmse_IRCNN(i) >= 0.3054 && rmse_IRCNN(i) < 0.4072
            color = cmap(4,:);
        elseif rmse_IRCNN(i) >= 0.4072 && rmse_IRCNN(i) < 0.5090
            color = cmap(5,:);
        elseif rmse_IRCNN(i) >= 0.5090 && rmse_IRCNN(i) < 0.6108
            color = cmap(6,:);
        elseif rmse_IRCNN(i) >= 0.6108 && rmse_IRCNN(i) < 0.7126
            color = cmap(7,:);
        elseif rmse_IRCNN(i) >= 0.7126 && rmse_IRCNN(i) < 0.8144
            color = cmap(8,:);
        elseif rmse_IRCNN(i) >= 0.8144 && rmse_IRCNN(i) < 0.9162
            color = cmap(9,:);
        elseif rmse_IRCNN(i) >= 0.9162 && rmse_IRCNN(i) < 1.0180
            color = cmap(10,:);
        elseif rmse_IRCNN(i) >= 1.0180 && rmse_IRCNN(i) < 1.12
            color = cmap(11,:);
        end
        axes(h(placeholder_EWT(i)));
        if placeholder_EWT(i) == 1 || placeholder_EWT(i) == 7 ||...
                placeholder_EWT(i) == 13 || placeholder_EWT(i) == 19 || placeholder_EWT(i) == 25
            x = [0,1,0];
            y = [0,0,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    for i = 1:length(placeholder_MSWD)
        if rmse_MSWD(i) < 0.1018
            color = cmap(1,:);
        elseif rmse_MSWD(i) >= 0.1018 && rmse_MSWD(i) < 0.2036
            color = cmap(2,:);
        elseif rmse_MSWD(i) >= 0.2036 && rmse_MSWD(i) < 0.3054
            color = cmap(3,:);
        elseif rmse_MSWD(i) >= 0.3054 && rmse_MSWD(i) < 0.4072
            color = cmap(4,:);
        elseif rmse_MSWD(i) >= 0.4072 && rmse_MSWD(i) < 0.5090
            color = cmap(5,:);
        elseif rmse_MSWD(i) >= 0.5090 && rmse_MSWD(i) < 0.6108
            color = cmap(6,:);
        elseif rmse_MSWD(i) >= 0.6108 && rmse_MSWD(i) < 0.7126
            color = cmap(7,:);
        elseif rmse_MSWD(i) >= 0.7126 && rmse_MSWD(i) < 0.8144
            color = cmap(8,:);
        elseif rmse_MSWD(i) >= 0.8144 && rmse_MSWD(i) < 0.9162
            color = cmap(9,:);
        elseif rmse_MSWD(i) >= 0.9162 && rmse_MSWD(i) < 1.0180
            color = cmap(10,:);
        elseif rmse_MSWD(i) >= 1.0180 && rmse_MSWD(i) < 1.12
            color = cmap(11,:);
        end
        axes(h(placeholder_MSWD(i)));
        xticklabels([])
        yticklabels([])
        if placeholder_MSWD(i) == 1 || placeholder_MSWD(i) == 7 ||...
                placeholder_MSWD(i) == 13 || placeholder_MSWD(i) == 19 || placeholder_MSWD(i) == 25
            hold on
            x = 0:0.01:1;
            y = 1-x;
            plot(x,y,'LineWidth',3.5,'Color','k')
            x = [1,1,0];
            y = [0,1,1];
        else
            x = [0,1,1,0];
            y = [0,0,1,1];
        end
        fill(x,y,color,'LineStyle','none');
        xticklabels([])
        yticklabels([])
        xticks([])
        yticks([])
    end
    c=colorbar(h(25),'Position',[0.91 0.03 0.022 0.94]); caxis([0 0.5])
    c.FontSize=12;
    if ~isempty(varargin)
        exportgraphics(gcf,[varargin{2},'_MSWD-IRCNN.png'],'Resolution',600)
    end
end

end



