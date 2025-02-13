function simulation12356Plots(cosdelphi_MEMD,cosdelphi_EWT,cosdelphi_IRCNN,cosdelphi_MVMD,cosdelphi_MSWD,t,varargin)

[cosdelphi_MEMD_mean,cosdelphi_MEMD_low,cosdelphi_MEMD_up] = confidence_interval(cosdelphi_MEMD);
[cosdelphi_EWT_mean,cosdelphi_EWT_low,cosdelphi_EWT_up] = confidence_interval(cosdelphi_EWT);
[cosdelphi_IRCNN_mean,cosdelphi_IRCNN_low,cosdelphi_IRCNN_up] = confidence_interval(cosdelphi_IRCNN);
[cosdelphi_MVMD_mean,cosdelphi_MVMD_low,cosdelphi_MVMD_up] = confidence_interval(cosdelphi_MVMD);
[cosdelphi_MSWD_mean,cosdelphi_MSWD_low,cosdelphi_MSWD_up] = confidence_interval(cosdelphi_MSWD);


figure
h = tight_subplot(5, 1, [0.06 0.02],[0.12 0.1],[0.05 0.01]);

axes(h(1))
plot(t,cosdelphi_MEMD_mean,"Color",'m','LineWidth',2)
hold on
plot(t,cosdelphi_MEMD_low,"Color",'m')
plot(t,cosdelphi_MEMD_up,"Color",'m')
grid on

patch([t fliplr(t)],[cosdelphi_MEMD_low fliplr(cosdelphi_MEMD_up)],'m', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
title('MEMD-PS','FontSize',12)
xlim([0,max(t)])
xticklabels([])
ylim([-1,1])
hold off

axes(h(2))
plot(t,cosdelphi_EWT_mean,"Color",'c','LineWidth',2)
hold on
plot(t,cosdelphi_EWT_low,"Color",'c')
plot(t,cosdelphi_EWT_up,"Color",'c')
grid on

patch([t fliplr(t)],[cosdelphi_EWT_low fliplr(cosdelphi_EWT_up)],'c', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
title('ETW-PS','FontSize',12)
xlim([0,max(t)])
xticklabels([])
ylim([-1,1])
hold off

axes(h(3))
plot(t,cosdelphi_IRCNN_mean,"Color",'y','LineWidth',2)
hold on
plot(t,cosdelphi_IRCNN_low,"Color",'y')
plot(t,cosdelphi_IRCNN_up,"Color",'y')
grid on

patch([t fliplr(t)],[cosdelphi_IRCNN_low fliplr(cosdelphi_IRCNN_up)],'y', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
title('IRCNN-PS','FontSize',12)
xlim([0,max(t)])
xticklabels([])
ylim([-1,1])
hold off

axes(h(4))
plot(t,cosdelphi_MVMD_mean,"Color",'b','LineWidth',2)
hold on
plot(t,cosdelphi_MVMD_low,"Color",'b')
plot(t,cosdelphi_MVMD_up,"Color",'b')
grid on

patch([t fliplr(t)],[cosdelphi_MVMD_low fliplr(cosdelphi_MVMD_up)],'b', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
title('MVMD-PS','FontSize',12)
xlim([0,max(t)])
xticklabels([])
ylim([-1,1])
hold off

axes(h(5))
plot(t,cosdelphi_MSWD_mean,"Color",'g','LineWidth',2)
hold on
plot(t,cosdelphi_MSWD_low,"Color",'g')
plot(t,cosdelphi_MSWD_up,"Color",'g')
grid on

patch([t fliplr(t)],[cosdelphi_MSWD_low fliplr(cosdelphi_MSWD_up)],'g', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
title('MSwD-PS','FontSize',12)
xlabel('Time (sec)','FontSize',16)
xlim([0,max(t)])
ylim([-1,1])
hold off

if ~isempty(varargin)
    exportgraphics(gcf,varargin{1},'Resolution',1000)
end
end