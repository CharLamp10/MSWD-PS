function simulation12356Plots(true_cosdelphi,cosdelphi_MVMD,cosdelphi_MSWD,t,varargin)

[cosdelphi_MVMD_mean,cosdelphi_MVMD_low,cosdelphi_MVMD_up] = confidence_interval(cosdelphi_MVMD);
[cosdelphi_MSWD_mean,cosdelphi_MSWD_low,cosdelphi_MSWD_up] = confidence_interval(cosdelphi_MSWD);


figure
h = tight_subplot(3, 1, [0.04 0.02],[0.12 0.1],[0.1 0.1]);

axes(h(1))
plot(t,true_cosdelphi,"Color",'r','LineWidth',1.5)
grid on
ax = gca;
ax.FontSize = 11; 
ylabel('True','FontSize',16)
xlim([0,max(t)])
xticklabels([])
ylim([-1,1])

axes(h(2))
plot(t,cosdelphi_MVMD_mean,"Color",'b','LineWidth',2)
hold on
plot(t,cosdelphi_MVMD_low,"Color",'b')
plot(t,cosdelphi_MVMD_up,"Color",'b')
grid on

patch([t fliplr(t)],[cosdelphi_MVMD_low fliplr(cosdelphi_MVMD_up)],'b', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
ylabel('MVMD-PS','FontSize',16)
xlim([0,max(t)])
xticklabels([])
ylim([-1,1])
hold off

axes(h(3))
plot(t,cosdelphi_MSWD_mean,"Color",'g','LineWidth',2)
hold on
plot(t,cosdelphi_MSWD_low,"Color",'g')
plot(t,cosdelphi_MSWD_up,"Color",'g')
grid on

patch([t fliplr(t)],[cosdelphi_MSWD_low fliplr(cosdelphi_MSWD_up)],'g', 'FaceAlpha',0.5, 'EdgeColor','none')
ax = gca;
ax.FontSize = 11; 
ylabel('MSwD-PS','FontSize',16)
xlabel('Time (sec)','FontSize',16)
xlim([0,max(t)])
ylim([-1,1])
hold off

if ~isempty(varargin)
    exportgraphics(gcf,varargin{1},'Resolution',1000)
end
end