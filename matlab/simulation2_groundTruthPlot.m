function simulation2_groundTruthPlot(cosdelphi,t,varargin)

figure
plot(t,cosdelphi,"Color",'b','LineWidth',2)
grid on

ylabel('CRP','FontSize',12)
xlabel('Time (sec)','FontSize',12)
xlim([0,max(t)])
ylim([-1,1])
hold off
sgtitle('Ground Truth')

if ~isempty(varargin)
    saveas(gcf,varargin{1})
end

end