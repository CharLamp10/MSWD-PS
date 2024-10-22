function simulation4Plots(mCorr_MVMD,midx_MVMD,mCorr_MSwD,midx_MSwD,grdidx,nS,N,varargin)

for i = 1:N
    for j = 1:nS 
        crpstate_MVMD{j}(:,:,i) = mCorr_MVMD{i}(:,:,j);
        crpstate_MSwD{j}(:,:,i) = mCorr_MSwD{i}(:,:,j);
    end
end

for i = 1:nS
    meanstate_MVMD{i} = mean(crpstate_MVMD{i},3,'omitnan');
    meanstate_MSwD{i} = mean(crpstate_MSwD{i},3,'omitnan');
end

placeholder = 1:2*nS;%[1 2 3;7 8 9];
figure;
h = tight_subplot(2, nS, [0.001 0.02],[.001 .001],[.05 .1]);
%set(h(4:6),'visible','off')
count = 1;
for i = 1:2
    for j = 1:nS
        if i == 1
            axes(h(placeholder(count)));gsplot(meanstate_MVMD{j});
        else
            axes(h(placeholder(count)));gsplot(meanstate_MSwD{j});
        end
        count = count + 1;
        axis square;
        set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
        c = get(gca, 'Children');
        set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w')
        title(strcat(['State ' num2str(j)]));
        % if j == 1
        %     ylabel(method,'interpreter','latex');
        % end
    end
end
colorbar(h(nS),'Position',[0.91 0.07 0.022 0.86],'FontSize',12); caxis([-1 1])
if ~isempty(varargin)
    exportgraphics(gcf,[varargin{1},'mean_states.png'],'Resolution',600)
end


for i = 1:N
    MIDX_MVMD(:,i) = midx_MVMD{i};
    MIDX_MSwD(:,i) = midx_MSwD{i};
end

for i = 1:size(MIDX_MVMD,1)
    flag_MVMD(i) = sum(MIDX_MVMD(i,:) == grdidx(i))/N;
    flag_MSwD(i) = sum(MIDX_MSwD(i,:) == grdidx(i))/N;
end

for j = 1:nS
    for i = 1:size(MIDX_MVMD,1)
        flags_MVMD(i,j) = sum(MIDX_MVMD(i,:) == j)/N;
        flags_MSwD(i,j) = sum(MIDX_MSwD(i,:) == j)/N;
    end
end


% figure;subplot(3,1,1);plot(grdidx,'r','LineWidth',2);ylim([0 4]);xlabel('t');
% ylabel('State #');title('Ground truth state')
% subplot(3,1,2);plot(flag,'g','LineWidth',2);
% % title([method,' :Accuracy of correctly classifying the state']);ylabel('Accuracy');xlabel('t')
% subplot(3,1,3);bar(flags,'stacked')
% title('Stacked bar of the classification accuracy across time');
% legend('State 1','State 2','State 3');
% ylabel('State Accuracy Proportion');xlabel('t')
% if ~isempty(varargin)
%     exportgraphics(gcf,[varargin{1},'state_accuracy.png'],'Resolution',600)
% end

customColormap = [1 1 0;    % Yellow for value 1
                  0 1 1;    % Cyan for value 2
                  1 0 1];   % Magenta for value 3

figure;
imagesc(MIDX_MVMD');colormap(customColormap);
alphaVal = 0.4;
imAlpha=alphaVal*ones(size(MIDX_MVMD'));
imAlpha(isnan(MIDX_MVMD'))=0;
imagesc(MIDX_MVMD','AlphaData',imAlpha);
% set(gca,'color',1*[1 1 1]);
% title(['State transitions - ',method,'-based-PS using tWPS: CIRC']);
xt = get(gca, 'XTick');
xtlbl = linspace(100, 1000, numel(xt));
set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',11)
xlabel('Time (sec)','FontSize',15);ylabel('Realizations','FontSize',15)
% cb = colorbar();
% drawnow
% % Set color labels (one for each row in RGB)
% label = 1:nS; 
% caxis([1,numel(label)])
% cb.YTick = 1 : nS;
% labelChar = label;
% cb.TickLabels = labelChar(1:end-1);
% cb.FontSize = 12; 
% cdata = cb.Face.Texture.CData;
% cdata(end,:) = uint8(alphaVal * cdata(end,:));
% cb.Face.Texture.ColorType = 'truecoloralpha';
% cb.Face.Texture.CData = cdata;
% drawnow
% cb.Face.ColorBinding = 'discrete';
if ~isempty(varargin)
    exportgraphics(gcf, [varargin{1},'state_transitions_MVMD.png'],"Resolution",600);
end

figure
imagesc(MIDX_MSwD');colormap(customColormap);
alphaVal = 0.4;
imAlpha=alphaVal*ones(size(MIDX_MSwD'));
imAlpha(isnan(MIDX_MSwD'))=0;
imagesc(MIDX_MSwD','AlphaData',imAlpha);
xt = get(gca, 'XTick');
xtlbl = linspace(100, 1000, numel(xt));
set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',11)
% set(gca,'color',1*[1 1 1]);
% title(['State transitions - ',method,'-based-PS using tWPS: CIRC']);
xlabel('Time (sec)','FontSize',15);ylabel('Realizations','FontSize',15)
% cb = colorbar();
% drawnow
% Set color labels (one for each row in RGB)
% label = 1:nS; 
% caxis([1,numel(label)])
% cb.YTick = 1 : nS;
% labelChar = label;
% cb.TickLabels = labelChar(1:end-1);
% cb.FontSize = 12; 
% cb.FontSize = 12; 
% cdata = cb.Face.Texture.CData;
% cdata(end,:) = uint8(alphaVal * cdata(end,:));
% cb.Face.Texture.ColorType = 'truecoloralpha';
% cb.Face.Texture.CData = cdata;
% drawnow
% cb.Face.ColorBinding = 'discrete';
if ~isempty(varargin)
    exportgraphics(gcf, [varargin{1},'state_transitions_MSwD.png'],"Resolution",600);
end

end
