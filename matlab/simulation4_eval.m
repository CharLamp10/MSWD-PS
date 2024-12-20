clear
clc
close all


%% Sim 4
save_plots = 1;
casee = 4;
TR = 2;                                                                     % Repetition time
fs = 1/TR;                                                                  % Sampling frequency
t = 0:1/fs:1000-1/fs;
f = 0.05;
N = 100;                                                 % number of repetition (realizations)
nS = 3;
freq_ranges = 0.01:0.02:0.1;
min_freq = 0;                                        % window sizes for the Windowed Phase Sync. Measures

true_plv1 = zeros(1,length(freq_ranges));
true_plv2 = zeros(1,length(freq_ranges));
true_plv3 = zeros(1,length(freq_ranges));

phi1 = pi*((t-100>0).*(t-250<0) + (t-300>0).*(t-500<0) + (t-600>0).*(t-800<0));
phi2 = pi*((t-100>0).*(t-250<0) +                      - (t-600>0).*(t-800<0));
phi3 = pi*(                     - (t-300>0).*(t-500<0) - (t-600>0).*(t-800<0));
x = cos(2*pi*f*t + phi1);    % first signal
y = cos(2*pi*f*t + phi2);    % second signal
z = cos(2*pi*f*t + phi3);    % thrid signal

pwelch_window = length(x);

phase_x = angle(hilbert(x)); 
phase_y = angle(hilbert(y));
phase_z = angle(hilbert(z));

plv1 = abs(mean(exp(1i*(phase_x - phase_y))));
plv2 = abs(mean(exp(1i*(phase_x - phase_z))));
plv3 = abs(mean(exp(1i*(phase_y - phase_z))));
pos = find(freq_ranges == f);
true_plv1(pos) = plv1; true_plv2(pos) = plv2; true_plv3(pos) = plv3;
true_plv = {true_plv1,true_plv2,true_plv3};

figure;subplot(3,1,1);plot(t,phi1, 'k', 'LineWidth', 1.5); ylim([-1.3*pi 1.3*pi]); xline(100,'-- k'); xline(250,'-- k'); xline(300,'-- k'); xline(500,'-- k'); xline(500,'-- k'); xline(600,'-- k'); xline(800,'-- k');
subplot(3,1,2);plot(t,phi2, 'r', 'LineWidth', 1.5);ylim([-1.3*pi 1.3*pi]); xline(100,'-- k'); xline(250,'-- k'); xline(300,'-- k'); xline(500,'-- k'); xline(500,'-- k');xline(600,'-- k'); xline(800,'-- k');
subplot(3,1,3);plot(t,phi3, 'b', 'LineWidth', 1.5);ylim([-1.3*pi 1.3*pi]); xline(100,'-- k'); xline(250,'-- k'); xline(300,'-- k'); xline(500,'-- k'); xline(500,'-- k');xline(600,'-- k'); xline(800,'-- k');
xlabel('Time (sec)','FontSize',12)
subplot(3,1,1);%title('Ground truth phases of the signals to generate multivariate signal')
exportgraphics(gcf, fullfile(pwd,'plots','Ground_truth_phases.png'),'Resolution',600)

figure; subplot(3,1,1);plot(t,cos(phi1-phi2), 'k', 'LineWidth', 1.5); hold on; ylim([-2 2]); %xline(100,'-- k'); xline(250,'-- k'); xline(300,'-- k'); xline(500,'-- k'); xline(500,'-- k');
x = [0 100 100 0]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none');
x = [100 250 250 100]; y = [-2 -2 2 2]; fill(x,y,'magenta','FaceAlpha',0.3,'LineStyle','none');
x = [250 300 300 250]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none');
x = [300 500 500 300]; y = [-2 -2 2 2]; fill(x,y,'yellow','FaceAlpha',0.3,'LineStyle','none');
x = [500 1000 1000 500]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none'); hold off
subplot(3,1,2);plot(t,cos(phi1-phi3), 'k', 'LineWidth', 1.5); hold on; ylim([-2 2]); %xline(100,'-- k'); xline(250,'-- k'); xline(300,'-- k'); xline(500,'-- k'); xline(500,'-- k');
x = [0 100 100 0]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none');
x = [100 250 250 100]; y = [-2 -2 2 2]; fill(x,y,'magenta','FaceAlpha',0.3,'LineStyle','none');
x = [250 300 300 250]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none');
x = [300 500 500 300]; y = [-2 -2 2 2]; fill(x,y,'yellow','FaceAlpha',0.3,'LineStyle','none');
x = [500 1000 1000 500]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none'); hold off
subplot(3,1,3);plot(t,cos(phi2-phi3), 'k', 'LineWidth', 1.5); hold on; ylim([-2 2]); %xline(100,'-- k'); xline(250,'-- k'); xline(300,'-- k'); xline(500,'-- k'); xline(500,'-- k');
x = [0 100 100 0]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none');
x = [100 250 250 100]; y = [-2 -2 2 2]; fill(x,y,'magenta','FaceAlpha',0.3,'LineStyle','none');
x = [250 300 300 250]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none');
x = [300 500 500 300]; y = [-2 -2 2 2]; fill(x,y,'yellow','FaceAlpha',0.3,'LineStyle','none');
x = [500 1000 1000 500]; y = [-2 -2 2 2]; fill(x,y,'cyan','FaceAlpha',0.3,'LineStyle','none'); hold off
xlabel('Time (sec)','FontSize',12)
subplot(3,1,1);%title('Ground truth phases of the signals to generate multivariate signal')
if save_plots
    exportgraphics(gcf, fullfile(pwd,'plots','Ground_truth_cosine_phase_differences.png'),'Resolution',600)
end

grdidx = zeros(size(t));
grdidx(1:50) = 2;
grdidx(125:150) = 2;
grdidx(250:500) = 2;
grdidx(50:125) = 3;
grdidx(150:250) = 1;
grdstate{1}(:,:,1) = 0.98*[1 -1 1;-1 1 -1;1 -1 1];
grdstate{1}(:,:,2) = 0.98*[1 1 1;1 1 1;1 1 1];
grdstate{1}(:,:,3) = 0.98*[1 1 -1;1 1 -1;-1 -1 1];

figure;
placeholder = 1:nS;
h = tight_subplot(1, nS, [0.001 0.02],[.001 .001],[.05 .1]);
labels = {'Ground Truth States'};
for j = 1:nS
    axes(h(placeholder(j)));gsplot(grdstate{1}(:,:,j));
    axis square;
    set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
    c = get(gca, 'Children');
    set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w')
    if j == 1
        ylabel(labels{j},'interpreter','latex');
    end
end
colorbar(h(3),'Position',[0.93 0.3 0.022 0.4],'FontSize',10); caxis([-1 1])
if save_plots
    saveas(gcf, fullfile(pwd,'plots','Ground_truth_states.png'))
end

indx = nchoosek(1:3,2);

load(fullfile("decomposed","sim4_MVMD.mat"))
load(fullfile("decomposed","sim4_MSWD.mat"))

%% MVMD
for m = 1:N
    imf = imfs_MVMD{m};
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis47(imf,"MVMD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI1_MVMD{m} = COSDELPHI11;
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end

for i = 1:size(indx,1)
    [plv_MVMD(:,:,i),counts_MVMD(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
    [mean_plv_MVMD(:,i),lower_plv_MVMD(:,i),upper_plv_MVMD(:,i)] = confidence_interval(plv_MVMD(:,:,i));
end

for i = 1:N
    % k-means clustering of the matrices of phase synch.
    [idx{i},Corr{i}] = mykmeans(COSDELPHI1_MVMD{i},nS,nS);
    [mCorr_MVMD{i},midx_MVMD{i}] = matchstatesClutster(grdstate{1},Corr{i},idx{i},nS);
end

disp('MVMD')
clear mfreq plv idx Corr midx mCorr


%% MSWD
for m = 1:N
    imf = imfs_MSWD{m};
    if m == 1
        plv_length = [length(freq_ranges),size(indx,1)];
    else
        plv_length = [length(plv(m-1,:,1)),size(indx,1)];
    end
    [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis47(imf,"MSWD",indx,...
     fs,f,min_freq,plv_length);
    COSDELPHI1_MSWD{m} = COSDELPHI11;
    plv(m,1:size(plv1,1),1:size(plv1,2)) = plv1;
    mfreq{m} = mfreq1;
end
for i = 1:size(indx,1)
    [plv_MSWD(:,:,i),counts_MSWD(:,:,i)] = ...
        sort_values(mfreq,freq_ranges,plv(:,:,i));
    [mean_plv_MSWD(:,i),lower_plv_MSWD(:,i),upper_plv_MSWD(:,i)] = confidence_interval(plv_MSWD(:,:,i));
end

for i = 1:N
    % k-means clustering of the matrices of phase synch.
    % if ~all(isnan(COSDELPHI1_MSWD{i}))
    [idx{i},Corr{i}] = mykmeans(COSDELPHI1_MSWD{i},nS,nS);
    [mCorr_MSwD{i},midx_MSwD{i}] = matchstatesClutster(grdstate{1},Corr{i},idx{i},nS);
    % else
    %     mCorr{1,i} = nan(nS,nS,size(indx,1)); mCorr{2,i} = nan(nS,nS,size(indx,1));
    %     midx{1,i} = nan(length(x),1); midx{2,i} = nan(length(x),1);
    % end
end
disp('MSWD')

%% PLOTS
filename = fullfile(pwd,'plots',['sim', num2str(casee), '_']);
if save_plots
    simulation4Plots(mCorr_MVMD,midx_MVMD,mCorr_MSwD,midx_MSwD,grdidx,nS,N,filename)
else
    simulation4Plots(mCorr_MVMD,midx_MVMD,mCorr_MSwD,midx_MSwD,grdidx,nS,N)
end



for i = 1:size(indx,1)
    data_plv(:,i) = [mean_plv_MVMD(:,i)',mean_plv_MSWD(:,i)'];

    E_plv_MVMD = mean_plv_MVMD(:,i) - lower_plv_MVMD(:,i);
    E_plv_MSWD = mean_plv_MSWD(:,i) - lower_plv_MSWD(:,i);
    
    
    plv_min = min([min(lower_plv_MVMD(:,i)),...
        min(lower_plv_MSWD(:,i)),min(true_plv{i})]);
    plv_max = max([max(upper_plv_MVMD(:,i)),...
        max(upper_plv_MSWD(:,i)),max(true_plv{i})]);


    figure
    errorbar(1:size(plv_MVMD,2),mean_plv_MVMD(:,i),E_plv_MVMD,'bo',"LineStyle","none",'LineWidth',1)
    hold on
    errorbar(size(plv_MVMD,2)+1:size(plv_MVMD,2)+size(plv_MSWD,2),mean_plv_MSWD(:,i),E_plv_MSWD,'go',"LineStyle","none",'LineWidth',1)
    plot([true_plv{i},true_plv{i}],'ro','LineWidth',1)
%     title(['Simulation 5: PLV (combination ', num2str(i),')'])
%     xlabel('Method')
    ylabel('PLV','FontSize',14)
    xlim([1,size(plv_MVMD,2)+size(plv_MSWD,2)])
    ylim([0,plv_max+0.02])
    start = (1 + length(freq_ranges))/2;
    stop = 2*length(freq_ranges) + 1 - start;
    ticks = start:length(freq_ranges):stop;
    xticks(ticks)
    xticklabels({'MVMD-PS','MSwD-PS'})
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',14)
    xline((length(freq_ranges):length(freq_ranges):length(freq_ranges))+0.5,'--')
    hold off
    if save_plots
        exportgraphics(gcf,fullfile(pwd,'plots',['sim',num2str(casee),'_comb',num2str(i),'_plv.png']),'Resolution',1000)
    end
end




