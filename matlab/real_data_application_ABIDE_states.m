% clear
clc
close all

signal_length = 200;
numICs = 100;
indx = nchoosek(1:numICs,2);
save_res = true;
plotStateTransitions = false;

%% Load TVPS data and concatenate them across subjects
load(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\COSDELPHI_MSWD_20_subjects.mat'])
load(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\COSDELPHI_MVMD_20_subjects.mat'])

for i = 1:length(COSDELPHI1_MSWD)
    if i == 1
        COSDELPHI_MSWD = COSDELPHI1_MSWD{i};
        COSDELPHI_MVMD = COSDELPHI1_MVMD{i};
    else
        COSDELPHI_MSWD = cat(1,COSDELPHI_MSWD,COSDELPHI1_MSWD{i});
        COSDELPHI_MVMD = cat(1,COSDELPHI_MVMD,COSDELPHI1_MVMD{i});
    end
end

%% Find optimal number of clusters
DBI_MVMD = 2;
DBI_MSWD = 2;

%% Perform clustering and get the centroids and the clustered data
if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\mean_C_MVMD_all.mat']) ||...
    ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\mean_C_MSWD_all.mat'])
    [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
    [mean_idx_MVMD_all{1},mean_C_MVMD_all{1}] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
    
    [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
    [mean_idx_MSWD_all{1},mean_C_MSWD_all{1}] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
    save(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\mean_C_MVMD_all.mat'],'mean_idx_MVMD_all','mean_C_MVMD_all');
    save(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\mean_C_MSWD_all.mat'],'mean_idx_MSWD_all',"mean_C_MSWD_all");
else
    load(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\mean_C_MVMD_all.mat']);
    load(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\mean_C_MSWD_all.mat']);
end


%% Count short state transitions
count = 1;
for i = 1:signal_length:length(mean_idx_MVMD_all{1})
    diff_states = diff(mean_idx_MVMD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MVMD(count) = sum((diff_poss>0).*(diff_poss<=10));
    count = count + 1;
end

count = 1;
for i = 1:signal_length:length(mean_idx_MSWD_all{1})
    diff_states = diff(mean_idx_MSWD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MSWD(count) = sum((diff_poss>0).*(diff_poss<=10));
    count = count + 1;
end

disp(['MVMD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MVMD))])
disp(['MSwD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MSWD))])

%% Plot the centroids/brain states for the two methods
figure;
h = tight_subplot(1, DBI_MVMD, [0.001 0.02],[.001 .001],[.05 .1]);
for i = 1:size(mean_C_MVMD_all{1},1)
    state = mean_C_MVMD_all{1}(i,:);
    state_new = eye(numICs,numICs);
    for n = 1:size(indx,1)
        state_new(indx(n,1),indx(n,2)) = state(n);
        state_new(indx(n,2),indx(n,1)) = state(n);
    end
    state_new = state_new(1:end-1,1:end-1);
    for j = 1:size(state_new,1)
        state_new(j,j) = 0;
    end
    axes(h(i));gsplot(state_new);
    set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
    c = get(gca, 'Children');
    axis square;
    if i == DBI_MVMD
        colorbar(h(i),'Position',[0.93 0.225 0.022 0.548]); caxis([-1 1])
    end
end
if save_res
    exportgraphics(gcf,fullfile(pwd,'plots_real_application',['ABIDE','_states_zero_diag_MVMD.png']),'Resolution',600)
end

figure;
h = tight_subplot(1, DBI_MSWD, [0.001 0.02],[.001 .001],[.05 .1]);
for i = 1:size(mean_C_MSWD_all{1},1)
    state = mean_C_MSWD_all{1}(i,:);
    state_new = eye(numICs,numICs);
    for n = 1:size(indx,1)
        state_new(indx(n,1),indx(n,2)) = state(n);
        state_new(indx(n,2),indx(n,1)) = state(n);
    end
    state_new = state_new(1:end-1,1:end-1);
    for j = 1:size(state_new,1)
        state_new(j,j) = 0;
    end
    axes(h(i));gsplot(state_new);
    set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])
    c = get(gca, 'Children');
    axis square;
    if i == DBI_MSWD
        colorbar(h(i),'Position',[0.93 0.225 0.022 0.548]); caxis([-1 1])
    end
end
if save_res
    exportgraphics(gcf,fullfile(pwd,'plots_real_application',['ABIDE','_states_zero_diag_MSWD.png']),'Resolution',600)
end