% In this script, the results of the time-varying phase synchronization
% analysis on rs-fMRI data from the HCP are extracted

clear
clc
close all

path_data = 'F:\MSWD_paper_new';

numSubjects = 50;
signal_length = 1200;
visit = 'rest1';
numICs = 100;
indx = nchoosek(1:numICs,2);
save_res = true;
plotStateTransitions = false;

%% Load TVPS data and concatenate them across subjects
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_50_subjects.mat']))
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_50_subjects.mat']))

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
if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']))
    E_MVMD = evalclusters(squeeze(COSDELPHI_MVMD),'kmeans','silhouette','klist',[2:5]);
    DBI_MVMD = E_MVMD.OptimalK;
end
if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']))
    E_MSWD = evalclusters(squeeze(COSDELPHI_MSWD),'kmeans','silhouette','klist',[2:5]);
    DBI_MSWD = E_MSWD.OptimalK;
end

%% Perform clustering and get the centroids and the clustered data
if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']))
    [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
    [mean_idx_MVMD_all{1},mean_C_MVMD_all{1}] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
    
    [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
    [mean_idx_MSWD_all{1},mean_C_MSWD_all{1}] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
    save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']),'mean_idx_MVMD_all','mean_C_MVMD_all');
    save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']),'mean_idx_MSWD_all',"mean_C_MSWD_all");
else
    load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']));
    load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']));
end


%% Count short state transitions
count = 1;
for i = 1:signal_length:length(mean_idx_MVMD_all{1})
    diff_states = diff(mean_idx_MVMD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MVMD(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
end

count = 1;
for i = 1:signal_length:length(mean_idx_MSWD_all{1})
    diff_states = diff(mean_idx_MSWD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MSWD(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
end
disp(['MVMD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MVMD))])
disp(['MSwD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MSWD))])


%% Bootstrap analysis to show that the centroids/brain states are robust
if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']))
    for i = 1:25
        N = randperm(numSubjects,10);
        mask = zeros(size(COSDELPHI_MSWD,1),1);
        for n = N
            start = (n-1)*signal_length+1;
            endd = signal_length*n;
            mask(start:endd) = 1;
        end
        COSDELPHI_MVMD_btsrp = COSDELPHI_MVMD(logical(mask),:);
        COSDELPHI_MSWD_btsrp = COSDELPHI_MSWD(logical(mask),:);
        [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD_btsrp),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_MVMD_bstrp{i},mean_C_MVMD_bstrp{i}] = kmeans(squeeze(COSDELPHI_MVMD_btsrp),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
        
        [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD_btsrp),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_MSWD_bstrp{i},mean_C_MSWD_bstrp{i}] = kmeans(squeeze(COSDELPHI_MSWD_btsrp),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
        disp(i)
    end
    save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat']),"mean_C_MVMD_bstrp")
    save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']),"mean_C_MSWD_bstrp")
else
    load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat']))
    load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']))
end

%% Correlation coefficient between the bootstrap centroids and the original centroids
for i = 1:size(mean_C_MSWD_bstrp,1)
    r1 = corrcoef(mean_C_MVMD_bstrp{i}(1,:),mean_C_MVMD_all{1}(1,:));
    r2 = corrcoef(mean_C_MVMD_bstrp{i}(2,:),mean_C_MVMD_all{1}(2,:));
    r3 = corrcoef(mean_C_MVMD_bstrp{i}(1,:),mean_C_MVMD_all{1}(2,:));
    r4 = corrcoef(mean_C_MVMD_bstrp{i}(2,:),mean_C_MVMD_all{1}(1,:));
    if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
        corrs_MVMD(i,1) = r1(1,2);
        corrs_MVMD(i,2) = r2(1,2);
    else
        corrs_MVMD(i,1) = r4(1,2);
        corrs_MVMD(i,2) = r3(1,2);
    end

    r1 = corrcoef(mean_C_MSWD_bstrp{i}(1,:),mean_C_MSWD_all{1}(1,:));
    r2 = corrcoef(mean_C_MSWD_bstrp{i}(2,:),mean_C_MSWD_all{1}(2,:));
    r3 = corrcoef(mean_C_MSWD_bstrp{i}(1,:),mean_C_MSWD_all{1}(2,:));
    r4 = corrcoef(mean_C_MSWD_bstrp{i}(2,:),mean_C_MSWD_all{1}(1,:));
    if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
        corrs_MSWD(i,1) = r1(1,2);
        corrs_MSWD(i,2) = r2(1,2);
    else
        corrs_MSWD(i,1) = r4(1,2);
        corrs_MSWD(i,2) = r3(1,2);
    end        
end
[mean_corr_MVMD,~,~,CI_MVMD] = confidence_interval(corrs_MVMD);
[mean_corr_MSWD,~,~,CI_MSWD] = confidence_interval(corrs_MSWD);

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
    exportgraphics(gcf,fullfile(pwd,'plots_real_application',['HCP_',visit,'_states_zero_diag_MVMD.png']),'Resolution',600)
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
    exportgraphics(gcf,fullfile(pwd,'plots_real_application',['HCP_',visit,'_states_zero_diag_MSWD.png']),'Resolution',600)
end