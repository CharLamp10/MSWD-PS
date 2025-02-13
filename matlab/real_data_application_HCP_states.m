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
save_res = false;

%% Load TVPS data and concatenate them across subjects
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MEMD_50_subjects.mat']))
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_EWT_50_subjects.mat']))
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_IRCNN_50_subjects.mat']))
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MVMD_50_subjects.mat']))
load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\COSDELPHI_MSWD_50_subjects.mat']))

for i = 1:length(COSDELPHI1_MSWD)
    if i == 1
        COSDELPHI_MEMD = COSDELPHI1_MEMD{i};
        COSDELPHI_EWT = COSDELPHI1_EWT{i};
        COSDELPHI_IRCNN = COSDELPHI1_IRCNN{i};
        COSDELPHI_MVMD = COSDELPHI1_MVMD{i};
        COSDELPHI_MSWD = COSDELPHI1_MSWD{i};
    else
        COSDELPHI_MEMD = cat(1,COSDELPHI_MEMD,COSDELPHI1_MEMD{i});
        COSDELPHI_EWT = cat(1,COSDELPHI_EWT,COSDELPHI1_EWT{i});
        COSDELPHI_IRCNN = cat(1,COSDELPHI_IRCNN,COSDELPHI1_IRCNN{i});
        COSDELPHI_MVMD = cat(1,COSDELPHI_MVMD,COSDELPHI1_MVMD{i});
        COSDELPHI_MSWD = cat(1,COSDELPHI_MSWD,COSDELPHI1_MSWD{i});
    end
end

%% Find optimal number of clusters
% We found that optimal number of clusters is 2. To validate this,
% uncomment the following lines of code
% E_MEMD = evalclusters(squeeze(COSDELPHI_MEMD),'kmeans','silhouette','klist',[2:5]);
% DBI_MEMD = E_MEMD.OptimalK;
% 
% E_EWT = evalclusters(squeeze(COSDELPHI_EWT),'kmeans','silhouette','klist',[2:5]);
% DBI_EWT = E_EWT.OptimalK;
% 
% E_IRCNN = evalclusters(squeeze(COSDELPHI_IRCNN),'kmeans','silhouette','klist',[2:5]);
% DBI_IRCNN = E_IRCNN.OptimalK;
% 
% E_MVMD = evalclusters(squeeze(COSDELPHI_MVMD),'kmeans','silhouette','klist',[2:5]);
% DBI_MVMD = E_MVMD.OptimalK;
% 
% E_MSWD = evalclusters(squeeze(COSDELPHI_MSWD),'kmeans','silhouette','klist',[2:5]);
% DBI_MSWD = E_MSWD.OptimalK;
DBI_MEMD = 2;
DBI_EWT = 2;
DBI_IRCNN = 2;
DBI_MVMD = 2;
DBI_MSWD = 2;

%% Perform clustering and get the centroids and the clustered data
if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_all.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_all.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_all.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']))
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_all.mat']))
        [~,C_init_MEMD,sumd_MEMD] = kmeans(squeeze(COSDELPHI_MEMD),DBI_MEMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_MEMD_all{1},mean_C_MEMD_all{1}] = kmeans(squeeze(COSDELPHI_MEMD),DBI_MVMD,'MaxIter',500,'Start',C_init_MEMD); %1000 default
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_all.mat']))
        [~,C_init_EWT,sumd_EWT] = kmeans(squeeze(COSDELPHI_EWT),DBI_EWT,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_EWT_all{1},mean_C_EWT_all{1}] = kmeans(squeeze(COSDELPHI_EWT),DBI_EWT,'MaxIter',500,'Start',C_init_EWT); %1000 default
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_all.mat']))
        [~,C_init_IRCNN,sumd_IRCNN] = kmeans(squeeze(COSDELPHI_IRCNN),DBI_IRCNN,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_IRCNN_all{1},mean_C_IRCNN_all{1}] = kmeans(squeeze(COSDELPHI_IRCNN),DBI_IRCNN,'MaxIter',500,'Start',C_init_IRCNN); %1000 default
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']))
        [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_MVMD_all{1},mean_C_MVMD_all{1}] = kmeans(squeeze(COSDELPHI_MVMD),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']))
        [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
        [mean_idx_MSWD_all{1},mean_C_MSWD_all{1}] = kmeans(squeeze(COSDELPHI_MSWD),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_all.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_all.mat']),'mean_idx_MEMD_all','mean_C_MEMD_all');
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_all.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_all.mat']),'mean_idx_EWT_all',"mean_C_EWT_all");
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_all.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_all.mat']),'mean_idx_IRCNN_all',"mean_C_IRCNN_all");
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']),'mean_idx_MVMD_all','mean_C_MVMD_all');
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']),'mean_idx_MSWD_all',"mean_C_MSWD_all");
    end
else
    if ~exist("mean_C_MEMD_all")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_all.mat']));
    end
    if ~exist("mean_C_EWT_all")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_all.mat']));
    end
    if ~exist("mean_C_IRCNN_all")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_all.mat']));
    end
    if ~exist("mean_C_MVMD_all")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_all.mat']));
    end
    if ~exist("mean_C_MSWD_all")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_all.mat']));
    end
end


%% Count short state transitions
count = 1;
for i = 1:signal_length:length(mean_idx_MEMD_all{1})
    diff_states = diff(mean_idx_MEMD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MEMD(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
    clear pos
end

count = 1;
for i = 1:signal_length:length(mean_idx_EWT_all{1})
    diff_states = diff(mean_idx_EWT_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_EWT(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
    clear pos
end

count = 1;
for i = 1:signal_length:length(mean_idx_IRCNN_all{1})
    diff_states = diff(mean_idx_IRCNN_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_IRCNN(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
    clear pos
end

count = 1;
for i = 1:signal_length:length(mean_idx_MVMD_all{1})
    diff_states = diff(mean_idx_MVMD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MVMD(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
    clear pos
end

count = 1;
for i = 1:signal_length:length(mean_idx_MSWD_all{1})
    diff_states = diff(mean_idx_MSWD_all{1}(i:i+signal_length-1));
    pos = find(diff_states ~= 0);
    pos = [1;pos;signal_length];
    diff_poss = diff(pos);
    count_0_15_sec_MSWD(count) = sum((diff_poss>0).*(diff_poss<=21));
    count = count + 1;
    clear pos
end
disp(['MEMD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MEMD))])
disp(['EWT-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_EWT))])
disp(['IRCNN-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_IRCNN))])
disp(['MVMD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MVMD))])
disp(['MSwD-PS total biologically implausible transitions: ', num2str(sum(count_0_15_sec_MSWD))])


%% Bootstrap analysis to show that the centroids/brain states are robust
if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_btsrp.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_btsrp.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_btsrp.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat'])) ||...
    ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']))
    for i = 1:25
        N = randperm(numSubjects,10);
        mask = zeros(size(COSDELPHI_MSWD,1),1);
        for n = N
            start = (n-1)*signal_length+1;
            endd = signal_length*n;
            mask(start:endd) = 1;
        end
        COSDELPHI_MEMD_btsrp = COSDELPHI_MEMD(logical(mask),:);
        COSDELPHI_EWT_btsrp = COSDELPHI_EWT(logical(mask),:);
        COSDELPHI_IRCNN_btsrp = COSDELPHI_IRCNN(logical(mask),:);
        COSDELPHI_MVMD_btsrp = COSDELPHI_MVMD(logical(mask),:);
        COSDELPHI_MSWD_btsrp = COSDELPHI_MSWD(logical(mask),:);
        if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_btsrp.mat']))
            [~,C_init_MEMD,sumd_MEMD] = kmeans(squeeze(COSDELPHI_MEMD_btsrp),DBI_MEMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_MEMD_bstrp{i},mean_C_MEMD_bstrp{i}] = kmeans(squeeze(COSDELPHI_MEMD_btsrp),DBI_MEMD,'MaxIter',500,'Start',C_init_MEMD); %1000 default
        end
        if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_btsrp.mat']))
            [~,C_init_EWT,sumd_EWT] = kmeans(squeeze(COSDELPHI_EWT_btsrp),DBI_EWT,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_EWT_bstrp{i},mean_C_EWT_bstrp{i}] = kmeans(squeeze(COSDELPHI_EWT_btsrp),DBI_EWT,'MaxIter',500,'Start',C_init_EWT); %1000 default
        end
        if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_btsrp.mat']))
            [~,C_init_IRCNN,sumd_IRCNN] = kmeans(squeeze(COSDELPHI_IRCNN_btsrp),DBI_IRCNN,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_IRCNN_bstrp{i},mean_C_IRCNN_bstrp{i}] = kmeans(squeeze(COSDELPHI_IRCNN_btsrp),DBI_IRCNN,'MaxIter',500,'Start',C_init_IRCNN); %1000 default
        end
        if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat']))
            [~,C_init_MVMD,sumd_MVMD] = kmeans(squeeze(COSDELPHI_MVMD_btsrp),DBI_MVMD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_MVMD_bstrp{i},mean_C_MVMD_bstrp{i}] = kmeans(squeeze(COSDELPHI_MVMD_btsrp),DBI_MVMD,'MaxIter',500,'Start',C_init_MVMD); %1000 default
        end
        if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']))
            [~,C_init_MSWD,sumd_MSWD] = kmeans(squeeze(COSDELPHI_MSWD_btsrp),DBI_MSWD,'MaxIter',150,'Start','sample','Replicates',10); % replicates default 200, MaxIter 150, replicates 10
            [mean_idx_MSWD_bstrp{i},mean_C_MSWD_bstrp{i}] = kmeans(squeeze(COSDELPHI_MSWD_btsrp),DBI_MSWD,'MaxIter',500,'Start',C_init_MSWD); %1000 default
        end
        disp(i)
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_btsrp.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_btsrp.mat']),"mean_C_MEMD_bstrp")
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_btsrp.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_btsrp.mat']),"mean_C_EWT_bstrp")
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_btsrp.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_btsrp.mat']),"mean_C_IRCNN_bstrp")
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat']),"mean_C_MVMD_bstrp")
    end
    if ~exist(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']))
        save(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']),"mean_C_MSWD_bstrp")
    end
else
    if ~exist("mean_C_MEMD_bstrp")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MEMD_btsrp.mat']))
    end
    if ~exist("mean_C_EWT_bstrp")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_EWT_btsrp.mat']))
    end
    if ~exist("mean_C_IRCNN_bstrp")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_IRCNN_btsrp.mat']))
    end
    if ~exist("mean_C_MVMD_bstrp")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MVMD_btsrp.mat']))
    end
    if ~exist("mean_C_MSWD_bstrp")
        load(fullfile(path_data,['MSWD_CL_paper_decomposed_HCP_',visit,'\mean_C_MSWD_btsrp.mat']))
    end
end

%% Correlation coefficient between the bootstrap centroids and the original centroids
for i = 1:size(mean_C_MSWD_bstrp,1)
    r1 = corrcoef(mean_C_MEMD_bstrp{i}(1,:),mean_C_MEMD_all{1}(1,:));
    r2 = corrcoef(mean_C_MEMD_bstrp{i}(2,:),mean_C_MEMD_all{1}(2,:));
    r3 = corrcoef(mean_C_MEMD_bstrp{i}(1,:),mean_C_MEMD_all{1}(2,:));
    r4 = corrcoef(mean_C_MEMD_bstrp{i}(2,:),mean_C_MEMD_all{1}(1,:));
    if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
        corrs_MEMD(i,1) = r1(1,2);
        corrs_MEMD(i,2) = r2(1,2);
    else
        corrs_MEMD(i,1) = r4(1,2);
        corrs_MEMD(i,2) = r3(1,2);
    end

    r1 = corrcoef(mean_C_EWT_bstrp{i}(1,:),mean_C_EWT_all{1}(1,:));
    r2 = corrcoef(mean_C_EWT_bstrp{i}(2,:),mean_C_EWT_all{1}(2,:));
    r3 = corrcoef(mean_C_EWT_bstrp{i}(1,:),mean_C_EWT_all{1}(2,:));
    r4 = corrcoef(mean_C_EWT_bstrp{i}(2,:),mean_C_EWT_all{1}(1,:));
    if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
        corrs_EWT(i,1) = r1(1,2);
        corrs_EWT(i,2) = r2(1,2);
    else
        corrs_EWT(i,1) = r4(1,2);
        corrs_EWT(i,2) = r3(1,2);
    end

    r1 = corrcoef(mean_C_IRCNN_bstrp{i}(1,:),mean_C_IRCNN_all{1}(1,:));
    r2 = corrcoef(mean_C_IRCNN_bstrp{i}(2,:),mean_C_IRCNN_all{1}(2,:));
    r3 = corrcoef(mean_C_IRCNN_bstrp{i}(1,:),mean_C_IRCNN_all{1}(2,:));
    r4 = corrcoef(mean_C_IRCNN_bstrp{i}(2,:),mean_C_IRCNN_all{1}(1,:));
    if r1(1,2) + r2(1,2) > r3(1,2) + r4(1,2)
        corrs_IRCNN(i,1) = r1(1,2);
        corrs_IRCNN(i,2) = r2(1,2);
    else
        corrs_IRCNN(i,1) = r4(1,2);
        corrs_IRCNN(i,2) = r3(1,2);
    end

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
[mean_corr_MEMD,~,~,CI_MEMD] = confidence_interval(corrs_MEMD);
[mean_corr_EWT,~,~,CI_EWT] = confidence_interval(corrs_EWT);
[mean_corr_IRCNN,~,~,CI_IRCNN] = confidence_interval(corrs_IRCNN);
[mean_corr_MVMD,~,~,CI_MVMD] = confidence_interval(corrs_MVMD);
[mean_corr_MSWD,~,~,CI_MSWD] = confidence_interval(corrs_MSWD);

%% Plot the centroids/brain states for the two methods
figure;
h = tight_subplot(1, DBI_MEMD, [0.001 0.02],[.001 .001],[.05 .1]);
for i = 1:size(mean_C_MEMD_all{1},1)
    state = mean_C_MEMD_all{1}(i,:);
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
    if i == DBI_MEMD
        colorbar(h(i),'Position',[0.93 0.225 0.022 0.548]); caxis([-1 1])
    end
end
if save_res
    exportgraphics(gcf,fullfile(pwd,'plots_real_application_new',['HCP_',visit,'_states_zero_diag_MEMD.png']),'Resolution',600)
end

figure;
h = tight_subplot(1, DBI_EWT, [0.001 0.02],[.001 .001],[.05 .1]);
for i = 1:size(mean_C_EWT_all{1},1)
    state = mean_C_EWT_all{1}(i,:);
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
    if i == DBI_EWT
        colorbar(h(i),'Position',[0.93 0.225 0.022 0.548]); caxis([-1 1])
    end
end
if save_res
    exportgraphics(gcf,fullfile(pwd,'plots_real_application_new',['HCP_',visit,'_states_zero_diag_EWT.png']),'Resolution',600)
end

figure;
h = tight_subplot(1, DBI_IRCNN, [0.001 0.02],[.001 .001],[.05 .1]);
for i = 1:size(mean_C_IRCNN_all{1},1)
    state = mean_C_IRCNN_all{1}(i,:);
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
    if i == DBI_IRCNN
        colorbar(h(i),'Position',[0.93 0.225 0.022 0.548]); caxis([-1 1])
    end
end
if save_res
    exportgraphics(gcf,fullfile(pwd,'plots_real_application_new',['HCP_',visit,'_states_zero_diag_IRCNN.png']),'Resolution',600)
end

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
    exportgraphics(gcf,fullfile(pwd,'plots_real_application_new',['HCP_',visit,'_states_zero_diag_MVMD.png']),'Resolution',600)
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
    exportgraphics(gcf,fullfile(pwd,'plots_real_application_new',['HCP_',visit,'_states_zero_diag_MSWD.png']),'Resolution',600)
end