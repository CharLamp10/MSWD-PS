clear
clc

save_res = 1;
numICs = 100;

path_tmaps = 'C:\Users\100063082\Desktop\MSWD_paper_files\dyn_tmap_comps_centered.nii';
path_hcp_average = 'C:\Users\100063082\Desktop\MSWD_paper_files\GICA\GICA_agg__component_ica_.nii';
path_data = 'C:\Users\100063082\Desktop\MSWD_paper_files\GICA';

tmaps = niftiread(path_tmaps);
hcp_average = niftiread(path_hcp_average);

SC = [57,8,13,15];
AUD = [99,98];
SM = [2,10,14,58,37,45,30,9];
VIS = [77,79,60,46,27,43,82,89,80,61];
CC = [67,93,56,90,71,76,47,88,21,59,50,29,66,81];
DM = [20,44,39,28,64,26,75,48,87];
CB = [12,7,32];

[~,~,~,nICs] = size(hcp_average);
[~,~,~,nICs1] = size(tmaps);

cov_term = zeros(nICs,nICs1);
for i = 1:nICs
    resizedVolume = imresize3(hcp_average(:,:,:,i), size(tmaps(:,:,:,1)));
    var_hcp = var(resizedVolume(:));
    for j = 1:nICs1
        tmap = tmaps(:,:,:,j);
        covariance_matrix = corrcoef(resizedVolume(:), tmap(:));
        cov_term(i,j) = abs(covariance_matrix(1, 2));%^2/var_hcp;
    end
end

cov_term_new_max = zeros(nICs,7);
cov_term_new_max(:,1) = max(cov_term(:,SC),[],2);
cov_term_new_max(:,2) = max(cov_term(:,AUD),[],2);
cov_term_new_max(:,3) = max(cov_term(:,SM),[],2);
cov_term_new_max(:,4) = max(cov_term(:,VIS),[],2);
cov_term_new_max(:,5) = max(cov_term(:,CC),[],2);
cov_term_new_max(:,6) = max(cov_term(:,DM),[],2);
cov_term_new_max(:,7) = max(cov_term(:,CB),[],2);

cov_term_new_sum = zeros(nICs,7);
cov_term_new_sum(:,1) = sum(maxk(cov_term(:,SC)',2));
cov_term_new_sum(:,2) = sum(maxk(cov_term(:,AUD)',2));
cov_term_new_sum(:,3) = sum(maxk(cov_term(:,SM)',2));
cov_term_new_sum(:,4) = sum(maxk(cov_term(:,VIS)',2));
cov_term_new_sum(:,5) = sum(maxk(cov_term(:,CC)',2));
cov_term_new_sum(:,6) = sum(maxk(cov_term(:,DM)',2));
cov_term_new_sum(:,7) = sum(maxk(cov_term(:,CB)',2));

maximums = max(cov_term_new_max,[],2);
for i = 1:size(cov_term_new_max,1)
    positions(i) = find(cov_term_new_max(i,:) == maximums(i));
end
[positions_sort,I] = sort(positions);

dir_data = dir(path_data);

fs = 1/0.72;
thresh_same_imf = 10;
min_freq = 0.01;

numSubjects = 20;
compStds = zeros(1,numSubjects);
winds = zeros(1,numSubjects);
minPeaks = zeros(1,numSubjects);
imfs_MSWD = cell(1,numSubjects);
imfs_MVMD = cell(1,numSubjects);

if exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MSWD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
    resMSWD = load(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MSWD_ABIDE_',num2str(numSubjects),'_subjects.mat']);
end
if exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MVMD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
    resMVMD = load(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MVMD_ABIDE_',num2str(numSubjects),'_subjects.mat']);
end

inds_MSWD = [1];

indx = nchoosek(1:numICs,2);

for i = 1:numSubjects
    name = ['GICA_ica_c',num2str(i),'-1.mat'];
    data = load(fullfile(path_data,name));
    data = data.tc;
    %% MSWD
    if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MSWD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
        data_MSWD = data(:,I);
        data_MSWD = data_MSWD./max(abs(data_MSWD));
        p_value = 1e-5;
        P_corr = 1;
        P_corr_imp = 0.02;
        wind = [];
        compStd = 0.05;
        param_struct  = struct('P_corr', P_corr, ...
            'P_corr_imp',   P_corr_imp,...
            'StD_th',       compStd, ...
            'Welch_window', wind, ...
            'p_value',      1e-5);
        imf = MSWD_CL(data_MSWD, param_struct);
        imfs_MSWD{i} = imf;
        disp(i)
    elseif ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MSWD_',num2str(numSubjects),'_subjects.mat'])
        imf = resMSWD.imfs_MSWD{i};
    end
    if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE','\COSDELPHI_MSWD_',num2str(numSubjects),'_subjects.mat'])
        COSDELPHI = phase_sync_analysis_HCP(imf,'MSWD',indx,...
             fs,inds_MSWD,min_freq);
        COSDELPHI1_MSWD{i} = COSDELPHI;
    end
end

for i = 1:numSubjects
    name = ['GICA_ica_c',num2str(i),'-1.mat'];
    data = load(fullfile(path_data,name));
    data = data.tc;
    %% MVMD
    if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MVMD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
        data_MVMD = data(:,I);
        data_MVMD = data_MVMD./max(abs(data_MVMD));
        tau = 0; DC = 1; init = 0; tol = 1e-9;
        K = 10; alpha = 1000;
        imf = MVMD_new(data_MVMD,alpha,tau,K,DC,init,tol);
        imfs_MVMD{i} = imf;
        disp(i)
    end
end
if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MVMD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
    resMVMD.imfs_MVMD = imfs_MVMD;
end
if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MVMD_',num2str(numSubjects),'.mat'])
    for i = 1:length(resMVMD.imfs_MVMD)
        for j = 1:size(resMVMD.imfs_MVMD{i},1)
            if i == 1 && j == 1
                [~,f_mvmd] = pwelch(squeeze(resMVMD.imfs_MVMD{i}(j,:,:)),200,[],[],fs);
            end
            Px(i,j,:) = sum(pwelch(squeeze(resMVMD.imfs_MVMD{i}(j,:,:)),200,[],[],fs),2);
        end
    end
    for i = 1:size(Px,1)
        for j = 1:size(Px,2)
            frequencies_MVMD_sep(i,j) = f_mvmd(Px(i,j,:) == max(Px(i,j,:)));
        end
    end
    Px = squeeze(sum(Px,1))';
    for i = 1:size(Px,2)
        frequencies_MVMD(i) = f_mvmd(Px(:,i) == max(Px(:,i)));
    end
    Pxs_MVMD = Px;
    [~,pos_maxs] = max(Pxs_MVMD);
    f_maxs = f_mvmd(pos_maxs);
    [~,inds_MVMD] = min(abs(f_maxs-0.05));
end
for i = 1:numSubjects
    if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MVMD_',num2str(numSubjects),'_subjects.mat'])
        COSDELPHI = phase_sync_analysis_HCP(resMVMD.imfs_MVMD{i},'MVMD',indx,...
             fs,inds_MVMD,min_freq);
        COSDELPHI1_MVMD{i} = COSDELPHI;
    end
end

if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MSWD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
    save(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MSWD_ABIDE_',num2str(numSubjects),'_subjects.mat'],"imfs_MSWD")
end

if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MSWD_',num2str(numSubjects),'_subjects.mat'])
    save(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MSWD_',num2str(numSubjects),'_subjects.mat'],'COSDELPHI1_MSWD','-v7.3')
end


if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MVMD_ABIDE_',num2str(numSubjects),'_subjects.mat'])
    save(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\MVMD_ABIDE_',num2str(numSubjects),'_subjects.mat'],"imfs_MVMD")
end

if ~exist(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MVMD_',num2str(numSubjects),'_subjects.mat'])
    save(['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE\COSDELPHI_MVMD_',num2str(numSubjects),'_subjects.mat'],'COSDELPHI1_MVMD','-v7.3')
end



