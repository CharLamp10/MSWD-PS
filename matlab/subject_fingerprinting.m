clear 
clc

casee = 'MSWD';
cat = '100';
save_results = 1;

path_data1 = ['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_HCP_rest1\',casee,'_HCP_98_subjects.mat'];
path_data2 = ['F:\MSWD_paper_new\MSWD_CL_paper_decomposed_HCP_rest2\',casee,'_HCP_98_subjects.mat'];
visit1 = load(path_data1);
visit2 = load(path_data2);

load("names_98_subjects.mat")
numSubjects = length(names);

if casee == "MVMD"
    imfs1 = visit1.imfs_MVMD;
    imfs2 = visit2.imfs_MVMD;
elseif casee == "MSWD"
    imfs1 = visit1.imfs_MSWD;
    imfs2 = visit2.imfs_MSWD;
end

if contains(pwd,'100063082')
    path_tmaps = 'C:\Users\100063082\Desktop\MSWD_paper_files\dyn_tmap_comps_centered.nii';
    path_hcp_average = ['C:\Users\100063082\Desktop\MSWD_paper_files\HCP_PTN1200\groupICA\groupICA_3T_HCP1200_MSMAll_d',cat,'.ica\melodic_IC_sum.nii'];
    path_data = ['C:\Users\100063082\Desktop\MSWD_paper_files\HCP_PTN1200\NodeTimeseries_3T_HCP1200_MSMAll_ICAd',cat,'_ts2','\node_timeseries','\3T_HCP1200_MSMAll_d',cat,'_ts2'];
else
    path_tmaps = 'C:\Users\Hp\Desktop\MSWD_paper_files\dyn_tmap_comps_centered.nii';
    path_hcp_average = ['C:\Users\Hp\Desktop\MSWD_paper_files\HCP_PTN1200\groupICA\groupICA_3T_HCP1200_MSMAll_d',cat,'.ica\melodic_IC_sum.nii'];
    path_data = ['C:\Users\Hp\Desktop\MSWD_paper_files\HCP_PTN1200\NodeTimeseries_3T_HCP1200_MSMAll_ICAd',cat,'_ts2','\node_timeseries','\3T_HCP1200_MSMAll_d',cat,'_ts2'];
end

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

for i = 1:length(names)
    name = names{i};
    data = load(fullfile(path_data,name));
    data1s{i} = data(1:1200,I);
    data1s{i} = data1s{i}./max(abs(data1s{i}));
    data2s{i} = data(2401:3600,I);
    data2s{i} = data2s{i}./max(abs(data2s{i}));
end


factors = [0.1:0.05:1,1.1:0.1:2,2.5:0.5:10];
correct_MSWD = zeros(1,length(factors));
correct = zeros(1,length(factors));
for k = 1:length(factors)
    for i = 1:length(imfs1)
        rec1 = sum(imfs1{i}(:,:,1:end-1),3);
        rec2 = sum(imfs2{i}(:,:,1:end-1),3);
        rec3 = imfs1{i}(:,:,end);
        rec4 = imfs2{i}(:,:,end);
        rec_v1 = rec1.*factors(k) + rec3;
        rec_v2 = rec2.*factors(k) + rec4;
        data1 = data1s{i};
        data2 = data2s{i};
        corr1 = corrcoef(data1);
        corr1 = triu(corr1,1);
        corr1(corr1 == 0) = [];
        corr2 = corrcoef(data2);
        corr2 = triu(corr2,1);
        corr2(corr2 == 0) = [];
        corr_dec1 = corrcoef(rec_v1);
        corr_dec1 = triu(corr_dec1,1);
        corr_dec1(corr_dec1 == 0) = [];
        corr_dec2 = corrcoef(rec_v2);
        corr_dec2 = triu(corr_dec2,1);
        corr_dec2(corr_dec2 == 0) = [];
        corrs_dec1(i,:) = corr_dec1;
        corrs_dec2(i,:) = corr_dec2;
        corrs1(i,:) = corr1;
        corrs2(i,:) = corr2;
    end
    
    for l = 1:length(imfs1)
        for n = 1:length(imfs1)
            temp = corrcoef(corrs_dec1(l,:),corrs_dec2(n,:));
            res(n) = temp(1,2);
        end
        [~,pos_max] = max(res);
        if pos_max == l
            correct_MSWD(k) = correct_MSWD(k) + 1;
        end
    end

    for l = 1:length(imfs1)
        for n = 1:length(imfs1)
            temp = corrcoef(corrs1(l,:),corrs2(n,:));
            res(n) = temp(1,2);
        end
        [~,pos_max] = max(res);
        if pos_max == l
            correct(k) = correct(k) + 1;
        end
    end
end

plot(log10(factors),correct_MSWD./98,'LineWidth',1.2)
hold on
plot(log10(factors),correct./98,'LineWidth',1.2)
ax = gca;
ax.FontSize = 11; 
xlim([log10(factors(1)),log10(factors(end))])
xlabel('Log_1_0(λ)','FontSize',14)
ylabel('Subject Identification Rate (%)','FontSize',14)
exportgraphics(gcf,'C:\Users\100063082\Desktop\Dissertation\MSWD_paper_new\paper_plots\SI_rate.png','Resolution',1000)
