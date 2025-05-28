clear
clc

data_path = 'C:\Users\100063082\Desktop\MSWD_PS_revision\rois_ho';
data_dir = dir(data_path);
save_path_corr = 'C:\Users\100063082\Desktop\MSWD_PS_revision\rois_ho_corr';
save_path_corr_MSWD = 'C:\Users\100063082\Desktop\MSWD_PS_revision\MSWD_rois_ho_corr';
save_MSWD_path = 'F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE_rois';
if ~exist(save_path_corr)
    mkdir(save_path_corr)
end
if ~exist(save_path_corr_MSWD)
    mkdir(save_path_corr_MSWD)
end
if ~exist(save_MSWD_path)
    mkdir(save_MSWD_path)
end

if exist(fullfile(save_MSWD_path,'MSWD_ABIDE_rois.mat'))
    load(fullfile(save_MSWD_path,'MSWD_ABIDE_rois.mat'))
else
    imfs_MSWD = cell(1,length(data_dir)-2);
    T_MSWD = zeros(1,length(data_dir)-2);
end

for i = 1:length(data_dir)-2
    name = data_dir(i+2).name;
    path = fullfile(data_path,name);
    data = importdata(path).data;
    if size(data,1) > size(data,2)
        %% Conventional PCC
        if ~exist(fullfile(save_path_corr,[name(1:end-10),'corr_ho.mat']))
            corr = corrcoef(data);
            save(fullfile(save_path_corr,[name(1:end-10),'corr_ho.mat']),'corr')
        end
        %% MSWD
        if isempty(imfs_MSWD{i})
            R = size(data,2);
            data_MSWD = data./max(abs(data));
            [r,c] = find(isnan(data_MSWD));
            C = [];
            for j = 1:R
                if sum(c == j) == 0
                    C = [C,j];
                end
            end
            data_MSWD(:,c) = [];
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
            tic
            imf = MSWD(data_MSWD, param_struct);
            if ~isempty(c)
                imf_temp = zeros(size(data_MSWD,1),R,size(imf,3));
                imf_temp(:,C,:) = imf;
                imf = imf_temp;
            end
            T_MSWD(i) = toc;
            imfs_MSWD{i} = imf;
        else
            imf = imfs_MSWD{i};
        end
         if ~exist(fullfile(save_path_corr_MSWD,[name(1:end-10),'corr_ho.mat'])) && ~all(imf==0,"all")
            global_x = sum(imf(:,:,1:end-1),3);
            local_x = imf(:,:,end);
            rec = global_x.*0.5 + local_x;
            corr_dec = corrcoef(rec);
            save(fullfile(save_path_corr_MSWD,[name(1:end-10),'corr_ho.mat']),'corr_dec')
         end
         disp(i)
    end
end

if ~exist(fullfile(save_MSWD_path,'MSWD_ABIDE_rois.mat'))
    save(fullfile(save_MSWD_path,'MSWD_ABIDE_rois.mat'),'imfs_MSWD','T_MSWD')
end