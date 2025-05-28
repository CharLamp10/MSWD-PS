clear
clc

path_rois = 'C:\Users\100063082\Desktop\MSWD_PS_revision\rois_ho';
dir_rois = dir(path_rois);
count  = 1;
for i = 1:length(dir_rois) - 2
    name = dir_rois(i+2).name;
    if ~contains(name,'OHSU') && contains(name,'.1D')
        path = fullfile(path_rois,name);
        % data = importdata(path).data;
        names(count) = convertCharsToStrings(name);
        % save([path(1:end-3),'.mat'], "data")
        count = count + 1;
    end
end

path_MSWD = 'F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE_rois\MSWD_ABIDE_rois.mat';
load(path_MSWD);
save_path_MSWD = 'C:\Users\100063082\Desktop\MSWD_PS_revision\MSWD_rois_ho';
if ~exist(save_path_MSWD)
    mkdir(save_path_MSWD)
end
for i = 1:length(imfs_MSWD)
    name = convertStringsToChars(names(i));
    imf = imfs_MSWD{i};
    global_x = sum(imf(:,:,1:end-1),3);
    local_x = imf(:,:,end);
    rec = global_x.*0.5 + local_x;
    % save(fullfile(save_path_MSWD,[name(1:end-3),'.mat']), "rec")
end