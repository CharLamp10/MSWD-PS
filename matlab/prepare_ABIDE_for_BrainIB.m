clear
clc

path_rois = 'C:\Users\100063082\Desktop\MSWD_PS_revision\rois_ho';
dir_rois = dir(path_rois);
phenotypic = readtable('C:\Users\100063082\Desktop\prepare_ADHD200_ABIDE_phenotypics\Phenotypic_V1_0b_preprocessed1.csv');
subid = table2array(phenotypic(:,3));
groups = table2array(phenotypic(:,8));

count = 1;
for i = 1:length(dir_rois) - 2
    name = dir_rois(i+2).name;
    if ~contains(name,'OHSU') && contains(name,'.1D')
        names(count) = convertCharsToStrings(name);
        count = count + 1;
    end
end
for n = 1:length(names)
    name = names(n);
    name = convertStringsToChars(name);
    pos = strfind(name,'00');
    ids(n) = str2num(name(pos+2:pos+6));
end
for g = 1:length(ids)
    id = ids(g);
    pos = find(subid == id);
    label(g) = groups(pos);
end

save_path = 'C:\Users\100063082\Desktop\MSWD_PS_revision';

for i = 1:length(names)
    res(i) = struct('edge', [], 'ROI', [], 'node_tags', [], 'adj', [], 'neighbor', []);
end
count  = 1;
for i = 1:length(dir_rois) - 2
    name = dir_rois(i+2).name;
    if ~contains(name,'OHSU') && contains(name,'.1D')
        path = fullfile(path_rois,name);
        data = importdata(path).data;
        corr = corrcoef(data);
        corr(isnan(corr)) = 0;
        abs_corr = abs(corr);
        adj = abs_corr .* (abs_corr > prctile(abs_corr(:),80));
        adj(adj == 1) = 0;
        adj(adj ~= 0) = 1;
        corr_triu = triu(corr,1);
        [r,c] = find(abs_corr > prctile(abs_corr(:),80));
        r = r-1;
        c = c-1;
        rows = [];
        columns = [];
        vals = [];
        for j = 1:size(corr,1)
            pos = find(r == j-1);
            rows = [rows; r(pos)];
            columns = [columns; c(pos)];
            temp = reshape(1:size(corr_triu,1)^2, size(corr_triu,1), size(corr_triu,1));
            linear_idx = sub2ind(size(corr_triu), r(pos)+1, c(pos)+1);
            vals = [vals; temp(linear_idx);];
            neighbors{j,1} = c(pos);
        end
        edge = [rows,columns,vals];
        corr(corr == 1) = 0;
        node_tags = 1:size(corr,1);
        res(count).edge = edge;
        res(count).ROI = corr;
        res(count).node_tags = node_tags;
        res(count).adj = adj;
        res(count).neighbor = neighbors;
        count = count + 1;
    end
end
graph_struct = res;
save(fullfile(save_path,'md_AAL_0.4.mat'), "graph_struct", "label")

path_MSWD = 'F:\MSWD_paper_new\MSWD_CL_paper_decomposed_ABIDE_rois\MSWD_ABIDE_rois.mat';
load(path_MSWD);
save_path_MSWD = 'C:\Users\100063082\Desktop\MSWD_PS_revision\MSWD_rois_ho';
if ~exist(save_path_MSWD)
    mkdir(save_path_MSWD)
end

for i = 1:length(names)
    res(i) = struct('edge', [], 'ROI', [], 'node_tags', [], 'adj', [], 'neighbor', []);
end
count = 1;
for i = 1:length(imfs_MSWD)
    name = convertStringsToChars(names(i));
    imf = imfs_MSWD{i};
    global_x = sum(imf(:,:,1:end-1),3);
    local_x = imf(:,:,end);
    rec = global_x.*0.5 + local_x;
    corr = corrcoef(rec);
    corr(isnan(corr)) = 0;
    abs_corr = abs(corr);
    adj = abs_corr .* (abs_corr > prctile(abs_corr(:),80));
    adj(adj == 1) = 0;
    adj(adj ~= 0) = 1;
    corr_triu = triu(corr,1);
    [r,c] = find(abs_corr > prctile(abs_corr(:),80));
    r = r-1;
    c = c-1;
    rows = [];
    columns = [];
    vals = [];
    for j = 1:size(corr,1)
        pos = find(r == j-1);
        rows = [rows; r(pos)];
        columns = [columns; c(pos)];
        temp = reshape(1:size(corr_triu,1)^2, size(corr_triu,1), size(corr_triu,1));
        linear_idx = sub2ind(size(corr_triu), r(pos)+1, c(pos)+1);
        vals = [vals; temp(linear_idx);];
        neighbors{j,1} = c(pos);
    end
    edge = [rows,columns,vals];
    corr(corr == 1) = 0;
    node_tags = 1:size(corr,1);
    res(count).edge = edge;
    res(count).ROI = corr;
    res(count).node_tags = node_tags;
    res(count).adj = adj;
    res(count).neighbor = neighbors;
    count = count + 1;
end
graph_struct = res;
save(fullfile(save_path,'MSWD_md_AAL_0.4.mat'), "graph_struct", "label")