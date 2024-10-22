function [plv_new,counts] = sort_values(mfreq,freq_ranges,plv)

plv_new = zeros(length(mfreq),length(freq_ranges));
counts = zeros(length(mfreq),length(freq_ranges));
for i = 1:length(mfreq)
    if ~isnan(mfreq{i})
        for j = 1:length(mfreq{i})
            [~,ind] = min(abs(freq_ranges - mfreq{i}(j)));
            plv_new(i,ind) = plv_new(i,ind) + plv(i,j);
            counts(i,ind) = counts(i,ind) + 1;
        end
    else
        plv_new(i,:) = zeros(1,length(freq_ranges)) - 121;
        counts(i,:) = counts(i,:) + 1;
    end
end
counts(counts == 0) = 1;
plv_new = plv_new./counts;
% plv_new(plv_new == 0) = nan;
% correlation_new(correlation_new == 0) = nan;
plv_new(plv_new == -121) = 0;
end