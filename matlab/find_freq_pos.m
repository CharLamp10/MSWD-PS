function positions = find_freq_pos(mfreq,f)

positions = nan(1,length(f));
for i = 1:length(f)
    for j = 1:length(mfreq)
        dist(i,j) = abs(f(i) - mfreq(j));
    end
end

while any(isnan(positions)) && ~all(all(dist == 1e5))
    [r,c] = find(dist == min(min(dist)));
    positions(r) = c;
    dist(r,:) = 1e5;
    dist(:,c) = 1e5;
end
end