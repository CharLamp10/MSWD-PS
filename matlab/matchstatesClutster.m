% This function matches each cluster obtained from K-means to its
% correponding ground truth, based on maximum similarity

function [mCorr,midx] = matchstatesClutster(grdstate,Corr,idx,nS)

if nS == 3
    for i = 1:size(grdstate,1)
        correlation(1,i) = sum(sum(abs(grdstate(:,:,1) - Corr(:,:,i))));
        correlation(2,i) = sum(sum(abs(grdstate(:,:,2) - Corr(:,:,i))));
        correlation(3,i) = sum(sum(abs(grdstate(:,:,3) - Corr(:,:,i))));
    end
    
    [pos1_val,pos1] = min(correlation(1,:));
    [pos2_val,pos2] = min(correlation(2,:));
    [pos3_val,pos3] = min(correlation(3,:));
    if (pos1 == pos2 || pos1 == pos3 || pos2 == pos3)
        warning('State not discovered')
        if pos1 == pos2 && pos1 == pos3 && pos2 == pos3
            [~,pos_min] = min([pos1_val,pos2_val,pos3_val]);
            if pos_min == 1
                if pos2_val < pos3_val
                    pos2 = min(correlation(2,setdiff([1,2,3],[pos1,pos2,pos3])));
                    pos3 = setdiff([1,2,3],[pos1,pos2]);
                else
                    pos3 = min(correlation(3,setdiff([1,2,3],[pos1,pos2,pos3])));
                    pos2 = setdiff([1,2,3],[pos1,pos3]);
                end
            elseif pos_min == 2
                if pos1_val < pos3_val
                    pos1 = min(correlation(1,setdiff([1,2,3],[pos1,pos2,pos3])));
                    pos3 = setdiff([1,2,3],[pos2,pos1]);
                else
                    pos3 = min(correlation(3,setdiff([1,2,3],[pos1,pos2,pos3])));
                    pos1 = setdiff([1,2,3],[pos2,pos3]);
                end
            elseif pos_min == 3
                if pos1_val < pos2_val
                    pos1 = min(correlation(1,setdiff([1,2,3],[pos1,pos2,pos3])));
                    pos2 = setdiff([1,2,3],[pos3,pos1]);
                else
                    pos2 = min(correlation(2,setdiff([1,2,3],[pos1,pos2,pos3])));
                    pos1 = setdiff([1,2,3],[pos2,pos3]);
                end
            end
        elseif pos1 == pos2 && pos1 ~= pos3 && pos2 ~= pos3
            if pos1_val < pos2_val
                pos2 = setdiff([1,2,3],[pos1,pos3]);
            else
                pos1 = setdiff([1,2,3],[pos2,pos3]);
            end
        elseif pos1 ~= pos2 && pos1 == pos3 && pos2 ~= pos3
            if pos1_val < pos3_val
                pos3 = setdiff([1,2,3],[pos1,pos2]);
            else
                pos1 = setdiff([1,2,3],[pos2,pos3]);
            end
        elseif pos1 ~= pos2 && pos1 ~= pos3 && pos2 == pos3
            if pos2_val < pos3_val
                pos3 = setdiff([1,2,3],[pos1,pos2]);
            else
                pos2 = setdiff([1,2,3],[pos1,pos3]);
            end
        end
        mCorr(:,:,1) = Corr(:,:,pos1);
        mCorr(:,:,2) = Corr(:,:,pos2);
        mCorr(:,:,3) = Corr(:,:,pos3);
    else
        mCorr(:,:,1) = Corr(:,:,pos1);
        mCorr(:,:,2) = Corr(:,:,pos2);
        mCorr(:,:,3) = Corr(:,:,pos3);
    end
    midx = nan(size(idx));
    midx(idx == pos1) = -1;
    midx(idx == pos2) = -2;
    midx(idx == pos3) = -3;
    
    midx(midx == -1) = 1;
    midx(midx == -2) = 2;
    midx(midx == -3) = 3;
elseif nS == 4
    for i = 1:size(grdstate,3)
        correlation(1,i) = sum(sum(abs(grdstate(:,:,1) - Corr(:,:,i))));
        correlation(2,i) = sum(sum(abs(grdstate(:,:,2) - Corr(:,:,i))));
        correlation(3,i) = sum(sum(abs(grdstate(:,:,3) - Corr(:,:,i))));
        correlation(4,i) = sum(sum(abs(grdstate(:,:,4) - Corr(:,:,i))));
    end
    
    [pos1_val,pos1] = min(correlation(1,:));
    [pos2_val,pos2] = min(correlation(2,:));
    [pos3_val,pos3] = min(correlation(3,:));
    [pos4_val,pos4] = min(correlation(4,:));
    if length(unique([pos1,pos2,pos3,pos4])) ~= 4
        for i = 1:4
            for j = 1:4
                for n = 1:4
                    for m = 1:4
                        if length(unique([i,j,n,m])) == 4
                            dist(i,j,n,m) = correlation(1,i) + correlation(2,j) + correlation(3,n) + correlation(4,m);
                        else
                            dist(i,j,n,m) = 1e10;
                        end
                    end
                end
            end
        end
        [minval, minidx] = min(dist(:));
        [pos1,pos2,pos3,pos4] = ind2sub(size(dist),minidx);
        mCorr(:,:,1) = Corr(:,:,pos1);
        mCorr(:,:,2) = Corr(:,:,pos2);
        mCorr(:,:,3) = Corr(:,:,pos3);
        mCorr(:,:,4) = Corr(:,:,pos4);
    else
        mCorr(:,:,1) = Corr(:,:,pos1);
        mCorr(:,:,2) = Corr(:,:,pos2);
        mCorr(:,:,3) = Corr(:,:,pos3);
        mCorr(:,:,4) = Corr(:,:,pos4);
    end
    midx = nan(size(idx));
    midx(idx == pos1) = -1;
    midx(idx == pos2) = -2;
    midx(idx == pos3) = -3;
    midx(idx == pos4) = -4;
    
    midx(midx == -1) = 1;
    midx(midx == -2) = 2;
    midx(midx == -3) = 3;
    midx(midx == -4) = 4;
end

end