% This function performs the phase synchronization analysis (static and
% time-varying) for simulations 4 and 7
function [COSDELPHI1,plv,mfreq] = phase_sync_analysis47(imf,casee,indx,fs,f,min_freq,length_plv)

if casee == "MSWD"
    if length(size(squeeze(imf))) == 3
        imf = permute(imf,[3,1,2]);
        % imf = zero_imfs(imf,thresh_same_imf);
        if length(size(imf)) == 3
            imf = imf(1:end-1,:,:);
            mfreq_temp = nan(size(imf,1),size(imf,3));
            for i = 1:size(imf,3)
                for j = 1:size(imf,1)
                    if ~all(squeeze(imf(j,:,i)) == 0)
                        mfreq_temp(j,i) = meanfreq(imf(j,:,i),fs);  
                    end
                end
            end
            mfreq_temp = mean(mfreq_temp,2,'omitnan');
            indeces = find(mfreq_temp > min_freq);
            imf = imf(indeces,:,:); mfreq_temp = mfreq_temp(indeces);
            mfreq = mfreq_temp;
            [~,pos_comp] = min(abs(mfreq - f));
            % Phase Synchronization analysis:
            for j = 1:size(indx,1)
                for i = 1:size(imf,1)
                    dat{j,i} = [imf(i,:,indx(j,1));imf(i,:,indx(j,2))];
                    H = hilbert(dat{j,i}');
                    sigphase = angle(H);
                    if i == pos_comp
                        if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                            DELPHI{j} = sigphase(:,1)-sigphase(:,2);
                            COSDELPHI1(:,j) = cos(DELPHI{j});
                        else
                            COSDELPHI1(:,j) = nan(1,size(imf,2));
                        end
                    end
                    if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                        dPhi = sigphase(:,1) - sigphase(:,2);
                        e = exp(1i*dPhi);
                        plv(i,j) = abs(mean(e));
                    else
                        plv(i,j) = 0;
                    end
                end
            end
        else
            imf = squeeze(imf);
            plv = zeros(length_plv);
            mfreq = nan;
            for j = 1:size(indx,1)
                COSDELPHI1(:,j) = nan(1,size(imf,1));%zeros(1,size(imf,1));
            end
        end
    else
        imf = squeeze(imf);
        plv = zeros(length_plv);
        mfreq = nan;
        for j = 1:size(indx,1)
            COSDELPHI1(:,j) = nan(1,size(imf,1));%zeros(1,size(imf,1));
        end
    end
elseif casee == "MVMD"
    % imf = zero_imfs(imf,thresh_same_imf);
    imf = squeeze(imf);
    if length(size(imf)) == 3
        mfreq_temp = nan(size(imf,1),size(imf,3));
    else
        mfreq_temp = nan(1,size(imf,2));
    end
    if length(size(imf)) == 3
        for i = 1:size(imf,3)
            for j = 1:size(imf,1)
                if ~all(squeeze(imf(j,:,i)) == 0)
                    mfreq_temp(j,i) = meanfreq(imf(j,:,i),fs);  
                end
            end
        end
    else
        for i = 1:length(mfreq_temp)
            if ~all(imf(:,i) == 0)
                mfreq_temp(i) = meanfreq(imf(:,i),fs);  
            end
        end
    end
    mfreq_temp = mean(mfreq_temp,2,'omitnan');
    if length(size(imf)) == 3
        indeces = find(mfreq_temp > min_freq);
        imf = imf(indeces,:,:); mfreq_temp = mfreq_temp(indeces);
    end
    mfreq = mfreq_temp;
    [~,pos_comp] = min(abs(mfreq - f));
    % Phase Synchronization analysis:
    if length(size(imf)) == 3
        for j = 1:size(indx,1)
            for i = 1:size(imf,1)
                dat{j,i} = [imf(i,:,indx(j,1));imf(i,:,indx(j,2))];
                H = hilbert(dat{j,i}');
                sigphase = angle(H);
                if i == pos_comp
                    if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                        DELPHI{j} = sigphase(:,1)-sigphase(:,2);
                        COSDELPHI1(:,j) = cos(DELPHI{j});
                    else
                        COSDELPHI1(:,j) = nan(1,size(imf,2));
                    end
                end
                if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                    dPhi = sigphase(:,1) - sigphase(:,2);
                    e = exp(1i*dPhi);
                    plv(i,j) = abs(mean(e));
                else
                    plv(i,j) = 0;
                end
            end
        end
    else
        for j = 1:size(indx,1)
            dat{j,1} = [imf(:,indx(j,1)),imf(:,indx(j,2))];
            H = hilbert(dat{j,1});
            sigphase = angle(H);
            if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                DELPHI{j} = sigphase(:,1)-sigphase(:,2);
                COSDELPHI1(:,j) = cos(DELPHI{j});
            else
                COSDELPHI1(:,j) = nan(1,size(imf,1));
            end
            if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                dPhi = sigphase(:,1) - sigphase(:,2);
                e = exp(1i*dPhi);
                plv(:,j) = [abs(mean(e)),zeros(1,size(length_plv,1)-1)];
            else
                plv(:,j) = [zeros(1,size(length_plv,1))];
            end
        end
    end
end

end