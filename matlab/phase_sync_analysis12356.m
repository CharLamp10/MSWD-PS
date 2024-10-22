% This function performs the phase synchronization analysis (static and
% time-varying) for simulations 1,2,3,4 and 6

function [COSDELPHI1,plv,mfreq] = phase_sync_analysis12356(imf,...
    casee,fs,min_freq,length_plv,f)

if casee == "MSWD"
    if length(size(squeeze(imf))) == 3
        imf = permute(imf,[3,1,2]);
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
            if length(f) == 1
                [~,positions] = min(abs(mfreq - f));
            else
                positions = find_freq_pos(mfreq,f);
            end
    %       Phase Synchronization analysis:
            c = 1;
            for p = 1:length(positions) 
                pos_comp = positions(p);
                for j = 1:size(imf,1)
                    dat{j} = [imf(j,:,1);imf(j,:,2)];
                    H = hilbert(dat{j}');
                    sigphase = angle(H);
                    if j == pos_comp
                        if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                            DELPHI = sigphase(:,1)-sigphase(:,2);
                            COSDELPHI1(:,c) = cos(DELPHI);
                        else
                            COSDELPHI1(:,c) = nan(1,size(imf,2));
                        end
                        c = c + 1;
                    end
                    if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                        dPhi = sigphase(:,1) - sigphase(:,2);
                        e = exp(1i*dPhi);
                        plv(j) = abs(mean(e));
                    else
                        plv(j) = 0;
                    end
                end
            end
        else
            imf = squeeze(imf);
            plv = zeros(1,length_plv);
            COSDELPHI1 = nan(1,size(imf,1));%zeros(1,size(imf,1));
            mfreq = nan;
            COSDELPHI1 = [];
        end
    else
        imf = squeeze(imf);
        plv = zeros(1,length_plv);
        COSDELPHI1 = nan(1,size(imf,1));%zeros(1,size(imf,1));
        mfreq = nan;
        COSDELPHI1 = [];
    end
elseif casee == "MVMD"
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
    if length(f) == 1
        [~,positions] = min(abs(mfreq - f));
    else
        positions = find_freq_pos(mfreq,f);
    end
    %Synchronization analysis:
    if length(size(imf)) == 3
        c = 1;
        for p = 1:length(positions)
            pos_comp = positions(p);
            for j = 1:size(imf,1)
                dat{j} = [imf(j,:,1);imf(j,:,2)];
                H = hilbert(dat{j}');
                sigphase = angle(H);
                if j == pos_comp
                    if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                        DELPHI = sigphase(:,1)-sigphase(:,2);
                        COSDELPHI1(:,c) = cos(DELPHI);
                    else
                        COSDELPHI1(:,c) = nan(1,size(imf,2));
                    end
                    c = c + 1;
                end
                if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
                    dPhi = sigphase(:,1) - sigphase(:,2);
                    e = exp(1i*dPhi);
                    plv(j) = abs(mean(e));
                else
                    plv(j) = 0;
                end
            end
        end
    else
        dat{1} = [imf(:,1),imf(:,2)];
        H = hilbert(dat{1});
        sigphase = angle(H);
        if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
            DELPHI = sigphase(:,1)-sigphase(:,2);
            COSDELPHI1 = cos(DELPHI);
        else
            COSDELPHI1 = nan(1,size(imf,1));
        end
        if ~all(sigphase(:,1) == 0) && ~all(sigphase(:,2) == 0)
            dPhi = sigphase(:,1) - sigphase(:,2);
            e = exp(1i*dPhi);
            plv = [abs(mean(e)),zeros(1,length_plv-1)];
        else
            plv = [zeros(1,length_plv)];
        end
    end
end

end