% This function performs the time-varying phase synchronization analysis 
% for the HCP data

function COSDELPHI1 = phase_sync_analysis_HCP(imf,casee,indx,f)

if casee == "MSWD"
    COSDELPHI1 = nan(size(imf,1),size(indx,1),length(f));
else
    COSDELPHI1 = nan(size(imf,2),size(indx,1),length(f));
end

if casee == "MSWD"
    if length(size(imf)) == 3
        imf = permute(imf,[3,1,2]);
        if length(size(imf)) == 3
            imf = imf(1:end-1,:,:);
            positions = f;
            for f1 = 1:length(f)
                pos_comp = positions(f1);
                % Phase Synchronization analysis:
                for j = 1:size(indx,1)
                    for i = 1:size(imf,1)
                        if i == pos_comp
                            dat{j,i} = [imf(i,:,indx(j,1));imf(i,:,indx(j,2))];
                            H = hilbert(dat{j,i}');
                            sigphase = angle(H);
                            DELPHI{j} = sigphase(:,1)-sigphase(:,2);
                            COSDELPHI1(:,j,f1) = cos(DELPHI{j});
                        end
                    end
                end
            end
        end
    end
elseif casee == "MVMD"
    positions = f;
    for f1 = 1:length(f)
        pos_comp = positions(f1);
        % Phase Synchronization analysis:
        for j = 1:size(indx,1)
            for i = 1:size(imf,1)
                if i == pos_comp
                    dat{j,i} = [imf(i,:,indx(j,1));imf(i,:,indx(j,2))];
                    H = hilbert(dat{j,i}');
                    sigphase = angle(H);
                    DELPHI{j} = sigphase(:,1)-sigphase(:,2);
                    COSDELPHI1(:,j,f1) = cos(DELPHI{j});
                end
            end
        end
    end
end

end