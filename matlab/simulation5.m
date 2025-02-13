clear
clc
close all


%% Sim 5
casee = 5;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
freq_ranges = 0.01:0.02:0.25;
min_freq = 0.01;
N = 100;                                                 % number of repetition (realizations)
SNR = [-6,0,3,6];

delphi = 2*pi./(1+exp(-0.01*(t-334)));
x = cos(2*pi*f*t);                                          % first signal
y = cos(2*pi*f*t + delphi);                                % second signal
pwelch_window = length(x);

if ~exist(fullfile("simulated","sim5_SNR-6.mat"))
    for nn = 1:length(SNR)
        for m = 1:N
            noise = mvnrnd([0 0],[1 0;0 1],length(t))';
            ex = noise(1,:);
            ey = noise(2,:);
            snr_x = snr(x,ex); snr_y = snr(y,ey);
            ex = ex*sqrt(10^((snr_x-SNR(nn))/10));
            ey = ey*sqrt(10^((snr_y-SNR(nn))/10));
            XN = x + ex;
            YN = y + ey;
            Data{m} = [XN;YN]';
            Data{m} = Data{m}./max(abs(Data{m}));
        end
        save(fullfile("simulated",['sim5_SNR',num2str(SNR(nn)),'.mat']),"Data")
    end
end

count = 1;
for nn = 1:length(SNR)
    minPeaks = zeros(1,N);
    compStds = zeros(1,N);
    winds = zeros(1,N);
    if exist(fullfile("simulated",['sim5_SNR',num2str(SNR(nn)),'.mat']))
        load(fullfile("simulated",['sim5_SNR',num2str(SNR(nn)),'.mat']))
    end

    %% NA-MEMD
    imfs_MEMD = cell(1,N);
    T_MEMD = zeros(1,N);
    if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MEMD_SNR',num2str(SNR(nn)),'.mat']))
        for m = 1:N
            stp_crit = 'stop';
            stp_vec = [0.3 0.3 0.3];
            mode = 'na_snr';
            intensity_noise = 0.75; 
            n_channel_na = size(Data{m},2);  
            ndir = 4*n_channel_na; 
            tic
            imfs = namemd(Data{m}, ndir, stp_crit, stp_vec, mode, intensity_noise, n_channel_na);
            for i = 1:length(imfs)
                imf(:,:,i) = imfs{i};
            end
            T_MEMD(m) = toc;
            imfs_MEMD{m} = imf;
            clear imf
        end
        save(fullfile("decomposed",['sim',num2str(casee),'_MEMD_SNR',num2str(SNR(nn)),'.mat']),"imfs_MEMD","T_MEMD")
        disp('MEMD')
    end
    
    %% EWT
    imfs_EWT = cell(1,N);
    T_EWT = zeros(1,N);
    if ~exist(fullfile("decomposed",['sim',num2str(casee),'_EWT_SNR',num2str(SNR(nn)),'.mat']))
        for m = 1:N
            tic
            for i = 1:size(Data{m},2)
                temp_imf = ewt(Data{m}(:,i))';
                mfreqs = zeros(1,size(temp_imf,1));
                for j = 1:size(temp_imf,1)
                    mfreqs(j) = meanfreq(temp_imf(j,:),fs);
                end
                [~,pos] = min(abs(mfreqs - f));
                temp = temp_imf(pos,:);
                temp_imf(pos,:) = temp_imf(1,:);
                temp_imf(1,:) = temp;
                imf(:,:,i) = temp_imf;
                clear temp_imf
            end
            T_EWT(m) = toc;
            imfs_EWT{m} = imf;
            clear imf
        end
        save(fullfile("decomposed",['sim',num2str(casee),'_EWT_SNR',num2str(SNR(nn)),'.mat']),"imfs_EWT","T_EWT")
        disp('EWT')
    end
    
    %% MVMD
    imfs_MVMD = cell(1,N);
    T_MVMD = zeros(1,N);
    if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MVMD_SNR',num2str(SNR(nn)),'.mat']))
        for m = 1:N
            tau = 0; DC = 1; init = 0; tol = 1e-20;
            K = 2; alpha = 1000;
            tic
            imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
            T_MVMD(m) = toc;
            imfs_MVMD{m} = imf;
        end
        disp('MVMD')
        save(fullfile("decomposed",['sim',num2str(casee),'_MVMD_SNR',num2str(SNR(nn)),'.mat']),"imfs_MVMD","T_MVMD")
    end
   
    
    %% MSWD
    imfs_MSWD = cell(1,N);
    corrs = cell(1,N);
    T_MSWD = zeros(1,N);
    if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MSWD_SNR',num2str(SNR(nn)),'.mat']))
        for m = 1:N
            p_value = 1e-5;
            wind = 1;
            compStd = 0.02;
            L = round(length(Data{m}(:,1))/wind); welch_window  = round(L);
            freq_interval = 1 / L;
            param_struct  = struct('P_corr', 1,... %0.95, ...
                'P_corr_imp',       0.02, ...
                'StD_th',       compStd, ...
                'Welch_window', welch_window, ...
                'Freq_int',     freq_interval, ...
                'p_value',      1e-5);
            tic
            imf = MSWD(Data{m}, param_struct);
            T_MSWD(m) = toc;
            imfs_MSWD{m} = imf;
        end
        disp('MSWD')
        save(fullfile("decomposed",['sim',num2str(casee),'_MSWD_SNR',num2str(SNR(nn)),'.mat']),"imfs_MSWD","T_MSWD")
    end
end

