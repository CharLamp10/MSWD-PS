clear
clc
close all


%% Sim 6
casee = 6;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f1 = 0.05;
f2 = 0.03;
N = 100;                                                 % number of repetition (realizations)
SNR = 0;

delphi1 = 2*pi./(1+exp(-0.01*(t-334)));
delphi2 = pi/4;
pwelch_window = length(delphi1);

if ~exist(fullfile("simulated","sim6.mat"))
    for m = 1:N
        a1 = 0.75 + (1-0.75)*rand(1); 
        a2 = 0.75 + (1-0.75)*rand(1);
        a3 = 0.75 + (1-0.75)*rand(1);
        a4 = 0.75 + (1-0.75)*rand(1);
        a5 = 0.75 + (1-0.75)*rand(1);
        a6 = 0.75 + (1-0.75)*rand(1);
        x{m} = a1*cos(2*pi*f2*t) + a2*cos(2*pi*f1*t) + a3*cos(2*pi*1.3*f1*t);
        y{m} = a4*cos(2*pi*f2*t + delphi2) + a5*cos(2*pi*f1*t + delphi1) + a6*cos(2*pi*0.3*f1*t);
    end
    
    for m = 1:N
        noise = mvnrnd([0 0],[1 0;0 1],length(t))';
        ex = noise(1,:);
        ey = noise(2,:);
        snr_x = snr(x{m},ex); snr_y = snr(y{m},ey);
        ex = ex*sqrt(10^((snr_x-SNR)/10));
        ey = ey*sqrt(10^((snr_y-SNR)/10));
        XN = x{m} + ex;
        YN = y{m} + ey;
        Data{m} = [XN;YN]';
        Data{m} = Data{m}./max(abs(Data{m}));
    end
    save(fullfile("simulated","sim6.mat"),"Data")
else
    load(fullfile("simulated","sim6.mat"))
end

minPeaks = zeros(1,N);
compStds = zeros(1,N);
winds = zeros(1,N);

%% NA-MEMD
imfs_MEMD = cell(1,N);
T_MEMD = zeros(1,N);
if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MEMD.mat']))
    for m = 1:N
        stp_crit = 'stop';
        stp_vec = [0.3 0.3 0.3];
        mode = 'na_snr';
        intensity_noise = 0.75; 
        n_channel_na = size(Data{m},2);  
        ndir = 8*n_channel_na; 
        tic
        imfs = namemd(Data{m}, ndir, stp_crit, stp_vec, mode, intensity_noise, n_channel_na);
        for i = 1:length(imfs)
            imf(:,:,i) = imfs{i};
        end
        T_MEMD(m) = toc;
        imfs_MEMD{m} = imf;
        clear imf
    end
    save(fullfile("decomposed",['sim',num2str(casee),'_MEMD.mat']),"imfs_MEMD","T_MEMD")
    disp('MEMD')
end

%% EWT
imfs_EWT = cell(1,N);
T_EWT = zeros(1,N);
if ~exist(fullfile("decomposed",['sim',num2str(casee),'_EWT.mat']))
    for m = 1:N
        tic
        for i = 1:size(Data{m},2)
            temp_imf = ewt(Data{m}(:,i))';
            mfreqs = zeros(1,size(temp_imf,1));
            for j = 1:size(temp_imf,1)
                mfreqs(j) = meanfreq(temp_imf(j,:),fs);
            end
            [val1,pos1] = min(abs(mfreqs - f1));
            [val2,pos2] = min(abs(mfreqs - f2));
            if pos1 ~= pos2
                temp1 = temp_imf(pos1,:);
                temp2 = temp_imf(pos2,:);
            else
                if val1 < val2
                    temp1 = temp_imf(pos1,:);
                    mfreqs(pos1) = 1e10;
                    [~,pos2] = min(abs(mfreqs - f2));
                    temp2 = temp_imf(pos2,:);
                else
                    temp2 = temp_imf(pos2,:);
                    mfreqs(pos2) = 1e10;
                    [~,pos1] = min(abs(mfreqs - f2));
                    temp1 = temp_imf(pos1,:);
                end
            end
            temp_imf(pos1,:) = temp_imf(1,:);
            temp_imf(pos2,:) = temp_imf(2,:);
            temp_imf(1,:) = temp1;
            temp_imf(2,:) = temp2;
            imf(1:size(temp_imf,1),1:size(temp_imf,2),i) = temp_imf;
            clear temp_imf
        end
        T_EWT(m) = toc;
        imfs_EWT{m} = imf;
        clear imf
    end
    save(fullfile("decomposed",['sim',num2str(casee),'_EWT.mat']),"imfs_EWT","T_EWT")
    disp('EWT')
end

%% MVMD
imfs_MVMD = cell(1,N);
T_MVMD = zeros(1,N);
if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MVMD.mat']))
    for m = 1:N
        tau = 0; DC = 1; init = 0; tol = 1e-9;
        K = 8; alpha = 1000; %K = 2 didn't work well for MVMD
        tic
        imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
        T_MVMD(m) = toc;
        imfs_MVMD{m} = imf;
    end
    disp('MVMD')
    save(fullfile("decomposed",['sim',num2str(casee),'_MVMD.mat']),"imfs_MVMD")
end


%% MSWD
imfs_MSWD = cell(1,N);
corrs = cell(1,N);
T_MSWD = zeros(1,N);
if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MSWD.mat']))
    for m = 1:N
        p_value = 1e-5;
        wind = 1;
        compStd = 0.002;
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
    save(fullfile("decomposed",['sim',num2str(casee),'_MSWD.mat']),"imfs_MSWD","T_MSWD")
end
