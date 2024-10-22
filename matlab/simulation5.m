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
    
    %% MVMD
    imfs_MVMD = cell(1,N);
    T_MVMD = zeros(1,N);
    if ~exist(fullfile("decomposed",['sim5_MVMD_SNR',num2str(SNR(nn)),'.mat']))
        for m = 1:N
            tau = 0; DC = 1; init = 0; tol = 1e-20;
            K = 2; alpha = 1000;
            tic
            imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
            T_MVMD(m) = toc;
            imfs_MVMD{m} = imf;
        end
        disp('MVMD')
        save(fullfile("decomposed",['sim5_MVMD_SNR',num2str(SNR(nn)),'.mat']),"imfs_MVMD","T_MVMD")
    end
   
    
    %% MSWD-CL
    imfs_MSWD_CL = cell(1,N);
    corrs = cell(1,N);
    T_MSWD_CL = zeros(1,N);
    if ~exist(fullfile("decomposed",['sim5_MSWD_SNR',num2str(SNR(nn)),'.mat']))
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
            [imf,corrs{m}] = MSWD(Data{m}, param_struct);
            T_MSWD_CL(m) = toc;
            imfs_MSWD_CL{m} = imf;
        end
        disp('MSWD')
        save(fullfile("decomposed",['sim5_MSWD_SNR',num2str(SNR(nn)),'.mat']),"imfs_MSWD_CL","corrs","T_MSWD_CL")
    end
end

