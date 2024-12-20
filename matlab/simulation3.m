clear
clc
close all


%% Sim 3
casee = 3;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
N = 100;                                                 % number of repetition (realizations)
SNR = -3;

delphi = 2*pi./(1+exp(-0.01*(t-334)));

x = cos(2*pi*f*t);                                          % first signal  (pure tone - ie monocomponent)
y = cos(2*pi*f*t + delphi) + cos(2*pi*f*1.1*t + delphi);    % second signal (multicomponent)
pwelch_window = length(x);

if ~exist(fullfile("simulated","sim3.mat"))
    for m = 1:N
        noise = mvnrnd([0 0],[1 0;0 1],length(t))';
        ex = noise(1,:);
        ey = noise(2,:);
        snr_x = snr(x,ex); snr_y = snr(y,ey);
        ex = ex*sqrt(10^((snr_x-SNR)/10));
        ey = ey*sqrt(10^((snr_y-SNR)/10));
        XN = x + ex;
        YN = y + ey;
        Data{m} = [XN;YN]';
        Data{m} = Data{m}./max(abs(Data{m}));
    end
    save(fullfile("simulated","sim3.mat"),"Data")
else
    load(fullfile("simulated","sim3.mat"))
end

minPeaks = zeros(1,N);
compStds = zeros(1,N);
winds = zeros(1,N);

if ~exist(fullfile("decomposed","sim3_MVMD.mat"))    
    %% na-MEMD
    numModes = zeros(1,100);
    for m = 1:N
        stp_crit = 'stop'; mode = 'na_fix';
        stp_vec = [0.3,0.3,0.3];
        n_channel_na = size(Data{m},2);  
        ndir = 8*n_channel_na;
        intensity_noise = 0.75;
        imf = namemd(Data{m}, ndir, stp_crit, stp_vec, mode, intensity_noise, n_channel_na);
        numModes(m) = size(imf{1},1);
    end

    %% MVMD
    imfs_MVMD = cell(1,N);
    T_MVMD = zeros(1,N);
    for m = 1:N
        tau = 0.01; DC = 1; init = 0; tol = 1e-9;
        K = numModes(m); alpha = 1000;
        tic
        imf = MVMD(Data{m},alpha,tau,K,DC,init,tol);
        T_MVMD(m) = toc;
        imfs_MVMD{m} = imf;
    end
    disp('MVMD')
    save(fullfile("decomposed","sim3_MVMD.mat"),"imfs_MVMD","T_MVMD")
end


%% MSWD-CL
imfs_MSWD = cell(1,N);
corrs = cell(1,N);
T_MSWD = zeros(1,N);
if ~exist(fullfile("decomposed","sim3_MSWD.mat"))
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
    save(fullfile("decomposed","sim3_MSWD.mat"),"imfs_MSWD","T_MSWD")
end
