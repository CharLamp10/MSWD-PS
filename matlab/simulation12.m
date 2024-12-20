clear
clc
close all


%% Sim 1&2
casee = 1;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
N = 100;                                                 % number of repetition (realizations)                                       % window sizes for the Windowed Phase Sync. Measures
SNR = -3;


if casee == 1
    delphi = 4*pi/334.*(t-334).*(t-334>=0);
elseif casee == 2
    delphi = 2*pi./(1+exp(-0.01*(t-334)));
end
pwelch_window = length(delphi);

x = cos(2*pi*f*t);                                          % first signal
y = cos(2*pi*f*t + delphi);                                % second signal

if ~exist("decomposed")
    mkdir("decomposed")
end
if ~exist("plots")
    mkdir("plots")
end
if ~exist("simulated")
    mkdir("simulated")
end

if ~exist(fullfile("simulated",['sim',num2str(casee),'.mat']))
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
    save(fullfile("simulated",['sim',num2str(casee),'.mat']),"Data")
else
    load(fullfile("simulated",['sim',num2str(casee),'.mat']))
end

minPeaks = zeros(1,N);
compStds = zeros(1,N);
winds = zeros(1,N);

if ~exist("decomposed")
    mkdir("decomposed")
end
if ~exist("plots")
    mkdir("plots")
end


%% MVMD
imfs_MVMD = cell(1,N);
T_MVMD = zeros(1,N);
if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MVMD.mat']))
    for m = 1:N
        tau = 0; DC = 1; init = 0; tol = 1e-20;
        K = 2; alpha = 1000;
        tic
        imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
        T_MVMD(m) = toc;
        imfs_MVMD{m} = imf;
    end
    save(fullfile("decomposed",['sim',num2str(casee),'_MVMD.mat']),"imfs_MVMD","T_MVMD")
    disp('MVMD')
end


%% MSWD-CL
imfs_MSWD = cell(1,N);
corrs = cell(1,N);
T_MSWD = zeros(1,N);
if ~exist(fullfile("decomposed",['sim',num2str(casee),'_MSWD.mat']))
    for m = 1:N
        p_value = 1e-5;
        wind = 1;
        compStd = 0.02;
        L = round(length(Data{m}(:,1))/wind); welch_window  = round(L);
        freq_interval = 1 / L;
        param_struct  = struct('P_corr', 1, ...%0.95, ...
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