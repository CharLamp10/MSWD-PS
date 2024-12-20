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

%% MVMD
imfs_MVMD = cell(1,N);
T_MVMD = zeros(1,N);
if ~exist(fullfile("decomposed","sim6_MVMD.mat"))
    for m = 1:N
        tau = 0; DC = 1; init = 0; tol = 1e-9;
        K = 8; alpha = 1000; %K = 2 didn't work well for MVMD
        tic
        imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
        T_MVMD(m) = toc;
        imfs_MVMD{m} = imf;
    end
    disp('MVMD')
    save(fullfile("decomposed","sim6_MVMD.mat"),"imfs_MVMD")
end


%% MSWD-CL
imfs_MSWD = cell(1,N);
corrs = cell(1,N);
T_MSWD = zeros(1,N);
if ~exist(fullfile("decomposed","sim6_MSWD.mat"))
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
    save(fullfile("decomposed","sim6_MSWD.mat"),"imfs_MSWD","T_MSWD")
end
