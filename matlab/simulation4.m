clear
clc
close all


%% Sim 4
casee = 4;
TR = 2;                                                                     % Repetition time
fs = 1/TR;                                                                  % Sampling frequency
t = 0:1/fs:1000-1/fs;
f = 0.05;
N = 100;                                                 % number of repetition (realizations) 
nS = 3;
SNR = -3;


phi1 = pi*((t-100>0).*(t-250<0) + (t-300>0).*(t-500<0) + (t-600>0).*(t-800<0));
phi2 = pi*((t-100>0).*(t-250<0) +                      - (t-600>0).*(t-800<0));
phi3 = pi*(                     - (t-300>0).*(t-500<0) - (t-600>0).*(t-800<0));
x = cos(2*pi*f*t + phi1);    % first signal
y = cos(2*pi*f*t + phi2);    % second signal
z = cos(2*pi*f*t + phi3);    % thrid signal
pwelch_window = length(x);

if ~exist(fullfile("simulated","sim4.mat"))
    for m = 1:N
        noise = mvnrnd([0 0 0],[1 0 0;0 1 0;0 0 1],length(t))';
        ex = noise(1,:);
        ey = noise(2,:);
        ez = noise(3,:);
        snr_x = snr(x,ex); snr_y = snr(y,ey); snr_z = snr(z,ez);
        ex = ex*sqrt(10^((snr_x-SNR)/10));
        ey = ey*sqrt(10^((snr_y-SNR)/10));
        ez = ez*sqrt(10^((snr_z-SNR)/10));
        XN = x + ex;
        YN = y + ey;
        ZN = z + ez;
        Data{m} = [XN;YN;ZN]';
        Data{m} = Data{m}./max(abs(Data{m}));
    end
    save(fullfile("simulated","sim4.mat"),"Data")
else
    load(fullfile("simulated","sim4.mat"))
end

minPeaks = zeros(1,N);
compStds = zeros(1,N);
winds = zeros(1,N);

indx = nchoosek(1:size(Data{1},2),2);

%% MVMD
imfs_MVMD = cell(1,N);
T_MVMD = zeros(1,N);
if ~exist(fullfile("decomposed","sim4_MVMD.mat"))
    for m = 1:N
        tau = 0.01; DC = 1; init = 0; tol = 1e-9;
        K = 4; alpha = 1000;
        tic
        imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
        T_MVMD(m) = toc;
        imfs_MVMD{m} = imf;
    end
    disp('MVMD')
    save(fullfile("decomposed","sim4_MVMD.mat"),"imfs_MVMD","T_MVMD")
end


%% MSWD-CL
imfs_MSWD = cell(1,N);
corrs = cell(1,N);
T_MSWD = zeros(1,N);
if ~exist(fullfile("decomposed","sim4_MSWD_new.mat"))
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
    save(fullfile("decomposed","sim4_MSWD_new.mat"),"imfs_MSWD","T_MSWD")
end