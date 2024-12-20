clear
clc
close all

%% Sim 8
casee = 8;
TR = 2;         % Repetition time
fs = 1/TR;      % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
N = 100;        % number of repetition (realizations)
SNR = -3;

x1 = cos(2*pi*f*t);
y1 = cos(2*pi*f*t - pi/4);
z1 = cos(2*pi*f*t - pi/2);
x2 = cos(2*pi*f*t + pi);
y2 = cos(2*pi*f*t + 3*pi/4);
z2 = cos(2*pi*f*t + pi/2);

pwelch_window = length(x1);

if ~exist(fullfile("simulated","sim8.mat"))
    for m = 1:N
        noise = mvnrnd([0 0 0,0,0,0],eye(6),length(t))';
        ex1 = noise(1,:); ex2 = noise(4,:);
        ey1 = noise(2,:); ey2 = noise(5,:);
        ez1 = noise(3,:); ez2 = noise(6,:);
        snr_x1 = snr(x1,ex1); snr_x2 = snr(x2,ex2);
        snr_y1 = snr(y1,ey1); snr_y2 = snr(y2,ey2);
        snr_z1 = snr(z1,ez1); snr_z2 = snr(z2,ez2);
        ex1 = ex1*sqrt(10^((snr_x1-SNR)/10));
        ex2 = ex2*sqrt(10^((snr_x2-SNR)/10));
        ey1 = ey1*sqrt(10^((snr_y1-SNR)/10));
        ey2 = ey2*sqrt(10^((snr_y2-SNR)/10));
        ez1 = ez1*sqrt(10^((snr_z1-SNR)/10));
        ez2 = ez2*sqrt(10^((snr_z2-SNR)/10));
        X1N = x1 + ex1; X2N = x2 + ex2;
        Y1N = y1 + ey1; Y2N = y2 + ey2;
        Z1N = z1 + ez1; Z2N = z2 + ez2;
        Data{m} = [X1N;Y1N;Z1N;X2N;Y2N;Z2N]';
        Data{m} = Data{m}./max(abs(Data{m}));
    end
    save(fullfile("simulated","sim8.mat"),"Data")
else
    load(fullfile("simulated","sim8.mat"))
end

indx = nchoosek(1:size(Data{1},2),2);
minPeaks = zeros(1,N);
compStds = zeros(1,N);
winds = zeros(1,N);

%% MVMD
imfs_MVMD = cell(1,N);
T_MVMD = zeros(1,N);
if ~exist(fullfile("decomposed","sim8_MVMD.mat"))
    for m = 1:N
        tau = 0; DC = 1; init = 0; tol = 1e-9;
        K = 2; alpha = 1000;
        tic
        imf = MVMD_new(Data{m},alpha,tau,K,DC,init,tol);
        T_MVMD(m) = toc;
        imfs_MVMD{m} = imf;
    end
    disp('MVMD')
    save(fullfile("decomposed","sim8_MVMD.mat"),"imfs_MVMD","T_MVMD")
end


%% MSWD-CL
imfs_MSWD = cell(1,N);
corrs = cell(1,N);
T_MSWD = zeros(1,N);
if ~exist(fullfile("decomposed","sim8_MSWD.mat"))
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
    save(fullfile("decomposed","sim8_MSWD.mat"),"imfs_MSWD","T_MSWD")
end
