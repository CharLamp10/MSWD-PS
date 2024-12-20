clear
clc
close all


%% Sim 5
save_plots = 1;
casee = 5;
TR = 2;                                                   % Repetition time 
fs = 1/TR;                                                % Sampling frequency
t = 0:1/fs:668-1/fs;
f = 0.05;
freq_ranges = 0.01:0.02:0.1;
min_freq = 0;
N = 100;                                                 % number of repetition (realizations)                                        % window sizes for the Windowed Phase Sync. Measures
SNR = [-6,0,3,6];

true_plv = zeros(1,length(freq_ranges));
delphi = 2*pi./(1+exp(-0.01*(t-334)));
x = cos(2*pi*f*t);                                          % first signal
y = cos(2*pi*f*t + delphi);                                % second signal
pwelch_window = length(x);
phase_x = angle(hilbert(x)); phase_y = angle(hilbert(y));
dPhi = phase_x - phase_y;
e = exp(1i*dPhi);
pos = find(freq_ranges == f);
true_plv(pos) = abs(mean(e));
true_cosdelphi = cos(phase_x - phase_y);
filename = fullfile(pwd,'plots',['sim',num2str(casee) , '_ground_truth.png']);
simulation2_groundTruthPlot(true_cosdelphi,t,filename)


for nn = 1:length(SNR)  
    load(fullfile("decomposed",['sim5_MVMD_SNR',num2str(SNR(nn)),'.mat']))
    load(fullfile("decomposed",['sim5_MSWD_SNR',num2str(SNR(nn)),'.mat']))
    %% MVMD
    for m = 1:N
        imf = imfs_MVMD{m};
        if m == 1
            plv_length = length(freq_ranges);
        else
            plv_length = length(plv(m-1,:));
        end
        [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
        "MVMD",fs,min_freq,plv_length,f);
        COSDELPHI1(:,m) = COSDELPHI11;
        mfreq{m} = mfreq1;
        plv(m,1:length(plv1)) = plv1;
    end
    [plv_MVMD,~] = sort_values(mfreq,freq_ranges,plv);
    [mean_plv_MVMD,lower_plv_MVMD,upper_plv_MVMD] = confidence_interval(plv_MVMD);
    
    for m = 1:N
        cosdelphi_mvmd(m,:) = COSDELPHI1(:,m);
    end
    disp('MVMD')
    clear COSDELPHI1 CCORSW mfreq plv
    
        
    %% MSWD
    for m = 1:N
        imf = imfs_MSWD{m};
        if m == 1
            plv_length = length(freq_ranges);
        else
            plv_length = length(plv(m-1,:));
        end
        [COSDELPHI11,plv1,mfreq1] = phase_sync_analysis12356(imf,...
        "MSWD",fs,min_freq,plv_length,f);
        COSDELPHI1(:,m) = COSDELPHI11;
        mfreq{m} = mfreq1;
        plv(m,1:length(plv1)) = plv1;
    end
    [plv_MSWD,~] = sort_values(mfreq,freq_ranges,plv);
    [mean_plv_MSWD,lower_plv_MSWD,upper_plv_MSWD] = confidence_interval(plv_MSWD);
    
    %% MSWD PLOTS
    for m = 1:N
        freqs = mfreq{m};
        cosdelphi_mswd(m,:) = COSDELPHI1(:,m);
    end
    
    disp('MSWD')
    clear COSDELPHI1 CCORSW mfreq plv

    filename = fullfile(pwd,'plots',['sim', num2str(casee), '_noisevar', num2str(SNR(nn)), '_TVPS.png']);
    if save_plots
        simulation12356Plots(true_cosdelphi,cosdelphi_mvmd,cosdelphi_mswd,t,filename)
    else
        simulation12356Plots(true_cosdelphi,cosdelphi_mvmd,cosdelphi_mswd,t)
    end
    
    E_plv_MVMD = mean_plv_MVMD - lower_plv_MVMD;
    E_plv_MSWD = mean_plv_MSWD - lower_plv_MSWD;
    
    
    plv_min = min([min(lower_plv_MVMD),...
        min(lower_plv_MSWD),min(true_plv)]);
    plv_max = max([max(upper_plv_MVMD),...
        max(upper_plv_MSWD),max(true_plv)]);
    
    figure
    errorbar(1:size(plv_MVMD,2),mean_plv_MVMD,E_plv_MVMD,'bo',"LineStyle","none",'LineWidth',1)
    hold on
    errorbar(size(plv_MVMD,2)+1:size(plv_MVMD,2)+size(plv_MSWD,2),mean_plv_MSWD,E_plv_MSWD,'go',"LineStyle","none",'LineWidth',1)
    plot([true_plv,true_plv],'ro','LineWidth',1)
%     title('Simulation 6: PLV')
%     xlabel('Method')
    ylabel('PLV','FontSize',16)
    ylim([0,plv_max+0.02])
    xlim([1,size(plv_MVMD,2)+size(plv_MSWD,2)])
    start = (1 + length(freq_ranges))/2;
    stop = 2*length(freq_ranges) + 1 - start;
    ticks = start:length(freq_ranges):stop;
    xticks(ticks)
    xticklabels({'MVMD-PS','MSwD-PS'})
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',16)
    xline((length(freq_ranges):length(freq_ranges):length(freq_ranges))+0.5,'--')
    hold off
    if save_plots
        exportgraphics(gcf,fullfile(pwd,'plots',['sim',num2str(casee), '_noisevar', num2str(SNR(nn)),'_plv.png']),'Resolution',1000)
    end
    
    clear plv_MVMD plv_MSWD... 
    ccorsw1_mvmd ccorsw1_mswd...
    ccorsw2_mvmd ccorsw2_mswd...
    ccorsw3_mvmd ccorsw3_mswd...
    cosdelphi_mvmd cosdelphi_mswd...
    minPeaks compStds winds
    clear Data
end

