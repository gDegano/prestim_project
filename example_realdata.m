clc
clear

warning off

prestim=load('eeg_sample_to_model.mat'); % data to model (only prestim 2 sec)
complete=load('eeg_sample.mat'); % complete data (3 sec)

% Each files contains the following data:
%   DATA = 3D matrix (chan x time x trl)
%   EEG  = structure containing Fs rate and time vector

n_trl=size(prestim.DATA,3);
n_chan=size(prestim.DATA,1);
central_freq=10;
num_cycles=5;

%% Compute with model

EEGopts.pnts=size(prestim.DATA,2);
EEGopts.srate=prestim.EEG.srate;
EEGopts.times=prestim.EEG.time;
EEGopts.model=2;

freqslideFilt=zeros(size(prestim.DATA,2),n_chan,n_trl);
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(prestim.DATA(k,:,j));
        temp_freq=CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        freqslideFilt(:,k,j) = temp_freq(1:EEGopts.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end

AVG_median=mean(mean(freqslideFilt,3),2);

%% Compute with real data

EEGopts_c.pnts=size(complete.DATA,2);
EEGopts_c.srate=complete.EEG.srate;
EEGopts_c.times=complete.EEG.time;
EEGopts_c.model=0;

freqslideFilt_c=zeros(size(complete.DATA,2),n_chan,n_trl);
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(complete.DATA(k,:,j));
        temp_freq=CCN_freq_slide(data2use,EEGopts_c,central_freq,num_cycles);
        freqslideFilt_c(:,k,j) = temp_freq(1:EEGopts_c.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end

AVG_median_c=mean(mean(freqslideFilt_c,3),2);

%% Compute with zero padding

EEGopts.pnts=size(prestim.DATA,2);
EEGopts.srate=prestim.EEG.srate;
EEGopts.times=prestim.EEG.time;
EEGopts.model=0;

freqslideFilt_zp=zeros(size(prestim.DATA,2),n_chan,n_trl);
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(prestim.DATA(k,:,j));
        temp_freq=CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        freqslideFilt_zp(:,k,j) = temp_freq(1:EEGopts.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end

AVG_median_zp=mean(mean(freqslideFilt_zp,3),2);

%% Plots

time_1sec=EEGopts_c.times(EEGopts_c.times>=-1 & EEGopts_c.times<=0);
freq_slide_1sec_AR=AVG_median(end-EEGopts.srate:end);
freq_slide_1sec_zp=AVG_median_zp(end-EEGopts.srate:end);
freq_slide_1sec_real=AVG_median_c(EEGopts_c.times>=-1 & EEGopts_c.times<=0);

figure
grid on
plot(time_1sec,freq_slide_1sec_AR)
hold on
plot(time_1sec,freq_slide_1sec_zp)
plot(time_1sec,freq_slide_1sec_real)
legend({'AR model';'Zero padding';'Real Data'},'Location','Best')
xlabel('Time')
ylabel('Frequency')

clc

% Samplig points to evaluate goodness of estimation
npoints=75;

% pearson correlation
disp(['Correlation coeff for last ',num2str(round(npoints/EEGopts.srate*1000,0)),' ms: '])
disp(['AR :',num2str(corr(freq_slide_1sec_AR(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
disp(['Zero padding :',num2str(corr(freq_slide_1sec_zp(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])

% mutual info
disp(['MI for last ',num2str(round(npoints/EEGopts.srate*1000,0)),' ms: '])
disp(['AR :',num2str(gcmi_cc(freq_slide_1sec_AR(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
disp(['Zero padding :',num2str(gcmi_cc(freq_slide_1sec_zp(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])

% RMSE sqrt(ei^2) with ei = yi - y~i
e_AR=freq_slide_1sec_AR(end-npoints:end)-freq_slide_1sec_real(end-npoints:end);
e_zp=freq_slide_1sec_zp(end-npoints:end)-freq_slide_1sec_real(end-npoints:end);
disp(['RMSE for last ',num2str(round(npoints/EEGopts.srate*1000,0)),' ms: '])
disp(['AR :',num2str(sqrt(mean(e_AR.^2)))])
disp(['Zero padding :',num2str(sqrt(mean(e_zp.^2)))])

% sMAPE = mean(200|ei|/(yi+y~i))
div_AR=freq_slide_1sec_AR(end-npoints:end)+freq_slide_1sec_real(end-npoints:end);
div_zp=freq_slide_1sec_zp(end-npoints:end)+freq_slide_1sec_real(end-npoints:end);
disp(['sMAPE for last ',num2str(round(npoints/EEGopts.srate*1000,0)),' ms: '])
disp(['AR :',num2str(mean(200*abs(e_AR)./div_AR))])
disp(['Zero padding :',num2str(mean(200*abs(e_zp)./div_zp))])



