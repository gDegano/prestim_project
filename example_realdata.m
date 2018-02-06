clc
clear

addpath('mutual_information_toolbox');

warning off

prestim=load('eeg_sample_to_model.mat'); % data to model (only prestim 2 sec)
complete=load('eeg_sample.mat'); % complete data (3 sec)

% Each files contains the following data:
%   DATA = 3D matrix (chan x time x trl)
%   EEG  = structure containing Fs rate and time vector

%% Debug settings
prestim.DATA  = prestim.DATA(1,:,:);
complete.DATA = complete.DATA(1,:,:);


n_trl=size(prestim.DATA,3);
n_chan=size(prestim.DATA,1);
central_freq=10;
num_cycles=5;


%% Compute with real data

EEGopts.pnts       = size(complete.DATA,2);
EEGopts.srate      = complete.EEG.srate;
EEGopts.times      = complete.EEG.time;
EEGopts.model      = []; % 1 = AR; 2 = ARMA
EEGopts.method     = 'hilbert';
EEGopts.plotfilt   = false;
EEGopts.filtwin    = [6 14];
EEGopts.transwidth = 0.15;
EEGopts.figures    = false;

freqslideFilt_c=zeros(size(complete.DATA,2),n_chan,n_trl);
dataAll_c = zeros(n_trl, n_chan, size(complete.DATA,2));
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(complete.DATA(k,:,j));
        [temp_freq, data_out] = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        dataAll_c(j,k,:) = data_out;
        freqslideFilt_c(:,k,j) = temp_freq(1:EEGopts.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end
AVG_median_c=mean(mean(freqslideFilt_c,3),2);



%% Compute with zero padding

EEGopts.pnts       = size(complete.DATA,2);
EEGopts.srate      = complete.EEG.srate;
EEGopts.times      = complete.EEG.time;
EEGopts.model      = []; % 1 = AR; 2 = ARMA
EEGopts.method     = 'hilbert';
EEGopts.plotfilt   = false;
EEGopts.filtwin    = [6 14];
EEGopts.transwidth = 0.15;
EEGopts.figures    = false;

freqslideFilt_zp=zeros(size(complete.DATA,2),n_chan,n_trl);
data2use = zeros(1,size(complete.DATA,2));
dataAll_zp = zeros(n_trl, n_chan, size(data2use,2));
for j=1:n_trl
    for k=1:n_chan
        data2use(1,1:size(prestim.DATA,2)) = squeeze(prestim.DATA(k,:,j));
        [temp_freq, data_out] = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        dataAll_zp(j,k,:) = data_out;
        freqslideFilt_zp(:,k,j) = temp_freq(1:EEGopts.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end
AVG_median_zp=mean(mean(freqslideFilt_zp,3),2);



%% Compute with mirror padding

EEGopts.pnts       = size(complete.DATA,2);
EEGopts.srate      = complete.EEG.srate;
EEGopts.times      = complete.EEG.time;
EEGopts.model      = []; % 1 = AR; 2 = ARMA
EEGopts.method     = 'hilbert';
EEGopts.plotfilt   = false;
EEGopts.filtwin    = [6 14];
EEGopts.transwidth = 0.15;
EEGopts.figures    = false;

freqslideFilt_mp=zeros(size(complete.DATA,2),n_chan,n_trl);
data2use = zeros(1,size(complete.DATA,2));
dataAll_mp = zeros(n_trl, n_chan, size(data2use,2));
for j=1:n_trl
    for k=1:n_chan
        pre_stim_temp = squeeze(prestim.DATA(k,:,j));
        data2use(1,1:size(prestim.DATA,2)) = pre_stim_temp;
        data2use(1,size(prestim.DATA,2)+1:end) = pre_stim_temp(1,end-1:-1:size(pre_stim_temp,2)-size(data2use(1,size(prestim.DATA,2)+1:end),2));
        [temp_freq, data_out] = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        dataAll_mp(j,k,:) = data_out;
        freqslideFilt_mp(:,k,j) = temp_freq(1:EEGopts.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end
AVG_median_mp=mean(mean(freqslideFilt_mp,3),2);



%% Compute with ARMA model

EEGopts_model.pnts       = size(prestim.DATA,2);
EEGopts_model.srate      = prestim.EEG.srate;
EEGopts_model.times      = prestim.EEG.time;
EEGopts_model.model      = 2; % 1 = AR; 2 = ARMA
EEGopts_model.method     = 'hilbert';
EEGopts_model.plotfilt   = false;
EEGopts_model.filtwin    = [6 14];
EEGopts_model.transwidth = 0.15;
EEGopts_model.figures    = false;

freqslideFilt=zeros(size(complete.DATA,2),n_chan,n_trl);
dataAll_arma = zeros(n_trl, n_chan, size(complete.DATA,2));
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(prestim.DATA(k,:,j));
        [temp_freq, data_out] = CCN_freq_slide(data2use,EEGopts_model,central_freq,num_cycles);
        dataAll_arma(j,k,1:length(data_out)) = data_out;
        freqslideFilt(1:length(data_out),k,j) = temp_freq;
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end
AVG_median_arma=mean(mean(freqslideFilt,3),2);



%% Compute with AR model

EEGopts_model.pnts       = size(prestim.DATA,2);
EEGopts_model.srate      = prestim.EEG.srate;
EEGopts_model.times      = prestim.EEG.time;
EEGopts_model.model      = 1; % 1 = AR; 2 = ARMA
EEGopts_model.method     = 'hilbert';
EEGopts_model.plotfilt   = false;
EEGopts_model.filtwin    = [6 14];
EEGopts_model.transwidth = 0.15;
EEGopts_model.figures    = false;

freqslideFilt=zeros(size(complete.DATA,2),n_chan,n_trl);
dataAll_ar = zeros(n_trl, n_chan, size(complete.DATA,2));
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(prestim.DATA(k,:,j));
        [temp_freq, data_out] = CCN_freq_slide(data2use,EEGopts_model,central_freq,num_cycles);
        dataAll_ar(j,k,1:length(data_out)) = data_out;
        freqslideFilt(1:length(data_out),k,j) = temp_freq;
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end
AVG_median_ar=mean(mean(freqslideFilt,3),2);



%% Plots

% size of prediction window for models
pw   = length(data_out);
time = EEGopts.times;

% plot input data
figure('color', 'w', 'position', [50 50 700 500])
plot(time, mean(squeeze(dataAll_zp(:,1,:))), 'r'); hold on      % Zero-padded
plot(time, mean(squeeze(dataAll_mp(:,1,:))), 'm');              % Mirror-padded
plot(time(1:pw), mean(squeeze(dataAll_ar(:,1,1:pw))), 'g');     % AR model
plot(time(1:pw), mean(squeeze(dataAll_arma(:,1,1:pw))), 'b');   % ARMA model
plot(time, mean(squeeze(dataAll_c(:,1,:))), 'k');               % ERP
xlim([-0.5, 0.25]); 
title('Input data'); xlabel('Time'); ylabel('Amplitude');
legend({'Zero-padded', 'Mirror-padded', 'AR', 'ARMA', 'ERP'}, 'location', 'NorthWest')
plotspecs; grid on


% plot output data (FS)
figure('color', 'w', 'position', [50 50 700 500])
plot(time(1:pw),AVG_median_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),AVG_median_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),AVG_median_ar(1:pw),  'g')              % AR model
plot(time(1:pw),AVG_median_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),AVG_median_c(1:pw),   'k')              % ERP
xlim([-0.5, 0.25]); ylim([7 11]); 
title('Output data'); xlabel('Time'); ylabel('Amplitude');
legend({'Zero-padded', 'Mirror-padded', 'AR', 'ARMA', 'ERP'}, 'location', 'NorthWest')
plotspecs; grid on



%% Goodness of fit 

time_1sec=EEGopts.times(EEGopts.times>=-1 & EEGopts.times<=0);
freq_slide_1sec_AR=AVG_median_ar(pw-EEGopts.srate:pw);
freq_slide_1sec_zp=AVG_median_zp(pw-EEGopts.srate:pw);
freq_slide_1sec_real=AVG_median_c(EEGopts.times>=-1 & EEGopts.times<=0);

% Samplig points to evaluate goodness of estimation
npoints=75;

% pearson correlation
disp(['Correlation coeff for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
disp(['AR :',num2str(corr(freq_slide_1sec_AR(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
disp(['Zero padding :',num2str(corr(freq_slide_1sec_zp(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])

% mutual info
disp(['MI for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
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



