% Comput frequency characteristics for actual EEG data using different
% preprocessing methods. We assume that EEG data consists of single trials
% with a stimulus onset at t=0, however, these data do not contain any
% stimulus at t=0. 
%
% The different methods are designed to address the problem of
% post-stimulus activity bleeding into the pre-stimulus analysis. The
% different methods are
% 1.) Real EEG data (simply include the post-stimulus data in the analysis)
% 2.) Zero-padding after t=0
% 3.) Mirror-padding after t=0
% 4.) Autoregressive model (AR), forecast data after t=0
% 5.) Autoregressive moving average model (ARMA), forecast data after t=0
%
% Freqeuency characteristics to be estimated are
% i.)   Freqeuency sliding (instantaneous frequency)
% ii.)  Instantaneous phase
% iii.) Inter-trial coherence
% iv.)  Power in time frequency domain
%
% ------
% FS code adapted from public online resources by Mike X Cohen
% www.mikexcohen.com/
% 
% Feb,2018 Giulio Degano & Steffen Buergers
% CCN LAB
% University of Birmingham
% GXD606@student.bham.ac.uk 
% SXB1173@student.bham.ac.uk 

clc
clear

warning off



%% Load data

% Each files contains the following data:
%   DATA = 3D matrix (chan x time x trl)
%   EEG  = structure containing Fs rate and time vector
prestim=load('eeg_sample_to_model.mat'); % data to model (only prestim 2 sec)
complete=load('eeg_sample.mat'); % complete data (3 sec)

% Add dependencies
addpath('mutual_information_toolbox');



%% Debug settings

% Select single channel to speed up processing
prestim.DATA  = prestim.DATA(1,:,:);
complete.DATA = complete.DATA(1,:,:);


% Useful variables
n_trl=size(prestim.DATA,3);
n_chan=size(prestim.DATA,1);

% Wavelet parameters (if you want to use method='wavelet')
central_freq=10;
num_cycles=5;



%% Compute with real data

disp('Compute with real data')

EEGopts.pnts       = size(complete.DATA,2);
EEGopts.srate      = complete.EEG.srate;
EEGopts.times      = complete.EEG.time;
EEGopts.model      = []; % 1 = AR; 2 = ARMA
EEGopts.method     = 'hilbert';
EEGopts.plotfilt   = false;
EEGopts.filtwin    = [6 14];
EEGopts.transwidth = 0.15;
EEGopts.figures    = false;

freqslideFilt_c = zeros(size(complete.DATA,2),n_chan,n_trl);
[dataAll_c, dataCmplx_c] = deal(zeros(n_trl, n_chan, size(complete.DATA,2)));
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(complete.DATA(k,:,j));
        [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        dataAll_c(j,k,:) = temp_out;
        dataCmplx_c(j,k,:) = temp_cmplx;
        freqslideFilt_c(:,k,j) = temp_freq(1:EEGopts.pnts);
    end
    disp(['Computing trial #',num2str(j)])
end
fslide_c = mean(mean(freqslideFilt_c,3),2);         % Inst. frequency GA
[ avgPhase_c, itc_c ] = CCN_getITC( dataCmplx_c );  % Inter-trial coherence and average phase
power_c = squeeze(mean(abs(dataCmplx_c),1));



%% Compute with zero padding

disp('Compute with zero padding')

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
[dataAll_zp, dataCmplx_zp] = deal(zeros(n_trl, n_chan, size(complete.DATA,2)));
for j=1:n_trl
    for k=1:n_chan
        data2use(1,1:size(prestim.DATA,2)) = squeeze(prestim.DATA(k,:,j));
        [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        dataAll_zp(j,k,:) = temp_out;
        dataCmplx_zp(j,k,:) = temp_cmplx;
        freqslideFilt_zp(:,k,j) = temp_freq(1:EEGopts.pnts);
    end
    disp(['Computing trial #',num2str(j)])
end
fslide_zp = mean(mean(freqslideFilt_zp,3),2);
[ avgPhase_zp, itc_zp ] = CCN_getITC( dataCmplx_zp );
power_zp = squeeze(mean(abs(dataCmplx_zp),1));



%% Compute with mirror padding

disp('Compute with mirror padding')

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
[dataAll_mp, dataCmplx_mp] = deal(zeros(n_trl, n_chan, size(complete.DATA,2)));
for j=1:n_trl
    for k=1:n_chan
        pre_stim_temp = squeeze(prestim.DATA(k,:,j));
        data2use(1,1:size(prestim.DATA,2)) = pre_stim_temp;
        data2use(1,size(prestim.DATA,2)+1:end) = pre_stim_temp(1,end-1:-1:size(pre_stim_temp,2)-size(data2use(1,size(prestim.DATA,2)+1:end),2));
        [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        dataAll_mp(j,k,:) = temp_out;
        dataCmplx_mp(j,k,:) = temp_cmplx;
        freqslideFilt_mp(:,k,j) = temp_freq(1:EEGopts.pnts);
    end
    disp(['Computing trial #',num2str(j)])
end
fslide_mp = mean(mean(freqslideFilt_mp,3),2);
[ avgPhase_mp, itc_mp ] = CCN_getITC( dataCmplx_mp );
power_mp = squeeze(mean(abs(dataCmplx_mp),1));



%% Compute with ARMA model

disp('Compute with ARMA model')

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
[dataAll_arma, dataCmplx_arma] = deal(zeros(n_trl, n_chan, size(complete.DATA,2)));
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(prestim.DATA(k,:,j));
        [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts_model,central_freq,num_cycles);
        dataAll_arma(j,k,1:length(temp_out)) = temp_out;
        dataCmplx_arma(j,k,1:length(temp_out)) = temp_cmplx;
        freqslideFilt(1:length(temp_out),k,j) = temp_freq;
    end
    disp(['Computing trial #',num2str(j)])
end
fslide_arma = mean(mean(freqslideFilt,3),2);
[ avgPhase_arma, itc_arma ] = CCN_getITC( dataCmplx_arma );
power_arma = squeeze(mean(abs(dataCmplx_arma),1));



%% Compute with AR models

orders_AR = [30, 128, 256, 512];

EEGopts_model.pnts       = size(prestim.DATA,2);
EEGopts_model.srate      = prestim.EEG.srate;
EEGopts_model.times      = prestim.EEG.time;
EEGopts_model.model      = 1; % 1 = AR; 2 = ARMA
EEGopts_model.method     = 'hilbert';
EEGopts_model.plotfilt   = false;
EEGopts_model.filtwin    = [6 14];
EEGopts_model.transwidth = 0.15;
EEGopts_model.figures    = false;
EEGopts_model.secPH      = 0.5; % forecast 0.5 seconds

[dataAll_ar, dataCmplx_ar, freqslideFilt, fslide_ar, avgPhase_ar, itc_ar, ...
    power_ar] = deal(cell(1,length(orders_AR)));

for io = 1:length(orders_AR)
    
    fprintf('Compute with ARMA model of order %i\n', orders_AR(io))

    EEGopts_model.order      = orders_AR(io);

    freqslideFilt{io}=zeros(size(complete.DATA,2),n_chan,n_trl);
    [dataAll_ar{io}, dataCmplx_ar{io}] = deal(zeros(n_trl, n_chan, size(complete.DATA,2)));
    for j=1:n_trl
        for k=1:n_chan
            data2use = squeeze(prestim.DATA(k,:,j));
            [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts_model,central_freq,num_cycles);
            dataAll_ar{io}(j,k,1:length(temp_out)) = temp_out;
            dataCmplx_ar{io}(j,k,1:length(temp_out)) = temp_cmplx;
            freqslideFilt{io}(1:length(temp_out),k,j) = temp_freq;
        end
        disp(['Computing trial #',num2str(j)])
    end
    fslide_ar{io} = mean(mean(freqslideFilt{io},3),2);
    [ avgPhase_ar{io}, itc_ar{io} ] = CCN_getITC( dataCmplx_ar{io} );
    power_ar{io} = squeeze(mean(abs(dataCmplx_ar{io}),1));

end




%% Figures

% size of prediction window for models
pw   = length(temp_out);
time = EEGopts.times;
xl = [-0.1, 0.25];

% colors for AR models
colors_AR = zeros(length(orders_AR), 3);
colors_AR(:,2) = linspace(0.2, 1, length(orders_AR));
AR_labels = cell(1, length(orders_AR));

% plot input data (single trial)
fh(1) = figure('color', 'w', 'position', [50 50 1000 600]);
plot(time, squeeze(dataAll_zp(1,1,:)), 'r'); hold on           % Zero-padded
plot(time, squeeze(dataAll_mp(1,1,:)), 'm');                   % Mirror-padded
plot(time(1:pw), squeeze(dataAll_arma(1,1,1:pw)), 'b');        % ARMA model
plot(time, squeeze(dataAll_c(1,1,:)), 'k');                    % ERP
for io = 1:length(orders_AR)                                   % AR models
    plot(time(1:pw), squeeze(dataAll_ar{io}(1,1,1:pw)), 'color', colors_AR(io,:)); 
    AR_labels{io} = sprintf('AR (%i)', orders_AR(io));
end
legend([{'Zero-padded', 'Mirror-padded', ...
    'ARMA', 'Data'}, AR_labels], 'location', 'NorthWest')
xlim(xl); 
title('Input data'); xlabel('Time'); ylabel('Amplitude');
plotspecs; grid on

% for other figures let's ignore what happens after t=0, we are not
% interested in that!
xl = [-0.5, 0];

% plot output data (FS)
fh(2) = figure('color', 'w', 'position', [50 50 1000 600]);
plot(time(1:pw),fslide_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),fslide_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),fslide_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),fslide_c(1:pw),   'k')              % ERP
for io = 1:length(orders_AR)                        % AR models
    plot(time(1:pw), fslide_ar{io}(1:pw), 'color', colors_AR(io,:)); 
end
xlim(xl);  
title('Output data (instantaneous frequency)'); xlabel('Time'); ylabel('Frequency');
plotspecs; grid on


% plot output data (POWER)
fh(3) = figure('color', 'w', 'position', [50 50 1000 600]);
plot(time(1:pw),power_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),power_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),power_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),power_c(1:pw),   'k')              % ERP
for io = 1:length(orders_AR)                       % AR models
    plot(time(1:pw), power_ar{io}(1:pw), 'color', colors_AR(io,:)); 
end
xlim(xl); 
title('Output data (power)'); xlabel('Time'); ylabel('Power');
plotspecs; grid on


% plot output data (ITC)
fh(4) = figure('color', 'w', 'position', [50 50 1000 600]);
plot(time(1:pw),itc_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),itc_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),itc_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),itc_c(1:pw),   'k')              % ERP
for io = 1:length(orders_AR)                     % AR models
    plot(time(1:pw), itc_ar{io}(1:pw), 'color', colors_AR(io,:)); 
end
xlim(xl); 
title('Output data (inter-trial coherence)'); xlabel('Time'); ylabel('ITC');
plotspecs; grid on


% SUMMARY FIGURE
% plot input data (single trial)
fh(5) = figure('color', 'w', 'position', [50 50 1000 650]);
subplot(221)
plot(time, squeeze(dataAll_zp(1,1,:)), 'r'); hold on           % Zero-padded
plot(time, squeeze(dataAll_mp(1,1,:)), 'm');                   % Mirror-padded
plot(time(1:pw), squeeze(dataAll_arma(1,1,1:pw)), 'b');        % ARMA model
plot(time, squeeze(dataAll_c(1,1,:)), 'k');                    % ERP
for io = 1:length(orders_AR)                                   % AR models
    plot(time(1:pw), squeeze(dataAll_ar{io}(1,1,1:pw)), 'color', colors_AR(io,:)); 
    AR_labels{io} = sprintf('AR (%i)', orders_AR(io));
end               % ERP
xl = [-0.5, 0.25]; xlim(xl); 
title('Input data (single trial)'); xlabel('Time'); ylabel('Amplitude');
legend([{'Zero-padded', 'Mirror-padded', ...
    'ARMA', 'Data'}, AR_labels], 'location', 'NorthWest')
plotspecs; grid on
% plot output data (FS)
subplot(222)
plot(time(1:pw),fslide_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),fslide_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),fslide_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),fslide_c(1:pw),   'k')              % ERP
for io = 1:length(orders_AR)                        % AR models
    plot(time(1:pw), fslide_ar{io}(1:pw), 'color', colors_AR(io,:)); 
end
xl = [-0.5, 0]; xlim(xl); 
title('Output data (instantaneous frequency)'); xlabel('Time'); ylabel('Frequency');
plotspecs; grid on
% plot output data (POWER)
subplot(223)
plot(time(1:pw),power_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),power_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),power_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),power_c(1:pw),   'k')              % ERP
for io = 1:length(orders_AR)                       % AR models
    plot(time(1:pw), power_ar{io}(1:pw), 'color', colors_AR(io,:)); 
end
xlim(xl); 
title('Output data (power)'); xlabel('Time'); ylabel('Power');
plotspecs; grid on
% plot output data (ITC)
subplot(224)
plot(time(1:pw),itc_zp(1:pw),  'r'); hold on     % Zero-padded
plot(time(1:pw),itc_mp(1:pw),  'm')              % Mirror-padded 
plot(time(1:pw),itc_arma(1:pw),'b')              % ARMA model
plot(time(1:pw),itc_c(1:pw),   'k')              % ERP
for io = 1:length(orders_AR)                     % AR models
    plot(time(1:pw), itc_ar{io}(1:pw), 'color', colors_AR(io,:)); 
end
xlim(xl); 
title('Output data (inter-trial coherence)'); xlabel('Time'); ylabel('ITC');
plotspecs; grid on




%% Goodness of fit 

% Needs to take into account multiple cells of AR models...


% time_1sec=EEGopts.times(EEGopts.times>=-1 & EEGopts.times<=0);
% freq_slide_1sec_AR=fslide_ar(pw-EEGopts.srate:pw);
% freq_slide_1sec_zp=fslide_zp(pw-EEGopts.srate:pw);
% freq_slide_1sec_real=fslide_c(EEGopts.times>=-1 & EEGopts.times<=0);
% 
% % Samplig points to evaluate goodness of estimation
% npoints=75;
% 
% % pearson correlation
% disp(['Correlation coeff for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
% disp(['AR :',num2str(corr(freq_slide_1sec_AR(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
% disp(['Zero padding :',num2str(corr(freq_slide_1sec_zp(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
% 
% % mutual info
% disp(['MI for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
% disp(['AR :',num2str(gcmi_cc(freq_slide_1sec_AR(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
% disp(['Zero padding :',num2str(gcmi_cc(freq_slide_1sec_zp(end-npoints:end),freq_slide_1sec_real(end-npoints:end)))])
% 
% % RMSE sqrt(ei^2) with ei = yi - y~i
% e_AR=freq_slide_1sec_AR(end-npoints:end)-freq_slide_1sec_real(end-npoints:end);
% e_zp=freq_slide_1sec_zp(end-npoints:end)-freq_slide_1sec_real(end-npoints:end);
% disp(['RMSE for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
% disp(['AR :',num2str(sqrt(mean(e_AR.^2)))])
% disp(['Zero padding :',num2str(sqrt(mean(e_zp.^2)))])
% 
% % sMAPE = mean(200|ei|/(yi+y~i))
% div_AR=freq_slide_1sec_AR(end-npoints:end)+freq_slide_1sec_real(end-npoints:end);
% div_zp=freq_slide_1sec_zp(end-npoints:end)+freq_slide_1sec_real(end-npoints:end);
% disp(['sMAPE for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
% disp(['AR :',num2str(mean(200*abs(e_AR)./div_AR))])
% disp(['Zero padding :',num2str(mean(200*abs(e_zp)./div_zp))])



% // eof








