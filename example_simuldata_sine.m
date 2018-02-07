% Comput frequency characteristics for actual EEG data using different
% preprocessing methods. We assume that EEG data consists of single trials
% with a stimulus onset at t=0, however, these data do not contain any
% stimulus at t=0. Within script a range of frequencies and phases for the
% simulation are defined and looped through. Output figures are saved in a
% directory specified by the user at the top of the script with sensible
% sub-folder names, which are created automatically. 
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


% Setup workspace

clc
clear

warning off

% Add dependencies
addpath('mutual_information_toolbox');

% Where do we want to save figures?
save_dir = 'H:\prestim_project_figures';



%% Creating dipole

% Do we want to plot our dipole model?
show_dipole_model = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General settings
trl_sec       = 3; % has to be 3 atm
EEGopts.srate = 256;
EEGopts.pnts  = EEGopts.srate*trl_sec';
EEGopts.times = (-trl_sec+1):1/EEGopts.srate:1;
n_trl         = 79;

% Load a canonical, realistic headmodel using FieldTrip template
headmodel = load('standard_bem.mat');

% Load standard electrode positions
elec_stnd = ft_read_sens('standard_1020.elc');
n_chan    = size(elec_stnd.label,1);
n_time    = size(EEGopts.times,2);

% Location of the dipole
ROI_dipole = [+20 -70 50]; % left posterior parietal cortex

% Relative noise level (SNR)
SNR = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating dipole with filter bank
phase_b=0:pi/4:pi;
freq_b = linspace(8,12,7);
raw = cell(length(phase_b),length(freq_b));
for m=1:length(phase_b)
    
    for n=1:length(freq_b)
        cfg = [];
        cfg.headmodel     = headmodel.vol;   % see above
        cfg.elec          = elec_stnd;       % see above
        cfg.dip.pos       = ROI_dipole;
        cfg.dip.mom       = [1 0 0]';        % note, it should be transposed
        cfg.dip.frequency = freq_b(n);
        cfg.dip.phase     = phase_b(m);
        cfg.ntrials       = n_trl;
        cfg.triallength   = trl_sec;
        cfg.fsample       = EEGopts.srate;
        cfg.relnoise      = 1/SNR;
        raw{m,n} = ft_dipolesimulation(cfg);
        
        disp(['Computing dipole at phase: ',num2str(m)])
        disp(['and frequency: ',num2str(n)])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit dipole of example simulation and visualise
if show_dipole_model
    
    % average over trials for example simulation
    avg1 = ft_timelockanalysis([], raw{m,n});
    
    % do a dipole fit of the simulated dataset
    cfg      = [];
    cfg.vol  = headmodel.vol;
    cfg.elec = elec_stnd;
    cfg.dip.pos = [0 0 0];  % initial search position
    cfg.gridsearch = 'no';
    dip1 = ft_dipolefitting(cfg, avg1);
    
    % Visualise dipole fit
    figure
    ft_plot_mesh(headmodel.vol.bnd(1), 'facecolor', 'r',    'surfaceonly', 'yes', 'facealpha', 0.1); hold on
    ft_plot_mesh(headmodel.vol.bnd(2), 'facecolor', 'g',    'surfaceonly', 'yes', 'facealpha', 0.1);
    ft_plot_mesh(headmodel.vol.bnd(3), 'facecolor', 'skin', 'surfaceonly', 'no',  'facealpha', 0.1);
    plot3(elec_stnd.chanpos(:,1), elec_stnd.chanpos(:,2), elec_stnd.chanpos(:,3), '*b', 'linewidth', 3);
    plot3(ROI_dipole(1), ROI_dipole(2), ROI_dipole(3), 'om', 'linewidth', 5)
    %plot3(dip1.dip.pos(1), dip1.dip.pos(2), dip1.dip.pos(3), 'oy', 'linewidth', 5)
    
end



%% Loop through phase and frequency simulations

for m=1:length(phase_b)
    
    for n=1:length(freq_b)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Convert simulated data to real data format
        n_chan = size(raw{m,n}.label,1);
        complete.EEG.srate = EEGopts.srate;
        complete.EEG.time  = EEGopts.times;
        prestim.EEG.srate  = EEGopts.srate;
        prestim.EEG.time   = EEGopts.times(1:2*EEGopts.srate+1);
        complete.DATA = nan(n_chan, n_time, n_trl);
        prestim.DATA = nan(n_chan, size(prestim.EEG.time,2), n_trl);
        for i = 1:n_trl
            complete.DATA(:,:,i) = [raw{m,n}.trial{i}, zeros(n_chan, 1)];
            prestim.DATA(:,:,i)  = raw{m,n}.trial{i}(:,1:2*EEGopts.srate+1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Select single channel to speed up processing
        prestim.DATA  = prestim.DATA(1,:,:);
        complete.DATA = complete.DATA(1,:,:);
        n_chan = 1;
        
        
        
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
                [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts);
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
                [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts);
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
                [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts);
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
                [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts_model);
                dataAll_arma(j,k,1:length(temp_out)) = temp_out;
                dataCmplx_arma(j,k,1:length(temp_out)) = temp_cmplx;
                freqslideFilt(1:length(temp_out),k,j) = temp_freq;
            end
            disp(['Computing trial #',num2str(j)])
        end
        fslide_arma = mean(mean(freqslideFilt,3),2);
        [ avgPhase_arma, itc_arma ] = CCN_getITC( dataCmplx_arma );
        power_arma = squeeze(mean(abs(dataCmplx_arma),1));
        
        
        
        %% Compute with AR model
        
        disp('Compute with AR model')
        
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
        [dataAll_ar, dataCmplx_ar] = deal(zeros(n_trl, n_chan, size(complete.DATA,2)));
        for j=1:n_trl
            for k=1:n_chan
                data2use = squeeze(prestim.DATA(k,:,j));
                [temp_freq, temp_out, temp_cmplx] = CCN_freq_slide(data2use,EEGopts_model);
                dataAll_ar(j,k,1:length(temp_out)) = temp_out;
                dataCmplx_ar(j,k,1:length(temp_out)) = temp_cmplx;
                freqslideFilt(1:length(temp_out),k,j) = temp_freq;
            end
            disp(['Computing trial #',num2str(j)])
        end
        fslide_ar = mean(mean(freqslideFilt,3),2);
        [ avgPhase_ar, itc_ar ] = CCN_getITC( dataCmplx_ar );
        power_ar = squeeze(mean(abs(dataCmplx_ar),1));
        
        
        
        %% Plots
        
        % Create figure save folder
        fig_dir = fullfile(save_dir, sprintf('simuldata_noise%i', round(1/SNR)), ...
            sprintf('simul_phase%i_freq%i', round(phase_b(m)*100), round(freq_b(n)*100)));
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir)
        end
        
        % size of prediction window for models
        pw   = length(temp_out);
        time = EEGopts.times;
        
        % plot input data
        fh(1) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time, squeeze(dataAll_zp(1,1,:)), 'r'); hold on      % Zero-padded
        plot(time, squeeze(dataAll_mp(1,1,:)), 'm');              % Mirror-padded
        plot(time(1:pw), squeeze(dataAll_ar(1,1,1:pw)), 'g');     % AR model
        plot(time(1:pw), squeeze(dataAll_arma(1,1,1:pw)), 'b');   % ARMA model
        plot(time, squeeze(dataAll_c(1,1,:)), 'k');               % ERP
        xlim([-0.5, 0.25]);
        title('Input data'); xlabel('Time'); ylabel('Amplitude');
        %legend({'Zero-padded', 'Mirror-padded', 'AR', 'ARMA', 'ERP'}, 'location', 'NorthWest')
        plotspecs; grid on
        saveas(fh(1), fullfile(fig_dir, 'input_data_1trial.emf'))
        close all
        
        % plot input data
        fh(1) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time, mean(squeeze(dataAll_zp(:,1,:))), 'r'); hold on      % Zero-padded
        plot(time, mean(squeeze(dataAll_mp(:,1,:))), 'm');              % Mirror-padded
        plot(time(1:pw), mean(squeeze(dataAll_ar(:,1,1:pw))), 'g');     % AR model
        plot(time(1:pw), mean(squeeze(dataAll_arma(:,1,1:pw))), 'b');   % ARMA model
        plot(time, mean(squeeze(dataAll_c(:,1,:))), 'k');               % ERP
        xlim([-0.5, 0.25]);
        title('Input data'); xlabel('Time'); ylabel('Amplitude');
        %legend({'Zero-padded', 'Mirror-padded', 'AR', 'ARMA', 'ERP'}, 'location', 'NorthWest')
        plotspecs; grid on
        saveas(fh(1), fullfile(fig_dir, 'input_data.emf'))
        
        
        % plot output data (FS)
        fh(2) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time(1:pw),fslide_zp(1:pw),  'r'); hold on     % Zero-padded
        plot(time(1:pw),fslide_mp(1:pw),  'm')              % Mirror-padded
        plot(time(1:pw),fslide_ar(1:pw),  'g')              % AR model
        plot(time(1:pw),fslide_arma(1:pw),'b')              % ARMA model
        plot(time(1:pw),fslide_c(1:pw),   'k')              % ERP
        xlim([-0.5, 0.25]); 
        title('Output data (instantaneous frequency)'); xlabel('Time'); ylabel('Frequency');
        plotspecs; grid on
        saveas(fh(2), fullfile(fig_dir, 'fslide.emf'))
        
        
        % plot output data (POWER)
        fh(3) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time(1:pw),power_zp(1:pw),  'r'); hold on     % Zero-padded
        plot(time(1:pw),power_mp(1:pw),  'm')              % Mirror-padded
        plot(time(1:pw),power_ar(1:pw),  'g')              % AR model
        plot(time(1:pw),power_arma(1:pw),'b')              % ARMA model
        plot(time(1:pw),power_c(1:pw),   'k')              % ERP
        xlim([-0.5, 0.25]);
        title('Output data (power)'); xlabel('Time'); ylabel('Amplitude');
        plotspecs; grid on
        saveas(fh(3), fullfile(fig_dir, 'power.emf'))
        
        
        % plot output data (ITC)
        fh(4) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time(1:pw),itc_zp(1:pw),  'r'); hold on     % Zero-padded
        plot(time(1:pw),itc_mp(1:pw),  'm')              % Mirror-padded
        plot(time(1:pw),itc_ar(1:pw),  'g')              % AR model
        plot(time(1:pw),itc_arma(1:pw),'b')              % ARMA model
        plot(time(1:pw),itc_c(1:pw),   'k')              % ERP
        xlim([-0.5, 0.25]);
        title('Output data (inter-trial coherence)'); xlabel('Time'); ylabel('ITC');
        plotspecs; grid on
        saveas(fh(4), fullfile(fig_dir, 'itc.emf'))
        
        close all
        
        
        % Same plots, but cut-off at t=0
        
        % plot input data
        fh(1) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time, squeeze(dataAll_zp(1,1,:)), 'r'); hold on      % Zero-padded
        plot(time, squeeze(dataAll_mp(1,1,:)), 'm');              % Mirror-padded
        plot(time(1:pw), squeeze(dataAll_ar(1,1,1:pw)), 'g');     % AR model
        plot(time(1:pw), squeeze(dataAll_arma(1,1,1:pw)), 'b');   % ARMA model
        plot(time, squeeze(dataAll_c(1,1,:)), 'k');               % ERP
        xlim([-0.5, 0]);
        title('Input data'); xlabel('Time'); ylabel('Amplitude');
        %legend({'Zero-padded', 'Mirror-padded', 'AR', 'ARMA', 'ERP'}, 'location', 'NorthWest')
        plotspecs; grid on
        saveas(fh(1), fullfile(fig_dir, 'zero_input_data_1trial.emf'))
        close all
        
        
        % plot input data
        fh(1) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time, mean(squeeze(dataAll_zp(:,1,:))), 'r'); hold on      % Zero-padded
        plot(time, mean(squeeze(dataAll_mp(:,1,:))), 'm');              % Mirror-padded
        plot(time(1:pw), mean(squeeze(dataAll_ar(:,1,1:pw))), 'g');     % AR model
        plot(time(1:pw), mean(squeeze(dataAll_arma(:,1,1:pw))), 'b');   % ARMA model
        plot(time, mean(squeeze(dataAll_c(:,1,:))), 'k');               % ERP
        xlim([-0.5, 0]);
        title('Input data'); xlabel('Time'); ylabel('Amplitude');
        %legend({'Zero-padded', 'Mirror-padded', 'AR', 'ARMA', 'ERP'}, 'location', 'NorthWest')
        plotspecs; grid on
        saveas(fh(1), fullfile(fig_dir, 'zero_input_data.emf'))
        
        
        % plot output data (FS)
        fh(2) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time(1:pw),fslide_zp(1:pw),  'r'); hold on     % Zero-padded
        plot(time(1:pw),fslide_mp(1:pw),  'm')              % Mirror-padded
        plot(time(1:pw),fslide_ar(1:pw),  'g')              % AR model
        plot(time(1:pw),fslide_arma(1:pw),'b')              % ARMA model
        plot(time(1:pw),fslide_c(1:pw),   'k')              % ERP
        xlim([-0.5, 0]); 
        title('Output data (instantaneous frequency)'); xlabel('Time'); ylabel('Frequency');
        plotspecs; grid on
        saveas(fh(2), fullfile(fig_dir, 'zero_fslide.emf'))
        
        
        % plot output data (POWER)
        fh(3) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time(1:pw),power_zp(1:pw),  'r'); hold on     % Zero-padded
        plot(time(1:pw),power_mp(1:pw),  'm')              % Mirror-padded
        plot(time(1:pw),power_ar(1:pw),  'g')              % AR model
        plot(time(1:pw),power_arma(1:pw),'b')              % ARMA model
        plot(time(1:pw),power_c(1:pw),   'k')              % ERP
        xlim([-0.5, 0]);
        title('Output data (power)'); xlabel('Time'); ylabel('Amplitude');
        plotspecs; grid on
        saveas(fh(3), fullfile(fig_dir, 'zero_power.emf'))
        
        
        % plot output data (ITC)
        fh(4) = figure('color', 'w', 'position', [50 50 500 250]);
        plot(time(1:pw),itc_zp(1:pw),  'r'); hold on     % Zero-padded
        plot(time(1:pw),itc_mp(1:pw),  'm')              % Mirror-padded
        plot(time(1:pw),itc_ar(1:pw),  'g')              % AR model
        plot(time(1:pw),itc_arma(1:pw),'b')              % ARMA model
        plot(time(1:pw),itc_c(1:pw),   'k')              % ERP
        xlim([-0.5, 0]);
        title('Output data (inter-trial coherence)'); xlabel('Time'); ylabel('ITC');
        plotspecs; grid on
        saveas(fh(4), fullfile(fig_dir, 'zero_itc.emf'))
        
        close all
        
        
        
        %% Goodness of fit
        
        time_1sec=EEGopts.times(EEGopts.times>=-1 & EEGopts.times<=0);
        freq_slide_1sec_AR=fslide_ar(pw-EEGopts.srate:pw);
        freq_slide_1sec_zp=fslide_zp(pw-EEGopts.srate:pw);
        freq_slide_1sec_real=fslide_c(EEGopts.times>=-1 & EEGopts.times<=0);
        
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
        disp(['RMSE for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
        disp(['AR :',num2str(sqrt(mean(e_AR.^2)))])
        disp(['Zero padding :',num2str(sqrt(mean(e_zp.^2)))])
        
        % sMAPE = mean(200|ei|/(yi+y~i))
        div_AR=freq_slide_1sec_AR(end-npoints:end)+freq_slide_1sec_real(end-npoints:end);
        div_zp=freq_slide_1sec_zp(end-npoints:end)+freq_slide_1sec_real(end-npoints:end);
        disp(['sMAPE for last ',num2str(round(npoints/EEGopts.srate*1000)),' ms: '])
        disp(['AR :',num2str(mean(200*abs(e_AR)./div_AR))])
        disp(['Zero padding :',num2str(mean(200*abs(e_zp)./div_zp))])
        
        
    end % frequency loop

end % phase loop

% // eof








