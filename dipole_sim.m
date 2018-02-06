clc
clear

warning off


%% CREATING DIPOLE

% Do we want to plot our dipole model?
show_dipole_model = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General settings
trl_sec       = 3; % has to be 3 atm  
EEGopts.srate = 256;
EEGopts.pnts  = EEGopts.srate*trl_sec';
EEGopts.times = -trl_sec+1/EEGopts.srate:1/EEGopts.srate:0;
n_trl         = 10;

% Load a canonical, realistic headmodel using FieldTrip template
headmodel = load('standard_bem.mat');

% Load standard electrode positions
elec_stnd = ft_read_sens('standard_1020.elc');
n_chan     = size(elec_stnd.label,1);

% Location of the dipole
ROI_dipole = [+20 -70 50]; % left posterior parietal cortex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating dipole with filter bank
phase_b=0:pi/4:pi;
raw = cell(length(phase_b),1);
for m=1:length(phase_b)
    
    cfg = [];
    cfg.headmodel     = headmodel.vol;   % see above
    cfg.elec          = elec_stnd;       % see above
    cfg.dip.pos       = ROI_dipole;
    cfg.dip.mom       = [1 0 0]';        % note, it should be transposed
    cfg.dip.frequency = 10;
    cfg.dip.phase     = phase_b(m);
    cfg.ntrials       = n_trl;
    cfg.triallength   = trl_sec;
    cfg.fsample       = EEGopts.srate;
    cfg.relnoise      = 0;
    raw{m} = ft_dipolesimulation(cfg);
    
    disp(['Computing dipole at phase: ',num2str(m)])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Fit dipole of example simulation and visualise
if show_dipole_model
    
    % average over trials for example simulation
    avg1 = ft_timelockanalysis([], raw{m});
    
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


%% MODELS

clc
AVG_median=nan(length(phase_b),length(EEGopts.times)-EEGopts.srate);
EEGopts.model = []; % 1 = AR; 2 = ARMA
EEGopts.method     = 'hilbert';
EEGopts.plotfilt   = true;
EEGopts.filtwin    = [6 14];
EEGopts.transwidth = 0.15;
EEGopts.figures    = true;
for m=1:length(phase_b)
   
    disp(['Computing phase: ',num2str(m)])
    % FREQ SLIDING
    freqslideFilt=zeros(length(EEGopts.times)-EEGopts.srate,n_chan,n_trl);
    for j=1:n_trl
        for k=1:n_chan
            data2use = squeeze(raw{m}.trial{j}(k,1:EEGopts.srate*2));
            central_freq=10;
            num_cycles=5;
            temp_freq=CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
            freqslideFilt(:,k,j) = temp_freq(1:EEGopts.srate*2);
            pause(0.5)
        end
        disp(['Computing trial #',num2str(j)])
    end
    AVG_median(m,:)=mean(mean(freqslideFilt,3),2);

end

% PLOTTING
figure 
for m=[1 2 3 5]
    LEGEND{m}=['Phase at: ',num2str(phase_b(m)),'rad'];
    plot(EEGopts.times(EEGopts.srate*+1:EEGopts.srate*2),AVG_median(m,end-EEGopts.srate:end))
    hold on 
end
legend(LEGEND{:},'Location','Best')


%% REAL DATA (GT)

AVG_median_real=nan(length(phase_b),length(EEGopts.times));
EEGopts.model=0;
for m=1:length(phase_b)
   
    disp(['Computing phase: ',num2str(m)])
    % FREQ SLIDING
    freqslideFilt=zeros(length(EEGopts.times),n_chan,n_trl);
    for j=1:n_trl
        for k=1:n_chan
            data2use = squeeze(raw{m}.trial{j}(k,:));
            central_freq=10;
            num_cycles=5;
            freqslideFilt(:,k,j) = CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        end
        disp(['Computing trial #',num2str(j)])
    end
    AVG_median_real(m,:)=mean(mean(freqslideFilt,3),2);

end

%plots
figure 
for m=[1 2 3 5]
    LEGEND{m}=['Phase at: ',num2str(phase_b(m)),'rad'];
    plot(EEGopts.times(EEGopts.times>-1 & EEGopts.times<=0),AVG_median_real(m,EEGopts.times>-1 & EEGopts.times<=0))
    hold on 
end
legend(LEGEND{:},'Location','Best')

% // eof













