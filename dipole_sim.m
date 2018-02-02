clc
clear

warning off

%% CREATING DIPOLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General settings
trl_sec=3;
EEGopts.srate=256;
EEGopts.pnts=EEGopts.srate*trl_sec';
EEGopts.times=-2:1/EEGopts.srate:1-1/EEGopts.srate;
n_trl=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Creating setup of volumes and electrodes
elec = [];
n_chan = 64;
elec.pnt = randn(n_chan,3); % 0 to 1 on three dim
dum = sqrt(sum(elec.pnt.^2,2));
elec.pnt = elec.pnt ./ [dum dum dum];  % scale them to a unit sphere
for i=1:n_chan
   elec.label{i} = sprintf('%03d', i);
end
vol = [];
vol.r = [0.88 0.92 1.00]; % radii of spheres
vol.cond = [1 1/80 1];       % conductivity
vol.o = [0 0 0];          % center of sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating dipole with filter bank
phase_b=0:pi/4:pi;
for m=1:length(phase_b)
    cfg         = [];
    cfg.vol     = vol;             % see above
    cfg.elec    = elec;            % see above
    cfg.dip.pos = [0 0.5 0.3];
    cfg.dip.mom = [1 0 0]';     % note, it should be transposed
    cfg.dip.frequency = 10;
    cfg.dip.phase = phase_b(m);
    cfg.ntrials = n_trl;
    cfg.triallength = trl_sec;
    cfg.fsample =  EEGopts.srate;
    cfg.relnoise = .8;
    raw{m}= ft_dipolesimulation(cfg);
end

%% MODELS

clc
AVG_median=nan(length(phase_b),length(EEGopts.times)-EEGopts.srate);
EEGopts.model=1;
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

















