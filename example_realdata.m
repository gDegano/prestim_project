clc
clear

warning off

prestim=load('eeg_sample_to_model.mat'); % data to model (only prestim 2 sec)
complete=load('eeg_sample.mat'); % complete data (3 sec)

% Each files contains the following data:
%   DATA = 3D matrix (chan x time x trl)
%   EEG  = structure containing Fs rate and time vector

%% Compute with model

n_trl=size(prestim.DATA,3);
n_chan=size(prestim.DATA,1);
EEGopts.pnts=size(prestim.DATA,2);
EEGopts.srate=prestim.EEG.srate;
EEGopts.times=prestim.EEG.time;
EEGopts.model=1;

freqslideFilt=zeros(size(prestim.DATA,2),n_chan,n_trl);
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(prestim.DATA(k,:,j));
        central_freq=3.125;
        num_cycles=5;
        temp_freq=CCN_freq_slide(data2use,EEGopts,central_freq,num_cycles);
        freqslideFilt(:,k,j) = temp_freq(1:EEGopts.pnts);
        pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end

AVG_median=mean(mean(freqslideFilt,3),2);

figure
plot(EEGopts.times(end-EEGopts.srate:end),AVG_median(end-EEGopts.srate:end))


%% Compute with model

EEGopts_c.pnts=size(complete.DATA,2);
EEGopts_c.srate=complete.EEG.srate;
EEGopts_c.times=complete.EEG.time;
EEGopts_c.model=0;

freqslideFilt_c=zeros(size(complete.DATA,2),n_chan,n_trl);
for j=1:n_trl
    for k=1:n_chan
        data2use = squeeze(complete.DATA(k,:,j));
        central_freq=3.125;
        num_cycles=5;
        temp_freq=CCN_freq_slide(data2use,EEGopts_c,central_freq,num_cycles);
        freqslideFilt_c(:,k,j) = temp_freq(1:EEGopts_c.pnts);
        %pause(0.5)
    end
    disp(['Computing trial #',num2str(j)])
end

AVG_median_c=mean(mean(freqslideFilt_c,3),2);

figure
plot(EEGopts_c.times(EEGopts_c.times>=-1 & EEGopts_c.times<=0),AVG_median_c(EEGopts_c.times>=-1 & EEGopts_c.times<=0))












