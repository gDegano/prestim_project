clc
clear

warning off

%% CREATING DIPOLE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General settings
trl_sec=2;
EEG.srate=500;
EEG.pnts=EEG.srate*trl_sec';
EEG.times=-2+1/EEG.srate:1/EEG.srate:0;
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
    cfg.fsample =  EEG.srate;
    cfg.relnoise = .1;
    raw1 = ft_dipolesimulation(cfg);
    
    %% FREQ SLIDING
    
    freqslideFilt=zeros(length(EEG.times),n_chan,n_trl);
    for j=1:n_trl
        for k=1:n_chan
            data2use = squeeze(raw1.trial{j}(k,:));
            freq2use = 10; % hz
            % define convolution parameters
            wavt  = -.25:1/EEG.srate:.25; % time vector for wavelet
            nData = EEG.pnts;
            nKern = length(wavt);
            nConv = nData+nKern-1;
            hwave = floor((length(wavt)-1)/2);
            s = 8 / (2*pi*freq2use); % for gaussian
            cmwX = fft( exp(1i*2*pi*freq2use*wavt) .* exp( -(wavt.^2)/(2*s^2) ) ,nConv);
            cmwX = cmwX ./ max(cmwX);
            % convolution and freq slid
            as = ifft( fft(data2use,nConv) .* cmwX );
            as = as(hwave+1:end-hwave);
            freqslide = EEG.srate*diff(unwrap(angle(as)))/(2*pi);
            % now apply median filter
            n_order = 10;
            orders  = linspace(10,400,n_order)/2;
            orders  = round( orders/(1000/EEG.srate) );
            phasedmed = zeros(length(orders),EEG.pnts);
            for oi=1:n_order
                for ti=1:EEG.pnts
                    temp = sort(freqslide(max(ti-orders(oi),1):min(ti+orders(oi),EEG.pnts-1)));
                    phasedmed(oi,ti) = temp(floor(numel(temp)/2)+1);
                end
            end
            % the final step is to take the mean of medians
            freqslideFilt(:,k,j) = mean(phasedmed,1);
        end
        
    end
    AVG_median(m,:)=mean(mean(freqslideFilt,3),2);
    disp(['Computing phase: ',num2str(m)])

end

%% PLOTTING

figure 
for m=1:length(phase_b)
    LEGEND{m}=['Phase at: ',num2str(phase_b(m)),'rad'];
    plot(EEG.times(EEG.times>-1),AVG_median(m,EEG.times>-1))
    hold on 
end
legend(LEGEND{:},'Location','Best')
















