function [ avg_phase, itc ] = CCN_getITC( stft )
%[ avg_phase, itc ] = CCN_getITC( stft )
% takes 'stft' matrix [TIME x CHANNELS x TRIAL] containing complex numbers
% from a frequency analysis (e.g. fft or wavelet).
%
% Note that in a complex vector in stft: -1.6734  +  0.0837i
%                                          ampl.    phase(rad)
%                                         |       .
%                                         |      .            
%                                         |     .
%                                         |   C.
%                                        B|   . ampl. = C, pow = C^2 (abs)
%                                         |  .
%                                         | .
%                                         |. Phi = phase angle (cmplx)
%                                         .________
%                                             A
%
% ITC explanation:
% intertrial coherence (ITC, also called phase-locking factor), is a measure 
% that tells you how often a certain phase occurs at a certain time-period. 
% If at time=x, the phase always is y, then the ITC will be 1. If y is always 
% different at time x, then the ITC at that time-point is 0. In a complex plane 
% all datapoints are vectors with a certain phase (the angle) and a certain 
% amplitude. If we normalize all amplitudes (set them to 1), then all 
% datapoints will be on the unit-circle. If we then compute the (complex) 
% mean, we will find a vector with amplitude somewhere in between 0 and 1, 
% with a phase that represents the average phase at that time-point. The 
% amplitude of this vector is the ITC. 
%
% -----
% Feb,2018 David Meijer & Steffen Buergers
% CCN LAB
% University of Birmingham
% DXM472@student.bham.ac.uk 
% SXB1173@student.bham.ac.uk 


n_chan = size(stft,2);
n_trls = size(stft,3);

% go through the stft and normalize the amplitudes (C = 1)
disp('Calculating inter-trial coherence...')
stft_norm = nan(size(stft));
for c = 1:n_chan
    for i = 1:n_trls
          stft_norm(:,c,i) = stft(:,c,i)./abs(stft(:,c,i));
    end
end

% get modulus and phase angle from average of normalized amplitude vectors
itc       = squeeze(   abs(mean(stft_norm,1))  );
avg_phase = squeeze( angle(mean(stft_norm,1))  );

% % show histogram of phase angles
% figure
% rose(angle(squeeze(stft_norm(:,400))),20); hold on
% rose(avg_phase(400), 'color', 'r')
% plotspecs

return 

% // eof



















