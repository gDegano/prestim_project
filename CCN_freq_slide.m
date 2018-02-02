function [ istan_freq ] = CCN_freq_slide(data2use,opts,central_freq,num_cycles)
% Frequency sliding extraction
% INPUT : data2use = 1xN time vector (channel/sensor/whatever)
%         opts.pnts = number of datapoints
%         opts.srate = sampling freq
%         central_freq = central frequency of the wavelet
%         num_cycles = number of cycles in the window (default 4)
%
%
% ------
% Jan,2018 Giulio Degano & Steffen Buergers
% CCN LAB
% University of Birmingham
% GXD606@student.bham.ac.uk 
% SXB1173@student.bham.ac.uk 

if nargin<4
    num_cycles=4;
end

%% Prediction method

if opts.model==1
    opts.secPH=.2;
    [yy] = CCN_AR_pred(data2use,opts.times(1:opts.srate*2),opts);
else
    yy=data2use;
end


%% FS computation
freq2use = central_freq; % hz
nData = length(yy);

% define convolution parameters
wavt  = -((1/freq2use)*num_cycles)/2:1/opts.srate:((1/freq2use)*num_cycles)/2; % time vector for wavelet
nKern = length(wavt);
nConv = nData+nKern-1;
hwave = floor((length(wavt)-1)/2);
s = 8 / (2*pi*freq2use); % for gaussian
cmwX = fft( exp(1i*2*pi*freq2use*wavt) .* exp( -(wavt.^2)/(2*s^2) ) ,nConv);
cmwX = cmwX ./ max(cmwX);
% convolution and freq slid
as = ifft( fft(yy,nConv) .* cmwX );
as = as(hwave+1:end-hwave);
freqslide = opts.srate*diff(unwrap(angle(as)))/(2*pi);

% close
% plot(freqslide)

% now apply median filter
n_order = 10;
orders  = linspace(10,400,n_order)/2;
orders  = round( orders/(1000/opts.srate) );
phasedmed = zeros(length(orders),nData);
for oi=1:n_order
    for ti=1:nData
        temp = sort(freqslide(max(ti-orders(oi),1):min(ti+orders(oi),nData-1)));
        phasedmed(oi,ti) = temp(floor(numel(temp)/2)+1);
    end
end
% the final step is to take the mean of medians
istan_freq = mean(phasedmed,1);


end

