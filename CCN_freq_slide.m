function [ instan_freq, yy, cmplx_data ] = CCN_freq_slide(data2use,opts,central_freq,num_cycles)
% Frequency sliding extraction following the method outlined in Cohen (2014):
%     Fluctuations in oscillation frequency control spike timing and coordinate 
%     neural networks. The Journal of Neuroscience, 2 July 2014, 34(27): 8988-8998; 
%     doi: 10.1523/JNEUROSCI.0261-14.2014
%
% Instantaneous frequency is defined as the change in phase per unit time, 
% i.e. the temporal derivative of the phase angle time-series. 
%
% INPUT : data2use     = 1xN time vector 
%         opts.pnts    = number of datapoints
%         opts.srate   = sampling freq
%         opts.model   = 1=AR (autoregressive model); 2=ARMA
%                        (autoregressive model with moving average
%                        component); Default is [], no modeling!
%         opts.method  = how to compute frequency information, can be
%                        'wavelet' or 'hilbert' (default = 'hilbert')
% If opts.method is 'hilbert', you are able to additionally specifiy 
% opts.freqwin (default = [6 14] Hz) for the filter. Plus, opts.transwidth,
% which specifies the transition width of the Tukey filter (default =
% 0.15). You can further specify if the filter response should be plotted
% with opts.plotfilt (default = false).
%
%         opts.figures = plot frequency sliding output (default = false)
%         central_freq = central frequency of the wavelet
%         num_cycles   = number of cycles in the window (default 4)
%
% OUTPUT: instan_freq  = 1xN time vector with instantaneous frequency
%         yy           = modeled data (corresponds to input data if
%                        opts.model is not in [1,2])
%         cmplx_data   = Complex signal after hilbert transform or wavelet
%                        convolution (used to compute instan_freq)
%
% ------
% FS code adapted from public online resources by Mike X Cohen
% www.mikexcohen.com/
% 
% Jan,2018 Giulio Degano & Steffen Buergers
% CCN LAB
% University of Birmingham
% GXD606@student.bham.ac.uk 
% SXB1173@student.bham.ac.uk 


%% Defaults

% method (for going to frequency domain)
if ~isfield(opts, 'method')
    opts.method = 'hilbert';
end

% set defaults for hilbert method
if strcmp(opts.method, 'hilbert')
    if ~isfield(opts, 'freqwin')  % Frequency bounds of filter
        opts.freqwin = [6 14];
    end
    if ~isfield(opts, 'transwidth') % Filter transition widths
        opts.transwidth = 0.15;
    end
    if ~isfield(opts, 'plotfilt') % Plot filter response?
        opts.plotfilt = false;
    end
end

% wavelet frequency
if nargin<3
    central_freq = 10; 
end

% wavelet cycles
if nargin<4
    num_cycles=4;
end

% plot figures or not
if ~isfield(opts, 'figures')
    opts.figures = false;
end

% length of window to predict (in seconds)
if ~isfield(opts, 'secPH')
    opts.secPH=.5; 
end

% Prediction method
if opts.model==1    
    [yy] = CCN_AR_pred(data2use,opts.times,opts);
elseif opts.model==2
    [yy] = CCN_ARMAX_pred(data2use,opts.times,opts);
else
    yy=data2use;
end


%% FS computation

% There are two approaches for calculating instantaneous frequency
% (aka frequency sliding) based on a wavelet approach or the Hilbert
% transform. The difference is that with the Hilbert transform the time
% series signal is filtered using a Tukey window with 15% transition zones,
% whereas the wavelet has a much smoother filter response that is simply
% broadened or widened by the number of specified cycles. For illustrative
% purpose the method does not really matter. But if you interested in e.g.
% alpha band activity and want to retain information throughout most of the
% frequency band of interest, the hilbert method should be preferred. Both
% pipelines ultimately yield a complex, band pass filtered signal centered
% around the frequency of interest.
%
% The derivative of the corresponding phase angle time series is unwrapped
% and normalised with 2*pi to obtain a time series of instantaneous 
% frequency (showing the dominant frequency over time). 

% create some useful variables to make code more readable
nData     = length(yy);
fs        = opts.srate;

if strcmp(opts.method, 'wavelet')
    
    % define convolution parameters
    freq2use  = central_freq; % hz
    wavt  = -((1/freq2use)*num_cycles)/2:1/fs:((1/freq2use)*num_cycles)/2; % time vector for wavelet
    nKern = length(wavt);
    nConv = nData+nKern-1;
    hwave = floor((length(wavt)-1)/2);
    s     = 8 / (2*pi*freq2use); % for gaussian
    cmwX  = fft( exp(1i*2*pi*freq2use*wavt) .* exp( -(wavt.^2)/(2*s^2) ) ,nConv);
    cmwX  = cmwX ./ max(cmwX);

    % convolution and freq slides
    cmplx_data = ifft( fft(yy,nConv) .* cmwX );
    cmplx_data = cmplx_data(hwave+1:end-hwave);
    freqslide = fs*diff(unwrap(angle(cmplx_data)))/(2*pi);

elseif strcmp(opts.method, 'hilbert')
    
    % apply a band-pass filter with 15% transition zones.
    fwin        = opts.freqwin;
    trns_wdth   = opts.transwidth;
    idealrsp    = [ 0 0 1 1 0 0 ];
    filtfrqbnds = [ 0 (1-trns_wdth)*fwin(1) fwin(1) fwin(2) fwin(2)*(1+trns_wdth) fs/2 ]/(fs/2);
    filt_ord    = round(3*(fs/fwin(1)));
    filtwghts   = firls(filt_ord,filtfrqbnds,idealrsp);
    
    if opts.plotfilt
        % plot filter response
        nyquist = fs/2;
        figure('color', 'w')
        plot(filtwghts,'r')
        xlabel('Time')
        title('Filter response');
        
        figure('color', 'w')
        plot(filtfrqbnds*nyquist,idealrsp,'r'); hold on
        fft_filtkern  = abs(fft(filtwghts));
        hz_filtkern   = linspace(0,nyquist,ceil(length(fft_filtkern)/2));
        fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
        plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')
        set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
        xlabel('Frequency');
        legend({'ideal';'applied'}); xlim([1 35])
    end
    
    % this part does the actual filtering
    filtdata = filtfilt(filtwghts, 1, double(yy));
    
    % hilbert transform
    cmplx_data = hilbert(filtdata);
    
    % FS
    freqslide = fs*diff(unwrap(angle(cmplx_data)))/(2*pi);
    
end

% now apply median filter
n_order = 10;
orders  = linspace(10,400,n_order)/2;
orders  = round( orders/(1000/fs) );
phasedmed = zeros(length(orders),nData);
for oi=1:n_order
    for ti=1:nData
        temp = sort(freqslide(max(ti-orders(oi),1):min(ti+orders(oi),nData-1)));
        phasedmed(oi,ti) = temp(floor(numel(temp)/2)+1);
    end
end
% the final step is to take the mean of medians
instan_freq = mean(phasedmed,1);


%% Optional figures

if opts.figures

    figure, clf
    subplot(311)
    plot(opts.times(1:length(yy)),yy);
    ylabel('Amplitude (\muV)')

    subplot(312)
    plot(opts.times(1:length(yy)),angle(cmplx_data))
    ylabel('Phase angles (rad.)')

    subplot(313)
    freqslide = fs*diff(unwrap(angle(cmplx_data)))/(2*pi);
    plot(opts.times(1:length(yy)-1),freqslide)
    ylabel('Frequency (Hz)')
    set(gca,'ylim',[freq2use*.75 freq2use*1.5])

    subplot(313)
    hold on
    plot(opts.times(1:length(yy)),instan_freq,'r')
    
end


return

% //eof
















