function [ yy,tt,model] = CCN_AR_pred(data2interp,time2interp,opts )
% Autoregressive model prediction (AR)
% 
% Takes a time series vector as input (1xN) and computes an
% autoregressive model specified in 'opts' forecasting 'time2interp'
% datapoints
%
% INPUT : data2interp  = 1xN time vector of data (e.g. EEG signal)
%         time2interp  = 1xN vector with time values (in seconds)
%         opts.srate   = sampling freq
%         opts.secPH   = length of the data segment to predict in seconds
%         opts.order   = number of coefficients of the AR model (default=7)
%
% OUTPUT: yy           = time course of original data plus forecast
%         tt           = timeline of yy
%         model        = AR model object
%
% ------
% Jan,2018 Giulio Degano & Steffen Buergers
% CCN LAB
% University of Birmingham
% GXD606@student.bham.ac.uk 
% SXB1173@student.bham.ac.uk 


%% defaults

if ~isfield(opts, 'order')
    p = 7;
else
    p = opts.order;
end


%% AR model

% Period
Ts=1/opts.srate;
% PH
numPeriods=round(opts.secPH*opts.srate);
ys=data2interp;

% To data obj
data = iddata(ys',[],Ts);
% Modelling....
model = ar(data, p, 'yw', 'Ts', Ts);
sigma2=model.NoiseVariance;
% Forecasting
yy = forecast(model,data,numPeriods);

tt_pred=time2interp(end)+Ts:Ts:opts.secPH;
yy_pred=yy.y';


%% Noise computation 

flag_noise=1;

if flag_noise==1
    ynoise = pinknoise(sqrt(sigma2),length(tt_pred));
    yy_pred_noise=yy_pred+ynoise;
else
    yy_pred_noise=yy_pred;
end


%% Concatenate 

yy=[data2interp yy_pred_noise];
tt=[time2interp tt_pred];

% close
% ts=time2interp;
% figure('Position', [0 0 1800 1000])
% subplot(211)
% plot(ts,ys,'r')
% hold on
% plot(tt_pred,yy_pred)
% subplot(212)
% plot(ts,ys,'r')
% hold on
% plot(tt_pred,yy_pred_noise)


end