function [ yy,tt] = CCN_AR_pred(data2interp,time2interp,opts )
% ARIMA prediction
%
%
% ------
% Jan,2018 Giulio Degano & Steffen Buergers
% CCN LAB
% University of Birmingham
% GXD606@student.bham.ac.uk 
% SXB1173@student.bham.ac.uk 


% Period
Ts=1/opts.srate;
% PH
numPeriods=round(opts.secPH*opts.srate);
ys=data2interp;
p=15;

% To data obj
data = iddata(ys',[],Ts);
% Modelling....
model = ar(data, p, 'yw', 'Ts', Ts);
sigma2=model.NoiseVariance;
% Forecasting
yy = forecast(model,data,numPeriods);

tt_pred=time2interp(end)+Ts:Ts:opts.secPH;
yy_pred=yy.y';

%% NOISE COMPUTATION...

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