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
ts=time2interp;
ys=data2interp;
p=25;

% To data obj
data = iddata(ys',[],Ts);
% Modelling....
model = ar(data, p, 'yw', 'Ts', Ts);
sigma2=model.NoiseVariance;
% Forecasting
yy = forecast(model,data,numPeriods);

tt_pred=time2interp(end)+Ts:Ts:opts.secPH;
yy_pred=yy.y';
ynoise = pinknoise(sqrt(sigma2),length(tt_pred));
yy_pred_noise=yy_pred+ynoise;

yy=[data2interp yy_pred_noise];
tt=[time2interp tt_pred];

close
figure('Position', [0 0 1800 1000])
subplot(211)
plot(ts,ys,'r')
hold on
plot(tt_pred,yy_pred)
subplot(212)
plot(ts,ys,'r')
hold on
plot(tt_pred,yy_pred_noise)


end