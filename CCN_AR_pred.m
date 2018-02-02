function [ yy,tt] = CCN_AR_pred(data2interp,time2interp,opts )
% AR interpolation as predictor

%Init data
ys=data2interp';
ts=time2interp';

% Prediciton horizon
PH=opts.secPH;
% Weights
mu=0.9;
Ts=1/opts.srate;

t=time2interp;
y=data2interp;
ns=length(y);
p=25;
P=eye(p);
ahat=(1/p)*ones(p,1);
psi=flipud(ys(1:p,1));

for k=p+1:ns
   
    y=ys(k);
    t=ts(k);
    
    % creating Kmatr
    K=(P*psi)/(mu+psi'*P*psi);
    % update P
    P=(1/mu)*(P-(P*psi*psi'*P))/(mu+psi'*P*psi);
    % compute y(N+1)
    y_pred = psi'*ahat;
    % predicion error
    e = y-y_pred;
    % estim param
    ahat = ahat + K*(e);
    % update psi
    psi=[y; psi(1:p-1,1)];
    
    % prediction
    tp=t+PH;
    psi_tmp=psi;
    for j=1:(PH/Ts)
       yp=psi_tmp'*ahat;
       psi_tmp=[yp; psi_tmp(1:length(psi_tmp)-1)];
    end
    
end

psi_tmp=psi;
for j=1:(PH/Ts)
    tt_pred(j)=t+Ts;
    yy_pred(j)=psi_tmp'*ahat;
    psi_tmp=[yy_pred(j); psi_tmp(1:length(psi_tmp)-1)];
    t=t+Ts;
end

%% Noise variance estimation

Ts=1/opts.srate;
N=length(ys);
% ar estim
model = ar(ys, p, 'yw', 'Ts', Ts);

sigma2=model.NoiseVariance;
ynoise = pinknoise(sqrt(sigma2),length(tt_pred));
yy_pred_noise=yy_pred+ynoise;

% close
% figure('Position', [0 0 1800 1000])
% subplot(211)
% plot(ts,ys)
% hold on
% plot(tt_pred,yy_pred)
% subplot(212)
% plot(ts,ys)
% hold on
% plot(tt_pred,yy_pred_noise)

%% Update

yy=[data2interp yy_pred_noise];
tt=[time2interp tt_pred];

end

