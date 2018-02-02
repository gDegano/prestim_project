function [ yy,tt] = CCN_param_estim(data2interp,time2interp,opts )
% Modelling of time series data

%Init data
ys=data2interp;
Ts=1/opts.srate;
N=length(ys);
p=25;
% ar estim
model = ar(ys, p, 'ls', 'Ts', Ts);
C=[model.A'; zeros(N-(p+1),1)];
R=[1 zeros(1,N-1)];
A=toeplitz(C,R);

sigma2=model.NoiseVariance;
B=inv(A'*A);
gamm=[];
m=4;

[uhat, resp, gamm, qgamm, itergamm]= CNN_discrep_methods(ys, B, sigma2, m,gamm);

% close
% figure('Position', [0 0 1800 1000]);
% plot(time2interp,data2interp)
% hold on
% plot(time2interp,uhat)
% title([num2str(gamm) ,' ',num2str(itergamm)])

end
