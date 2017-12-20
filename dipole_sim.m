clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create a set of electrodes, randomly placed on the sphere
elec = [];
n_chan = 64;
elec.pnt = randn(n_chan,3); % 0 to 1 on three dim
dum = sqrt(sum(elec.pnt.^2,2));
elec.pnt = elec.pnt ./ [dum dum dum];  % scale them to a unit sphere
for i=1:n_chan
   elec.label{i} = sprintf('%03d', i);
end

% create a concentric 3-sphere volume conductor, the radius is the same as for the electrodes
vol = [];
vol.r = [0.88 0.92 1.00]; % radii of spheres
vol.cond = [1 1/80 1];       % conductivity
vol.o = [0 0 0];          % center of sphere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% create a dipole simulation with one dipole and a 10Hz sine wave
cfg      = [];
cfg.vol  = vol;             % see above
cfg.elec = elec;            % see above
cfg.dip.pos = [0 0.5 0.3];
cfg.dip.mom = [1 0 0]';     % note, it should be transposed
cfg.dip.frequency = 10;
cfg.ntrials = 100;
cfg.triallength = 1;        % seconds
cfg.fsample = 500;          % Hz
cfg.relnoise = .5;

raw1 = ft_dipolesimulation(cfg);

cfg = [];
cfg.channel = 'all';
cfg.viewmode   = 'vertical';
cfg.continuous = 'no';
ft_databrowser(cfg,raw1)


























