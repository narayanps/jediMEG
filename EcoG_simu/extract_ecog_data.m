addpath('/scratch/nbe/braintrack/ecog/SubjectNY394')
addpath('/m/nbe/scratch/braintrack/fieldtrip')
ft_defaults
%% load electrode locations
fid = fopen('NY394_MRI_coor.txt');
elec_info = textscan(fid,'%s %f %f %f %s');
fclose(fid);

% create FieldTrip electrode structure
% elec       = [];
elec.label   = elec_info{1};
elec.elecpos = [elec_info{2} elec_info{3} elec_info{4}];
elec.unit    = 'mm';

%% load pial surface
load('NY394_MRI_rh_pial_surface.mat');

% create FieldTrip surface mesh structure
mesh      = [];
mesh.pos  = surface.pos;
mesh.tri  = surface.tri;
mesh.unit = 'mm';

%% plot surface and electrodes
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none')
view([90 25])
lighting gouraud
material shiny
camlight

% plot electrodes
hs = ft_plot_sens(elec, 'style', 'ko', 'label', 'on');
set(hs, 'MarkerFaceColor', 'k', 'MarkerSize', 6);


% load trial info
load('NY394_trl.mat');

% load and segment data
cfg            = [];
cfg.dataset    = 'NY394_VisualLoc_R1.edf';
cfg.trl        = trl; % from NY394_trl.mat
cfg.continuous = 'yes';
epoch_data = ft_preprocessing(cfg);


cfg         = [];
cfg.channel = 'EEG*'; % select 'EEG' channles
epoch_data = ft_selectdata(cfg,epoch_data);

cfg         = [];
cfg.method  = 'channel'; % browse through channels
cfg.channel = 'all';
epoch_data_clean_chan = ft_rejectvisual(cfg, epoch_data);

cfg         = [];
cfg.method  = 'summary'; % summary statistics across channels and trials
cfg.channel = 'all';
epoch_data_clean = ft_rejectvisual(cfg, epoch_data_clean_chan);

cfg = [];
cfg.resamplefs  = 250;
cfg.demean      = 'no';
cfg.detrend     = 'no';
epoch_data_clean = ft_resampledata(cfg, epoch_data_clean);


% calculate ERPs
cfg                  = [];
cfg.keeptrials       = 'yes'; % keep trials for statistics
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq   = 40;    % smooth ERP with low-pass filter
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 0.5;     % reduce slow drifts
cfg.preproc.detrend  = 'yes';


cfg.trials = find(epoch_data_clean.trialinfo == 3); % select only 'object' trials (event code 3)
ERP_object = ft_timelockanalysis(cfg, epoch_data_clean);

cfg.trials = find(epoch_data_clean.trialinfo == 7); % select only 'face' trials (event code 7)
ERP_face   = ft_timelockanalysis(cfg, epoch_data_clean);

% baseline correction
cfg          = [];
cfg.baseline = [-.3 -.05];

ERP_object_bl = ft_timelockbaseline(cfg,ERP_object);
ERP_face_bl   = ft_timelockbaseline(cfg,ERP_face);

cfg            = [];
cfg.avgoverrpt = 'yes';
ERP_object_avg = ft_selectdata(cfg, ERP_object_bl);
ERP_face_avg   = ft_selectdata(cfg, ERP_face_bl);


cfg           = [];
cfg.parameter = 'trial';
cfg.xlim      = [-.3 .6];
cfg.channel   = 'EEG PT_04-REF'; % other responsive channels: 'EEG PT_04-REF', 'EEG IO_02-REF', 'EEG IO_04-REF', 'EEG SO_01-REF', 'EEG SO_02-REF''EEG SO_03-REF'

figure, ft_singleplotER(cfg,ERP_face_avg)

