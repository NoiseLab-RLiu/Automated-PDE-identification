%% Automated Partial Differential Equation Identification
%% Ruixian Liu, 07/26/2021
%% Synthetic waves
close all
clear all
load('highfreq_non_dis_non_att.mat') % non-dispersive non-attenuated waves
dispersive = 0;
% load('highfreq_dis_non_att.mat') % dispersive non-atteunated waves
% dispersive = 1;
% load('highfreq_non_dis_att.mat')  % non-dispersive attenuated waves
% dispersive = 0;

Mused = 100;
dx = 1e-3;
dy = 1e-3;
dt = 1e-6;
df = freq_src(2)-freq_src(1); % assume multi-frequencies

% the region that should be dropped due to the source
source_x=[46,55];
source_y=[46,55];

WaveEq_LoopTheta
Helmholtz_LoopTheta % work for non-attenuated waves
%% Vibrating plate data
close all
clear all
%load('Xsus.mat')
load('VP_Uf.mat')
%get_VP_Params

% the region that should be dropped due to the source
source_x=[];
source_y=[];
freq_src = freq_src(21:5:71);
%freq_src = freq_src(31:15:61);
dispersive = 1;

classic_speed_circwave % classic wavenumber extraction
WaveEq_LoopTheta

%% Burgers Equations identification
close all
clear all
dataset_index = 1; % choose the burgers equation dataset, ranged from 1 to 3
Burgers_identi_1D