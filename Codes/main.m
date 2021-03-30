close all
clear all
%% Synthetic data
% load('highfreq_non_dis_non_att.mat')
% dispersive = 0;
load('highfreq_dis_non_att.mat')
dispersive = 1;
% load('highfreq_non_dis_att.mat')
% dispersive = 0;

Mused = 100;
dx = 1e-3;
dy = 1e-3;
dt = 1e-6;
df = freq_src(2)-freq_src(1); % assume multi-frequencies

% the region that should be dropped due to the source
source_x=[46,55];
source_y=[46,55];
%% Vibrating plate data
load('Xsus.mat')
get_VP_Params
%var_gauss=1e-4;
%add_noise
% the region that should be dropped due to the source
source_x=[];
source_y=[];
freq_src = freq_src(21:5:71);
dispersive = 1;
%% Classic wavnumber extraction method; PDE recovery in time and freq domain
classic_speed_circwave
WaveEq_LoopTheta
Helmholtz_LoopTheta

