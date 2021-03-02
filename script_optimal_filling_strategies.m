% Script for optimization of informed strategies.

% NB. this is a demo version of the code used to produce the results
% displayed in the paper. In the following code, the employed optimization 
% algorithm is the Matlab built-in gamultiobj, while we employed BorgMOEA, 
% to produce the paper results. BorgMOEA is available at borgmoea.org

clc
clear

%% load variables for objectives computation
global sim reservoir_specifics SPEI policy sys_param spei_date
load Hsdp.mat
load optimization_data.mat
load SPEIforecast.mat
addpath('SDP')
addpath('data')
addpath('utils')

n_drought_categories = 3; 
% drought categories for SPEI6
% 1: spei > 0.5          wet 
% 2: -0.5 < spei < 0.5   normal
% 3: spei < -0.5  dry

%% optimization via gamultiobj

Nvar = n_drought_categories; 
lb = zeros(1,n_drought_categories);
ub = ones(1,n_drought_categories);
options = optimoptions('gamultiobj', 'Generations', 20); % demo optimization

[k, J] = gamultiobj(@obj_func, Nvar,[],[],[],[],lb,ub, options);
figure; scatter(J(:,2), -J(:,1), -J(:,3)/2, 'filled')
xlabel('environmental devaition [(m^3/s)^2]')
ylabel('Hydropower Production [MWh/year]')
title('Optimization results')
ylim([0,2300])


