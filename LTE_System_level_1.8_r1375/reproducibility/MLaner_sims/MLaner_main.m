warning off;
% close all force;
clc;
clear all
clear global;
clear classes;

%% Load parameters. Now done outside
sim_num = 1;
MLaner_params;
% LTE_load_params_LLvsSL;
print_log(1,'Loaded configuration file\n');
MLaner_sim;
% LTE_sim_main_LLvsSL;

