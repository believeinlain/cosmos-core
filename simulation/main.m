% Main SPH code:

clear; clc




%%%%%%%%%%%%%
% SPH setup %
%%%%%%%%%%%%%


%SPH simulation parameters
param = struct(...
'ndim',{2},... %dimension of the simulation (2 or 3)
'gain',struct('sph',{1},'ext',{.25},'drag',{10}),... %gain coefs
'accel',struct('veh',{1},'obs',{1},'rd',{.25}),... %scales accel due to vehicles/obstacles/reduced density particles
'Re',{10},... %Reynolds number
'dt', 0.01 ... %timestep for SPH simulation
);


%The groups of vehicle, obstacles, and reduced density particles
%(dimensional parameters)
group_conf = struct(...
'num_veh',{15},... % Array containing the number of vehicles in each group
'veh_init',struct('x',{0},... %initial positions for the veh. groups
                  'y',{0},...
                  'z',{0},...
                  'u',{0},... %initial velocities for the veh. groups
                  'v',{0},...
                  'w',{0}),...
'veh_h',{2},... % Smoothing width for each group
'veh_limits',struct('vmin',{3},... %limits for speed and acceleration
                    'vmax',{6},...
                    'turning_radius',{.1}),...
'num_obs',{5},...    % total number of obstacle particles
'obs_h',{[2 2 2 2 2 2]},...    % .5*size of obstacle particles
'obs_init',struct('x',{[7 12 16 9 22]},... %positions for the obstacles
                  'y',{[0 4 2 -2 -6]},... - SYNC WITH OBX
                  'z',{[0 0 0 0 0 0]}),...- SYNC WITH OBX
'num_rd',{0},...     % total number of reduced density particles
'rd_group',{1},...% which group does each red. density particle belong to?
...                  % group number corresponds to array index for num_veh,
...                  % 0 means not active
'rd_init',struct('x',{0},... %initial positions for the rd 
                  'y',{0},...
                  'z',{0},...
                  'u',{0},... %initial velocities for the rd
                  'v',{0},...
                  'w',{0}),...
'rd_h',{30},...
'num_loiter',{1},...     % total number of loiter circles
'loiter_group',{1}...% which group does each loiter circle belong to?
...                  % group number corresponds to array index for num_veh,
...                  % 0 means not active
);

%initialize the SPH simulation
SPH = sph_sim(param,group_conf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
