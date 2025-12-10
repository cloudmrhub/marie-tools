% ===============================================================
% FILE: CO_SIMULATION_runner.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Performs a standalone co-simulation optimization using precomputed 
%   Y-parameters and stored EM fields. Intended for users who already 
%   have EM field data and network parameters and wish to perform 
%   circuit-level optimization only.
%
% SYNTAX:
%   Run directly in MATLAB:
%       >> CO_SIMULATION_runner
%
% DESCRIPTION:
%   This script executes the co-simulation and circuit optimization stage:
%       1. Load precomputed network parameters (Y-matrix)
%       2. Read coil and lumped-element configurations
%       3. Identify coil excitation type (Tx, Rx, TxRx, hybrid)
%       4. Perform PSO-based tuning, matching, and decoupling optimization
%       5. Update lumped-element values
%       6. Store optimized network parameters and co-simulation results
%
% VERSION:
%   MARIE 3.0 (2025)
%
% NOTES:
%   - Requires precomputed EM fields and Y-parameters.
%   - The user must specify:
%       B0             - Main magnetic field strength [T]
%       YPn_mat        - .mat file containing the network parameters (Y-matrix)
%       wire_file_inp  - JSON coil file with lumped-element definitions
%   - Output includes tuned and matched network parameter matrices:
%       YP_s, ZP_s, SP_s, SP_p, 
%     A transformation matrix for the EM fields
%       M_cal, 
%     and a struct for plotting the network parameters for a frequency
%     sweep
%       figure_freq
%   - The optimization is performed using Particle Swarm Optimization (PSO)
%     as defined in co_simulation_optimization_settings().
%
% ---------------------------------------------------------------
% STEP 0: Initialize Environment
% ---------------------------------------------------------------
close all
clearvars
clc

% ---------------------------------------------------------------
% STEP 1: User Inputs
% ---------------------------------------------------------------
B0            = 7.0037;                          % Main magnetic field (T)
YPn_mat       = 'YPn_mat.mat';              % Y-parameter file
coil_file_inp = 'Head_Receive_Array_31_Coils_CST.json';  % Coil JSON file

% ---------------------------------------------------------------
% STEP 2: Load Constants and Network Parameters
% ---------------------------------------------------------------
emc = em_constants(B0,'1H');  % Physical constants based on B0
load(fullfile('./data/coils/external_coils', YPn_mat), 'YPn'); 

% ---------------------------------------------------------------
% STEP 3: Read Lumped Elements
% ---------------------------------------------------------------
coil_file = fullfile('./data/coils/external_coils', coil_file_inp);
RLC_org   = geo_scoil_lumped_elements(coil_file, 0);  % Initial elements
RLC_tmd   = geo_scoil_lumped_elements(coil_file, 1);  % Variable elements

% ---------------------------------------------------------------
% STEP 4: Load PSO Optimization Settings
% ---------------------------------------------------------------
[options_PSO_1, options_PSO_2] = co_simulation_optimization_settings();

% ---------------------------------------------------------------
% STEP 5: Identify Coil Excitation Type
% ---------------------------------------------------------------
coil_types = arrayfun(@(x) x.excitation.TxRx, RLC_tmd, 'UniformOutput', false);
coil_types = string(coil_types);

coil_type_Rx                 = all(coil_types == "Rx");
coil_type_Tx                 = all(coil_types == "Tx");
coil_type_Tx_and_Rx          = all(ismember(coil_types, ["Rx", "Tx"])) && all(ismember(["Rx", "Tx"], coil_types));
coil_type_TxRx               = all(coil_types == "TxRx");
coil_type_TxRx_and_Rx        = all(ismember(coil_types, ["Rx", "TxRx"])) && all(ismember(["Rx", "TxRx"], coil_types));
coil_type_TxRx_and_Tx        = all(ismember(coil_types, ["TxRx", "Tx"])) && all(ismember(["TxRx", "Tx"], coil_types));
coil_type_TxRx_and_Tx_and_Rx = all(ismember(coil_types, ["TxRx", "Rx", "Tx"])) && all(ismember(["TxRx", "Rx", "Tx"], coil_types));

% ---------------------------------------------------------------
% STEP 6: Co-Simulation Optimization
% ---------------------------------------------------------------
fprintf('Performing Co-simulation ...\n');

if coil_type_Rx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_Rx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);

elseif coil_type_Tx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_Tx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);

elseif coil_type_Tx_and_Rx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_Rx_and_Tx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);

elseif coil_type_TxRx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_TxRx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);

elseif coil_type_TxRx_and_Rx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_TxRx_and_Rx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);

elseif coil_type_TxRx_and_Tx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_TxRx_and_Tx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);

elseif coil_type_TxRx_and_Tx_and_Rx
    [return_elements, M_cal, YPm, ZPm, SPm, SP_det, figure_freq] = ...
        co_simulation_TxRx_and_Rx_and_Tx_optim(YPn, emc, RLC_org, RLC_tmd, options_PSO_1, options_PSO_2);
end

% ---------------------------------------------------------------
% STEP 7: Update Lumped Elements
% ---------------------------------------------------------------
N_ports = geo_scoil_print_lumped_elements(coil_file, RLC_org, RLC_tmd, return_elements);

% ---------------------------------------------------------------
% STEP 8: Store Co-Simulation Results
% ---------------------------------------------------------------
MREDM.fields.netp.YP_s         = YPm;
MREDM.fields.netp.ZP_s         = ZPm;
MREDM.fields.netp.SP_s         = SPm;
MREDM.fields.netp.SP_p         = SP_det;
MREDM.fields.netp.figure_freq  = figure_freq;
MREDM.fields.M_cal             = M_cal;

fprintf('Storing Co-simulation Results ...\n');
solution = MREDM.fields;
save(fullfile('./data/solutions/', ...
     ['CoSimulation_', coil_file_inp(1:end-4), '.mat']), 'solution', '-v7.3');

fprintf('\n[MARIE 3.0] Standalone Co-simulation completed successfully.\n');
fprintf('[Saved] - ./data/solutions/CoSimulation_%s.mat\n', coil_file_inp(1:end-4));

% ============================= END OF FILE =============================
