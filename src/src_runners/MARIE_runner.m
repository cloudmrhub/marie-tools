% ===============================================================
% FILE: MARIE_runner.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Main execution script for the WSVIE in MARIE 3.0. This runner executes 
%   the full WSVIE simulation and co-simulation for the input MRI setup.
%
% SYNTAX:
%   Run directly in MATLAB:
%       >> MARIE_runner
%
% DESCRIPTION:
%   This script sequentially executes all major stages of the MARIE workflow:
%       1. Load user-defined input file
%       2. Build geometry and discretizations
%       3. Assemble WSVIE operators
%       4. Solve system of equations
%       5. Perform co-simulation (Tuning/Matching/Decoupling/Detuning)
%          to update variable lumped elements and matching networks
%       6. Compute electromagnetic fields
%       7. Store simulation results
%
% VERSION:
%   MARIE 3.0 (2025)
%
% NOTES:
%   - The input file specifies:
%     1. the main magnetic field strength (mandatory)
%     2. the VIE basis functions (PWC or PWL) (mandatory)
%     3. the body model (mandatory)
%     4. the coil model (can be left empty)
%     5. the wire model (can be left empty)
%     5. the shield model (can be left empty)
%     6. the substrate that acts as support of the EM basis (ignored)
%     7. the co-simulation activation flag
%     8. the basis file to store or load the EM basis and the MRGF 
%     ignored if a coil, wire, or shield model is provided, 
%     otherwise it is used as external excitation
%   - The final solution is stored under ./data/solutions/

% ---------------------------------------------------------------
% STEP 0: Initialize Environment
% ---------------------------------------------------------------
close all
clearvars
clc

% ---------------------------------------------------------------
% STEP 1: Specify Input File
% ---------------------------------------------------------------
input_file = 'inp_Duke_StadiumTriangular.json';  % Input file name

% ---------------------------------------------------------------
% STEP 2: Load Inputs
% ---------------------------------------------------------------
fprintf('Loading inputs ...\n');
MREDM = load_inputs(input_file);  % Load simulation parameters

% ---------------------------------------------------------------
% STEP 3: Geometry Construction
% ---------------------------------------------------------------
fprintf('Geometry Construction ...\n');
MREDM = geo_assembly(MREDM, 1);  % Build coil/body discretizations

% ---------------------------------------------------------------
% STEP 4: Assemble Integral Equation Operators
% ---------------------------------------------------------------
fprintf('Assembling Integral Equation Operators ...\n');
MREDM = wsvie_assembly(MREDM, 1);  % Assemble WSVIE system matrices

% ---------------------------------------------------------------
% STEP 5: System Solution
% ---------------------------------------------------------------
fprintf('Solving System ...\n');
MREDM = solver_wsvie(MREDM);  % Solve linear system iteratively

% ---------------------------------------------------------------
% STEP 6: Co-Simulation
% ---------------------------------------------------------------
fprintf('Performing Co-simulation ...\n');
MREDM = co_simulation(MREDM);  % Perform circuit co-simulation

% ---------------------------------------------------------------
% STEP 7: Compute Electromagnetic Fields
% ---------------------------------------------------------------
fprintf('Computing electromagnetic fields ...\n');
MREDM = em_ehfield_wsvie(MREDM);  % Compute EM fields, SNR, and TXE

% ---------------------------------------------------------------
% STEP 8: Store Results
% ---------------------------------------------------------------
fprintf('Storing Fields ...\n');
solution = MREDM.fields;

save(fullfile('./data/solutions/', ...
     ['MARIE_', num2str(MREDM.inputs.tmd), '_', input_file(1:end-4), '.mat']), ...
     'solution', '-v7.3');

fprintf('\n[MARIE 3.0] Simulation completed successfully.\n');
fprintf('[Saved] - ./data/solutions/MARIE_%d_%s.mat\n', ...
        MREDM.inputs.tmd, input_file(1:end-4));

% ============================= END OF FILE =============================