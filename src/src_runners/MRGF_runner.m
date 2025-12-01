% ===============================================================
% FILE: MRGF_runner.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Main execution script for the MRGF workflow in MARIE 3.0.
%   Executes the full WSVIE simulation using the MRGF compressed 
%   operators and precomputed EM basis.
%
% SYNTAX:
%   Run directly in MATLAB:
%       >> MRGF_runner
%
% DESCRIPTION:
%   This script sequentially executes all major stages of the MRGF workflow:
%       1. Load user-defined input file
%       2. Build geometry and discretizations
%       3. Load precomputed basis functions and MRGF operators
%       4. Assemble compressed WSVIE operators
%       5. Solve system of equations
%       6. Perform co-simulation (Tuning/Matching/Decoupling/Detuning)
%          to update variable lumped elements and matching networks
%       7. Compute electromagnetic fields
%       8. Store simulation results
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
%     at least one of the coil/wire/shield models is needed
%     6. the substrate that acts as support of the EM basis (ignored)
%     7. the co-simulation activation flag
%     8. the basis file to store or load the EM basis and the MRGF
%   - The final solution is stored under ./data/solutions/
%   - The MRGF runner assumes that the EM basis and operator compression
%     matrices are already generated and stored in the input .mat file.

function MRGF_runner(input_file)

% ---------------------------------------------------------------
% STEP 0: Initialize Environment
% ---------------------------------------------------------------
close all
clc

% ---------------------------------------------------------------
% STEP 1: Specify Input File
% ---------------------------------------------------------------
% Parse command-line arguments
if nargin < 1 || ~isfile(input_file)
    fprintf('[ERROR] Usage: MARIE_runner(''input_file.json'')\n');
    fprintf('[ERROR] Example: MARIE_runner(''./data/inputs/inp_Duke_StadiumTriangular.json'')\n');
    error('MARIE:FileNotFound', 'Valid input file is required. Exiting.');
end

% ---------------------------------------------------------------
% STEP 2: Load Inputs
% ---------------------------------------------------------------
fprintf('Loading inputs ...\n');
MREDM = load_inputs(input_file);  % Load simulation parameters

% ---------------------------------------------------------------
% STEP 3: Geometry Construction
% ---------------------------------------------------------------
fprintf('Geometry Construction ...\n');
MREDM = geo_assembly(MREDM, 4);  % Geometry setup for MRGF mode

% ---------------------------------------------------------------
% STEP 4: Load Basis Functions and Compressed Operators
% ---------------------------------------------------------------
fprintf('Reading compressed EM basis ...\n');

Ue        = h5read(MREDM.inputs.basis_file, '/BASIS/Ue');
Ub        = h5read(MREDM.inputs.basis_file, '/BASIS/Ub');
U_hat_inv = h5read(MREDM.inputs.basis_file, '/BASIS/U_hat_inv');
S_hat_inv = h5read(MREDM.inputs.basis_file, '/BASIS/S_hat_inv');
V_hat_inv = h5read(MREDM.inputs.basis_file, '/BASIS/V_hat_inv');
X         = h5read(MREDM.inputs.basis_file, '/BASIS/X');

% Convert from real/imaginary structs to complex matrices
Ue        = Ue.real + 1i * Ue.imag;
Ub        = Ub.real + 1i * Ub.imag;
U_hat_inv = U_hat_inv.real + 1i * U_hat_inv.imag;
V_hat_inv = V_hat_inv.real + 1i * V_hat_inv.imag;
X         = X.real + 1i * X.imag;

% ---------------------------------------------------------------
% STEP 5: Assemble Integral Equation Operators
% ---------------------------------------------------------------
fprintf('Assembling Integral Equation Operators ...\n');
MREDM = wsvie_assembly_mrgf(MREDM);  % Assemble WSVIE matrices using MRGF

% ---------------------------------------------------------------
% STEP 6: System Solution
% ---------------------------------------------------------------
fprintf('Solving System ...\n');
MREDM = ie_solver_svie_mrgf(MREDM, U_hat_inv, S_hat_inv, V_hat_inv, X); 
% Solve linear system directly

% ---------------------------------------------------------------
% STEP 7: Co-Simulation
% ---------------------------------------------------------------
fprintf('Performing Co-simulation ...\n');
MREDM = co_simulation(MREDM);  % Perform circuit co-simulation

% ---------------------------------------------------------------
% STEP 8: Compute Electromagnetic Fields
% ---------------------------------------------------------------
fprintf('Computing electromagnetic fields ...\n');
MREDM = em_ehfield_wsvie_mrgf(MREDM, Ue, Ub);  % Compute EM fields, SNR, TXE

% ---------------------------------------------------------------
% STEP 9: Store Results
% ---------------------------------------------------------------
fprintf('Storing Fields ...\n');
solution = MREDM.fields;

[~, file_name, ext] = fileparts(input_file);
save(fullfile('./data/solutions/', ...
     ['MRGF_', file_name, '.mat']), 'solution', '-v7.3');

fprintf('\n[MARIE 3.0] MRGF simulation completed successfully.\n');
fprintf('[Saved] - ./data/solutions/MRGF_%s.mat\n', file_name);

end
% ============================= END OF FILE =============================