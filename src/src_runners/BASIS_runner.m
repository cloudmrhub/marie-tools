% ===============================================================
% FILE: BASIS_runner.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Generates and stores the electromagnetic basis used by the
%   MRGF solver in MARIE 3.0. Executes a full VIE-based simulation for a
%   dipole- or surface-supported excitation basis.
%
% SYNTAX:
%   Run directly in MATLAB:
%       >> BASIS_runner
%
% DESCRIPTION:
%   This script executes all major stages of the BASIS generation workflow:
%       1. Load user-defined input file
%       2. Build geometry and discretizations
%       3. Assemble VIE operators
%       4. Generate the incident field basis
%       5. Solve the VIE for each basis vector
%       6. Compute total EM field basis
%       7. Store basis and MRGF operators for future simulations
%
% VERSION:
%   MARIE 3.0 (2025)
%
% NOTES:
%   - The input file specifies:
%     1. the main magnetic field strength (mandatory)
%     2. the VIE basis functions (PWC or PWL) (mandatory)
%     3. the body model (mandatory)
%     4. the coil model (ignored)
%     5. the wire model (ignored)
%     5. the shield model (ignored)
%     6. the substrate that acts as support of the EM basis (only used if 
%     geo_flag = 2)
%     7. the co-simulation activation flag (ignored)
%     8. the basis file to store the EM basis (mandatory)
%   - The BASIS generator produces the fields and MRGF operators used for
%     MRGF-accelerated simulations.
%   - The variable 'geo_flag' selects the basis type:
%       2 - Surface-supported basis
%       3 - Dipole-supported basis
%   - The stored .mat file includes:
%       Total Electric Field Basis
%       Total Magnetic Field Basis
%       MRGF operators
%       UISNR
%       UITXE
%       Ideal Current Patterns only for surface-supported basis

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
geo_flag   = 2;                                 % 2 = Surface-supported basis

% ---------------------------------------------------------------
% STEP 2: Load Inputs
% ---------------------------------------------------------------
fprintf('Loading inputs ...\n');
MREDM = load_inputs(input_file);  % Load simulation parameters

% ---------------------------------------------------------------
% STEP 3: Geometry Construction
% ---------------------------------------------------------------
fprintf('Geometry Construction ...\n');
MREDM = geo_assembly(MREDM, geo_flag);  % Build geometry for basis support

% ---------------------------------------------------------------
% STEP 4: Assemble VIE Operators
% ---------------------------------------------------------------
fprintf('Assembling VIE Operators ...\n');
MREDM = vie_assembly(MREDM, geo_flag);  % Assemble VIE operators

% ---------------------------------------------------------------
% STEP 5: Generate Incident Field Basis
% ---------------------------------------------------------------
fprintf('Generating Incident Field Basis ...\n');
BASIS = IncidentField_Basis(MREDM, geo_flag);  % Create EM basis excitations

if geo_flag == 2
    BASIS.Substrate = MREDM.SIE.coil;  % Store substrate for surface-supported basis
end

% ---------------------------------------------------------------
% STEP 6: Solve VIE System
% ---------------------------------------------------------------
fprintf('Solving VIE ...\n');
BASIS = ie_solver_vie_basis(MREDM, BASIS);  % Solve for basis current distributions

% ---------------------------------------------------------------
% STEP 7: Compute Total Field Basis
% ---------------------------------------------------------------
fprintf('Generating Total Field Basis ...\n');
BASIS = em_ehfield_vie(MREDM, BASIS);  % Compute total EM field basis

% ---------------------------------------------------------------
% STEP 8: Store Basis
% ---------------------------------------------------------------
fprintf('Storing Basis ...\n');
save(MREDM.inputs.basis_file, 'BASIS', '-v7.3');

fprintf('\n[MARIE 3.0] Basis generation completed successfully.\n');
fprintf('[Saved] - %s\n', MREDM.inputs.basis_file);

% ============================= END OF FILE =============================
