% ===============================================================
% FILE: DLDG_runner.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Execution script for a series of WSVIE simulations with MARIE 3.0.
%   sharing the same coil configurations and the same VIE domain.
%   This runner executes the full WSVIE simulation and co-simulation.
%   All body models should be placed in the exact same domain.
%
% SYNTAX:
%   Run directly in MATLAB:
%       >> DLDG_runner
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
%       8. Repeats the WSVIE solution for multiple body models
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
%   - The user needs to specify the directory where the solution will be
%   stored
%   - The user needs to specify the directory with the body models

% ---------------------------------------------------------------
% STEP 0: Initialize Environment
% ---------------------------------------------------------------
close all
clearvars
clc

% ---------------------------------------------------------------
% STEP 1: Specify Input File
% ---------------------------------------------------------------
input_file = 'inp_PregnantModel_SKYRA.json';  % Input file name

% ---------------------------------------------------------------
% STEP 2: Specify Directory for Body Models and Solutions
% ---------------------------------------------------------------
body_model_directory = './data/bodies/RHBM/4mm/Pregnant_Fetal_Models/Model_v1/';
solution_directory   = './data/solutions/Pregnant_Fetal_Model_v1/';
all_body_files       = dir(fullfile(body_model_directory, '*.mat'));

% ---------------------------------------------------------------
% STEP 3: Load Inputs
% ---------------------------------------------------------------
fprintf('Loading inputs ...\n');
MREDM = load_inputs(input_file); % Load simulation parameters

% ---------------------------------------------------------------
% STEP 4: Update to Master Domain
% ---------------------------------------------------------------
fprintf('Create master domain ...\n');
MREDM = update_RHBM(MREDM); % Build master domain

% ---------------------------------------------------------------
% STEP 5: Geometry Construction
% ---------------------------------------------------------------
fprintf('Geometry Construction ...\n');
MREDM = geo_assembly(MREDM, 1); % Build coil/body discretizations

% ---------------------------------------------------------------
% STEP 6: Assemble Integral Equation Operators
% ---------------------------------------------------------------
fprintf('Assembling Integral Equation Operators ...\n');
MREDM = wsvie_assembly(MREDM, 1); % Assemble WSVIE system matrices

% ---------------------------------------------------------------
% STEP 7: Assemble Coil Preconditioners
% ---------------------------------------------------------------
fprintf('Assembling Coil Preconditioners ...\n');
MREDM = prec_wsvie(MREDM); % Assemble WSVIE system preconditioners

% ---------------------------------------------------------------
% STEP 8: Set Master Operatros
% ---------------------------------------------------------------
fprintf('Store Master Operators ...\n');
PS_master      = MREDM.operators.pfft_PS;
Z_B2C_N_master = MREDM.operators.pfft_Z_bc_N;
Z_B2C_K_master = MREDM.operators.pfft_Z_bc_K;

% ---------------------------------------------------------------
% STEP 9: Repeat WSVIE Simulations per Body
% ---------------------------------------------------------------
fprintf('Solving WSVIE ...\n\n');
for i = 1:numel(all_body_files)
    
    fprintf('\tBody Model %i ...',i);

    % =========================== START OF LOOP ===========================
    % Update body from the train of files
    MREDM.inputs.rhbm = fullfile(all_body_files(i).folder, all_body_files(i).name);

    if strcmp(all_body_files(i).name,'Master_Domain.mat')
        fprintf(' Skipping Master Domain.\n');
        continue;
    end
    if exist(strcat(solution_directory,all_body_files(i).name),'file')
        fprintf(' Solution already exists.\n');
        continue;
    end
    fprintf('\n');
    
    % Update EP and VIE Domain Parameters
    [MREDM,N_coil] = geo_recompute_domain(MREDM,1);
    
    % Get dimensions
    ql     = MREDM.dimensions.ql;
    N_scat = MREDM.dimensions.N_scat;

    % Get PFFT Body-Specific Operators
    if ~isempty(PS_master)
        idx_solution                = MREDM.dimensions.idx_solution;
        pfft_idx_domain             = MREDM.dimensions.pfft_idx_domain;
        MREDM.operators.pfft_Z_bc_K = Z_B2C_K_master(idx_solution,:);
        MREDM.operators.pfft_Z_bc_N = Z_B2C_N_master(idx_solution,:);
        MREDM.operators.pfft_PS     = PS_master(pfft_idx_domain,[1:N_coil N_coil+idx_solution']);
    end
    
    % Update RHS
    MREDM.solver.rhs = rhs_assembly(MREDM,N_scat*ql); 
    
    % Update VIE preconditioner
    MREDM = prec_vie(MREDM);          

    % Perform wsvie inversion
    try 
        MREDM.fields.Jcb = MREDM.functions.ie_solver_wsvie(MREDM,0);  
    catch
        warning('Out of GPU memory. Running in CPU.');
        MREDM.fields.Jcb = MREDM.functions.ie_solver_wsvie(MREDM,5);  
    end

    % Compute Network Parameters
    MREDM = np_compute(MREDM);  

    % Perform circuit co-simulation
    MREDM = co_simulation(MREDM);      

    % Compute EM fields and SNR
    MREDM = em_ehfield_wsvie(MREDM); 
    % ============================ END OF LOOP ============================

    % Save fields
    solution = MREDM.fields;                                                    
    save(strcat(solution_directory,all_body_files(i).name), 'solution', '-v7.3');
end
% ============================= END OF FILE =============================