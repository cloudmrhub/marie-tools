% ===============================================================
% FILE: load_inputs.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Reads the input JSON file, parses simulation parameters, initializes
%   tolerances and quadrature settings, and sets up the MREDM structure.
%
% SYNTAX:
%   MREDM = load_inputs(input_file)
%
% INPUTS:
%   input_file - (string) Path to the JSON input file located in ./data/inputs/
%
% OUTPUTS:
%   MREDM      - Struct containing simulation parameters.
%
% DESCRIPTION:
%   This function reads the user-defined input file (JSON format) and
%   initializes the primary simulation structure (MREDM) with:
%       - Model configuration (coil, body, wire, shield, etc.)
%       - EM constants (from B0)
%       - Numerical tolerances and quadrature settings
%       - Function handles for geometry, assembly, and solver routines
%       - Optional precomputed EM basis paths (if defined)
%
% VERSION:
%   MARIE 3.0 (2025)
%
% NOTES:
%   - Automatically starts a parallel pool if none is active.
%   - The input JSON must define fields:
%       B0, Basis_Functions_VIE, BodyFile, CoilFile, WireFile,
%       ShieldFile, SurfaceBasisSupportFile, TMD, BasisFile.
%   - All paths are relative to ./data/
%
% ===============================================================

function [MREDM] = load_inputs(input_file)
    
    if isempty(gcp('nocreate'))
        parpool('IdleTimeout', 240);  
    end
    
    fid  = fopen(input_file, 'r');
    raw  = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    inputs_json = jsondecode(raw);
    
    inputs.B0            = double(inputs_json.B0);                  % Magnetic field strength (T)
    inputs.nucleus       = inputs_json.Nucleus;                     % Magnetic field strength (T)
    inputs.PWX           = double(inputs_json.Basis_Functions_VIE); % Basis type (PWC/PWL)
    inputs.rhbm          = inputs_json.BodyFile;                    % Body model
    inputs.coil          = inputs_json.CoilFile;                    % Coil geometry
    inputs.wire          = inputs_json.WireFile;                    % Wire geometry
    inputs.shield        = inputs_json.ShieldFile;                  % Shield geometry
    inputs.basis_support = inputs_json.SurfaceBasisSupportFile;     % Basis surface support (if any)
    inputs.tmd           = double(inputs_json.TMD);                 % Co-simulation flag
    inputs.basis_file    = inputs_json.BasisFile;                   % Optional external EM basis

    if ~isempty(inputs.basis_file)
        inputs.basis_file = fullfile('./data/bases/', inputs.basis_file);
    end
    
    inputs.tol              = 1e-5;             % GMRES
    inputs.tol_HOSVD        = inputs.tol*1e-2;  % HOSVD on Zbb (lower than GMRES)
    inputs.tol_HOSVD_couple = 1e-12;            % HOSVD on Zbc (much lower than GMRES)
    inputs.tol_TT           = inputs.tol*1e+2;  % TT on Zbc (1e-3 is fine for most applications)
    inputs.tol_ACA          = inputs.tol;       % As accurate as GMRES (Lower is fine but needs more time and memory)
    inputs.tol_rSVD         = inputs.tol*1e+2;  % rSVD: Lower than GMRES
    inputs.tol_sSVD         = inputs.tol*1e+2;  % sSVD: Lower than GMRES
    inputs.tol_DEIM         = inputs.tol*1e+1;  % DEIM: between rSVD/sSVD and GMRES
    
    inputs.Quad_order_wie      = 6;  % WIE
    inputs.Np_1D_far_V         = 4;  % VIE non-self term
    inputs.Np_1D_medium_V      = 8;  % VIE adjacent terms
    inputs.Np_1D_near_S        = 15; % VIE self-terms
    inputs.Quad_order_wie_coup = 4;  % WVIE and WSIE
    inputs.Quad_order_sie      = 5;  % SVIE and WSIE
    inputs.Quad_order_vie      = 2;  % WVIE and SVIE
    inputs.N_NS                = 4;  % WVIE and SVIE
    inputs.N_ST_psi            = 10; % SIE self-terms
    inputs.N_EA_theta          = 6;  % SIE edge-terms
    inputs.N_EA_psi            = 6;  % SIE edge-terms
    inputs.N_VA_theta_p        = 6;  % SIE vertex-terms
    inputs.N_VA_theta_q        = 6;  % SIE vertex-terms
    inputs.N_VA_psi            = 6;  % SIE vertex-terms
    
    inputs.Basis_distance  = 0.02;   % Offset distance for dipole basis support
    inputs.Basis_thickness = 3;      % Number of voxel layer thickness
    inputs.rSVD_blocksize  = 1000;   % Block size for randomized SVD
    
    inputs.pFFT_flag = 1;  % Enable precorrected FFT acceleration
    
    MREDM.inputs    = inputs;                   % Store all parameters
    MREDM.emc       = em_constants(inputs.B0,inputs.nucleus);  % Electromagnetic constants
    
end
