% ===============================================================
% FILE: em_constants.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Generates fundamental electromagnetic and MRI-related physical
%   constants for a given static magnetic field B0.
%
% SYNTAX:
%   emc = em_constants(B0)
%
% INPUTS:
%   B0  - Static magnetic field strength [Tesla]
%
% OUTPUTS:
%   emc - Struct containing electromagnetic and MRI constants:
%         - Fundamental physical constants 
%         - Derived EM quantities
%         - Material parameters (
%         - MRI-specific constants 
%
% DESCRIPTION:
%   Used throughout MARIE to compute electromagnetic quantities
%   and normalization factors that depend on the main magnetic
%   field strength. Includes both universal and MRI-specific
%   constants (e.g., gyromagnetic ratio, magnetization).
%
% VERSION:
%   MARIE 3.0 (2025)
%
% NOTES:
%   - Magnetization is computed assuming pure water proton density.
%
% ===============================================================

function [emc] = em_constants(B0)

    emc.gamma0      = 42577478.518;        % Gyromagnetic ratio / (2*pi) [Hz/T]
    emc.gamma0_2pi  = emc.gamma0 * 2 * pi; % Gyromagnetic ratio [rad/T*s]
    emc.mu0         = 4e-7 * pi;           % Magnetic permeability of free space [H/m]
    emc.c0          = 299792458;           % Speed of light [m/s]
    emc.Nproton     = 6.691e28;            % Proton density in water [1/m^3]
    emc.h           = 6.626e-34;           % Planck constant [J*s]
    emc.Ispin       = 1/2;                 % Proton spin
    emc.k_B         = 1.3806503e-23;       % Boltzmann constant [J/K]
    emc.se_copper   = 5.96e7;              % Conductivity of copper [S/m]
    emc.re_wire     = 1.2579e-08;          % Resistivity of copper [\Omega*m]
    emc.T0          = 310;                 % Body temperature [K]

    emc.freq   = emc.gamma0 * B0;                 % Larmor frequency [Hz]
    emc.e0     = 1 / (emc.c0^2 * emc.mu0);        % Permittivity of free space [F/m]
    emc.eta0   = emc.mu0 * emc.c0;                % Impedance of free space [\Omega]
    emc.omega  = 2 * pi * emc.freq;               % Angular frequency [rad/s]
    emc.lamda  = emc.c0 / emc.freq;               % Free-space wavelength [m]
    emc.k0     = 2 * pi / emc.lamda;              % Wavenumber [rad/m]

    emc.thick_wire = 0.0005 * 2;                                          % Wire thickness [m]
    emc.skin_depth = sqrt(2) / sqrt(emc.omega * emc.mu0 * emc.se_copper); % Skin depth [m]
    emc.rho_s      = 1 / (emc.se_copper * emc.skin_depth);                % Surface resistance [\Omega]

    emc.hbar       = emc.h / (2 * pi);                                % Reduced Planck constant [J*s]
    emc.ce         = 1i * emc.omega * emc.e0;                         % Electric-field scaling constant
    emc.cm         = 1i * emc.omega * emc.mu0;                        % Magnetic-field scaling constant
    emc.Z0         = 50;                                              % Characteristic impedance [\Omega]
    emc.Preamp_res = 1500;                                            % Preamplifier decoupling resistance [\Omega]
    emc.Detune_res = 1e5;                                             % Detuning resistance [\Omega]

    % Magnetization per unit volume (water protons, thermal equilibrium)
    emc.M0 = emc.Nproton * ((emc.gamma0_2pi * emc.hbar)^2) * ...
              emc.Ispin * (emc.Ispin + 1) * B0 / (3 * emc.k_B * emc.T0);

end
