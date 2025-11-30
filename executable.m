% ===============================================================
% FILE: executable.m
% PROJECT: MARIE 3.0 | Magnetic Resonance Integral Equation Suite
% YEAR: 2025
% LICENSE: MIT License
% ===============================================================
%
% PURPOSE:
%   Compiles a MARIE runner into a standalone executable using MATLAB Compiler
%
% SYNTAX:
%   executable(runner_name)
%
% INPUTS:
%   runner_name - (string) Name of the runner to compile (without .m extension)
%                 Examples: 'MARIE_runner', 'BASIS_runner', 'MRGF_runner', 'CO_SIMULATION_runner'
%
% OUTPUTS:
%   Compiled executable in ./bin/ directory
%
% EXAMPLES:
%   executable('MARIE_runner')
%   executable('BASIS_runner')
%
% ===============================================================

function executable(runner_name)

% Validate input
if nargin < 1 || isempty(runner_name)
    error('MARIE:NoRunner', 'Runner name is required.\nUsage: executable(''MARIE_runner'')');
end

% Remove .m extension if provided
if endsWith(runner_name, '.m')
    runner_name = runner_name(1:end-2);
end

% Define paths
runner_dir = './src/src_runners/';
runner_path = fullfile(runner_dir, [runner_name, '.m']);
output_dir = './bin';

% Get list of available runners
available_runners = dir(fullfile(runner_dir, '*_runner.m'));
runner_list = {available_runners.name};

% Check if runner file exists
if ~exist(runner_path, 'file')
    fprintf('ERROR: Runner file not found: %s\n', runner_path);
    fprintf('\nAvailable runners in %s:\n', runner_dir);
    for i = 1:length(runner_list)
        fprintf('  - %s\n', runner_list{i}(1:end-2));  % Remove .m extension
    end
    error('MARIE:RunnerNotFound', 'Runner "%s" not found.', runner_name);
end

fprintf('===============================================================\n');
fprintf('Compiling %s executable...\n', runner_name);
fprintf('===============================================================\n');

try
    fprintf('  Adding source directories to path...\n');
    addpath(genpath('./src'));

    fprintf('  Compiling %s...\n', runner_name);
    fprintf('  This may take several minutes...\n\n');

    % Compile the runner
    % Note: -R -singleCompThread forces serial execution in deployed mode
    mcc('-m', ...
        runner_path, ...
        '-o', runner_name, ...
        '-d', output_dir, ...
        '-a', './src', ...
        '-v');

    fprintf('\n===============================================================\n');
    fprintf('SUCCESS: Executable compiled successfully!\n');
    fprintf('===============================================================\n');
    fprintf('  Executable: %s\n', fullfile(output_dir, runner_name));
    fprintf('  Runner script: %s\n', fullfile(output_dir, ['run_' runner_name '.sh']));
    fprintf('\n  To run (with MATLAB installed):\n');
    fprintf('    %s\n', fullfile(output_dir, runner_name));
    fprintf('\n  To run (with MATLAB Runtime):\n');
    fprintf('    %s <mcr_path> <args>\n', fullfile(output_dir, ['run_' runner_name '.sh']));
    fprintf('===============================================================\n');

catch ME
    fprintf('\n===============================================================\n');
    fprintf('ERROR: Failed to compile executable\n');
    fprintf('===============================================================\n');
    fprintf('  Message: %s\n', ME.message);
    fprintf('  Error details:\n%s\n', ME.getReport());
    fprintf('===============================================================\n');
    error('MARIE:CompilationFailed', 'Compilation failed. See error details above.');
end

end