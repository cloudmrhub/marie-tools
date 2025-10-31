clc
clearvars
close all

FLAGS=cell(0,1);
FLAGS{end+1}='-funroll-loops -fopenmp -std=c++17 -fPIC -lgomp';
FLAGS{end+1}='-fconstexpr-depth=100000 -ftemplate-depth=100000';
FLAGS{end+1}='--param large-function-growth=2000 --param inline-unit-growth=2000 -finline-limit=3000';
FLAGS{end+1}='-fno-signed-zeros -fno-signaling-nans -fno-trapping-math -fassociative-math';
FLAGS{end+1}='-Wall -Wopenmp-simd -Wvector-operation-performance -Winline';
FLAGS{end+1}='-Wunused-const-variable';
FLAGS{end+1}='-march=native -O3 -DUSE_SIMD -DASSEMBLE_BOTH';
FLAGS=strjoin(FLAGS);

CXXFLAGS      = ['CXXFLAGS="',FLAGS,'"'];
LDFLAGS       = ['LDFLAGS="',FLAGS,' -t'];
LDFLAGS       = [LDFLAGS,'"'];
CXXOPTIMFLAGS ='CXXOPTIMFLAGS="-DNDEBUG"';

% WVIE PATHS
SRCDIR_WVIE   = './src/src_integral_equations/src_wvie/Cpp_Assembly/';
BINDIR_WVIE   = fullfile(SRCDIR_WVIE,'bin');
IPATH_WVIE    = fullfile(SRCDIR_WVIE,'include');
% SVIE PATHS
SRCDIR_SVIE   = './src/src_integral_equations/src_svie/Cpp_Assembly/';
BINDIR_SVIE   = fullfile(SRCDIR_SVIE,'bin');
IPATH_SVIE    = fullfile(SRCDIR_SVIE,'include');
% SIE PATHS
SRCDIR_SIE    = './src/src_integral_equations/src_sie/src_operators_sie/Cpp_Assembly/'; 
BINDIR_SIE    = fullfile(SRCDIR_SIE,'bin');
IPATH_SIE     = fullfile(SRCDIR_SIE,'include');
% VIE PATHS
SRCDIR_VIE    = './src/src_integral_equations/src_vie/src_operators_vie/Cpp_Assembly/'; 
BINDIR_VIE    = fullfile(SRCDIR_VIE,'bin');
IPATH_VIE     = fullfile(SRCDIR_VIE,'include');
% mex Arguments
mexArgs_wvie  = {CXXFLAGS,CXXOPTIMFLAGS,LDFLAGS,'-outdir',BINDIR_WVIE,['-I',IPATH_WVIE]};
mexArgs_svie  = {CXXFLAGS,CXXOPTIMFLAGS,LDFLAGS,'-outdir',BINDIR_SVIE,['-I',IPATH_SVIE]};
mexArgs_sie   = {CXXFLAGS,CXXOPTIMFLAGS,LDFLAGS,'-outdir',BINDIR_SIE,['-I',IPATH_SIE]};

%% Compile Volume-Wire Integral Equation Operators
files      = dir(fullfile(strcat(SRCDIR_WVIE,'src/'), '*.*'));
fileNames  = {files.name}';
fileNames  = fileNames(~ismember(fileNames, {'.', '..'}));
mexFiles   = fileNames(startsWith(fileNames, 'mex'));  
otherFiles = fileNames(~startsWith(fileNames, 'mex'));
for i = 1:length(otherFiles)
    fileName = otherFiles{i}; 
    [~,name,~] = fileparts(fileName); 
    mex(mexArgs_wvie{:},'-R2018a', '-output', fullfile(SRCDIR_WVIE,'src/',name), fullfile(SRCDIR_WVIE,'src/',mexFiles{i}), fullfile(SRCDIR_WVIE,'src/', otherFiles{i}));
end

%% Compile Volume-Surface Integral Equation Operators
files      = dir(fullfile(strcat(SRCDIR_SVIE,'src/'), '*.*'));
fileNames  = {files.name}';
fileNames  = fileNames(~ismember(fileNames, {'.', '..'}));
mexFiles   = fileNames(startsWith(fileNames, 'mex'));  
otherFiles = fileNames(~startsWith(fileNames, 'mex'));
for i = 1:length(otherFiles)
    fileName = otherFiles{i}; 
    [~,name,~] = fileparts(fileName); 
    mex(mexArgs_svie{:},'-R2018a', '-output', fullfile(SRCDIR_SVIE,'src/',name), fullfile(SRCDIR_SVIE,'src/',mexFiles{i}), fullfile(SRCDIR_SVIE,'src/', otherFiles{i}));
end

%% Complile Surface Integral Equation Operators
files      = dir(fullfile(strcat(SRCDIR_SIE,'src/'), '*.*'));
fileNames  = {files.name}';
fileNames  = fileNames(~ismember(fileNames, {'.', '..'}));
mexFiles   = fileNames(endsWith(fileNames, '_mex.cpp'));  
otherFiles = fileNames(~endsWith(fileNames, '_mex.cpp'));
for i = 1:length(otherFiles)
    fileName = otherFiles{i}; 
    [~,name,~] = fileparts(fileName); 
    mex(mexArgs_sie{:},'-R2017b', '-output', fullfile(SRCDIR_SIE,'src/',name), fullfile(SRCDIR_SIE,'src/',mexFiles{i}), fullfile(SRCDIR_SIE,'src/', otherFiles{i}));
end

%% Complile Volume Integral Equation Operators
str_mex_command_with_opts = ['mex "-I' IPATH_VIE '" -c -outdir "' SRCDIR_VIE '"'];
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/GL_1D.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_EA/A_functions_ws_ea.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_EA/N_functions_ws_ea.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_EA/PSI_limits_ws_ea.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_EA/Simplex_ws_ea.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_EA/THETA_limits_ws_ea.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_EA/quadric_ws_ea.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/a_functions_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/lambda_limits_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/n_functions_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/psi_limits_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/quadric_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/subtriangles_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_ST/u_limits_ws_st.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_VA/rho_limits_ws_va.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_VA/theta_p_limits_ws_va.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_VA/theta_q_limits_ws_va.cpp')]);
eval([str_mex_command_with_opts, ' ', fullfile(SRCDIR_VIE,'src/src_cpp/WS_VA/quadric_ws_va.cpp')]);
s1 = strcat(SRCDIR_VIE,'GL_1D.o');
s2 = strcat(SRCDIR_VIE,'A_functions_ws_ea.o');
s3 = strcat(SRCDIR_VIE,'N_functions_ws_ea.o');
s4 = strcat(SRCDIR_VIE,'PSI_limits_ws_ea.o');
s5 = strcat(SRCDIR_VIE,'Simplex_ws_ea.o');
s6 = strcat(SRCDIR_VIE,'THETA_limits_ws_ea.o');
s7 = strcat(SRCDIR_VIE,'quadric_ws_ea.o');
s8 = strcat(SRCDIR_VIE,'a_functions_ws_st.o');
s9 = strcat(SRCDIR_VIE,'lambda_limits_ws_st.o');
s10 = strcat(SRCDIR_VIE,'n_functions_ws_st.o');
s11 = strcat(SRCDIR_VIE,'psi_limits_ws_st.o');
s12 = strcat(SRCDIR_VIE,'quadric_ws_st.o');
s13 = strcat(SRCDIR_VIE,'subtriangles_ws_st.o');
s14 = strcat(SRCDIR_VIE,'u_limits_ws_st.o');
s15 = strcat(SRCDIR_VIE,'rho_limits_ws_va.o');
s16 = strcat(SRCDIR_VIE,'theta_p_limits_ws_va.o');
s17 = strcat(SRCDIR_VIE,'theta_q_limits_ws_va.o');
s18 = strcat(SRCDIR_VIE,'quadric_ws_va.o');
str_obj = [' ' s1 ' ' s2 ' ' s3 ' ' s4 ' ' s5 ' ' s6 ' ' s7 ' ' s8 ' ' s9 ' ' s10 ' ' s11 ' ' s12 ' ' s13 ' ' s14 ' ' s15 ' ' s16 ' ' s17 ' ' s18];
mexArgs_vie = [CXXFLAGS, CXXOPTIMFLAGS, LDFLAGS, '-outdir', BINDIR_VIE, ['-I', IPATH_VIE], strsplit(str_obj, ' '), fullfile(SRCDIR_VIE, 'src/src_cpp/Kernels.cpp')];
mex(mexArgs_vie{:}, '-R2017b', '-output', fullfile(SRCDIR_VIE, 'src/src_cpp/solve_ea'), fullfile(SRCDIR_VIE, 'src/solve_ea_mex.cpp'), fullfile(SRCDIR_VIE, 'src/src_cpp/create_EA.cpp'));
mex(mexArgs_vie{:}, '-R2017b', '-output', fullfile(SRCDIR_VIE, 'src/src_cpp/solve_st'), fullfile(SRCDIR_VIE, 'src/solve_st_mex.cpp'), fullfile(SRCDIR_VIE, 'src/src_cpp/create_ST.cpp'));
mex(mexArgs_vie{:}, '-R2017b', '-output', fullfile(SRCDIR_VIE, 'src/src_cpp/solve_va'), fullfile(SRCDIR_VIE, 'src/solve_va_mex.cpp'), fullfile(SRCDIR_VIE, 'src/src_cpp/create_VA.cpp'));
delete(fullfile(SRCDIR_VIE, '*.o'))

clc
clearvars

fprintf('All C++ files of the MARIE suite are installed. \nHappy simulations!\n')