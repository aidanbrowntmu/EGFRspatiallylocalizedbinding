clc
clear
close all

% perform PRCC on both the localized and uniform model for SS signaling
% October 3, 2026

% Define the uniform models DEs
R_uni_func  = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);

RL_uni_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);

Rc_uni_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

% Define the localized (tetraspanin) models DEs
R_loc_func  = @(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki) kp ./(kbi + kT + kL - ((kT_neg.*kT)./(kT_neg+kLT)) + ((kL_neg .*(((kLT.*kT) ./ (kT_neg+kLT)) +kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))));

T_loc_func  = @(R,kT,kT_neg,kLT) (kT.*R) ./ (kT_neg+kLT);

RL_loc_func = @(R,kT,kT_neg,kLT,kL,kL_neg,ki,kc,kc_neg) -R.*(((kLT.*kT) ./ (kT_neg+kLT)) +kL)./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);

Rc_loc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

% Define the parameters across both models

param_names = {'k_p', 'k_{bi}', 'k_{-L}', 'k_{+c}', 'k_{-c}', 'k_i', 'k_{+T}', 'k_{-T}', 'k_L', 'k_{L,slow}'};

num_params  = length(param_names);

% Default parameters (estimate of experimental conditions)
kp_base      = 1.0;
kbi_base     = 3*10^(-4);
kL_neg_base  = 0.05;
kc_base      = 0.1;
kc_neg_base  = 0.05;
ki_base      = 0.05;
kT_base      = 1/30;
kT_neg_base  = 0.12;
kL_base      = 10^-2;
kL_slow_base = 0.1 * kL_base;

% define vector of default parameters
p_base = [kp_base, kbi_base, kL_neg_base, kc_base, kc_neg_base, ki_base, kT_base, kT_neg_base, kL_base, kL_slow_base];

% Define upper and lower bounds for sampling parameter values
% Vary up and down by two orders of magnitude from default 
lower = p_base * 0.01;
upper = p_base * 100;

% Use latin hypercube function for random sampling 
N = 1000; % Have 1000 iterations
rng(1);  

% 
LHS = lhsdesign(N, num_params);
X   = lower + LHS .* (upper - lower); % Rescaled matrix (N x 10)

% Record steady state signaling levels for both models and the ratio of across both models
% initialize empty vectors 

Y_uniform   = zeros(N, 1); % Uniform steady state concentrations (RL+Rc)
Y_localized   = zeros(N, 1); % Localized steady state concentration (RL+Rc)
Y_ratio = zeros(N, 1); % Ratio of Localized/Uniform 

for i = 1:N

    % iterate through each parameter
    kp      = X(i, 1);
    kbi     = X(i, 2);
    kL_neg  = X(i, 3);
    kc      = X(i, 4);
    kc_neg  = X(i, 5);
    ki      = X(i, 6);
    kT      = X(i, 7);
    kT_neg  = X(i, 8);
    kL      = X(i, 9);
    kL_slow = X(i, 10);
    
    % Uniform Model
    R_u  = R_uni_func(kp, kbi, kL, kL_neg, kc, kc_neg, ki);
    RL_u = RL_uni_func(R_u, kL, kL_neg, ki, kc, kc_neg);
    Rc_u = Rc_uni_func(kc, kc_neg, ki, RL_u);
    
    Y_uniform(i) = RL_u + Rc_u; % store steady state signaling value for uniform
    
    % Localized Model 
    R_l  = R_loc_func(kp, kbi, kT, kL_slow, kT_neg, kL, kL_neg, kc, kc_neg, ki);
    T_l  = T_loc_func(R_l, kT, kT_neg, kL);
    RL_l = RL_loc_func(R_l, kT, kT_neg, kL, kL_slow, kL_neg, ki, kc, kc_neg);
    Rc_l = Rc_loc_func(kc, kc_neg, ki, RL_l);
    
    Y_localized(i) = RL_l + Rc_l; % store steady state signaling value for localized
    
    % Ratio of SS signaling levels
    Y_ratio(i) = Y_localized(i) / Y_uniform(i);
end

% Rank in order the (N=1000) parameter sets
X_rank = tiedrank(X);

% Use partial correlation for parameter sets and steady state signaling levels
[R_u_mat, ~] = partialcorr([X_rank, tiedrank(Y_uniform)]);
[R_l_mat, ~] = partialcorr([X_rank, tiedrank(Y_localized)]);
[R_s_mat, ~] = partialcorr([X_rank, tiedrank(Y_ratio)]);

% store PRCC values from partial correlation
prcc_uniform = R_u_mat(1:num_params, end); 
prcc_localized = R_l_mat(1:num_params, end); 
prcc_ratio = R_s_mat(1:num_params, end); 

% Ignore parameters absent in uniform model: tetraspanin rates and slow ligand binding
% k_T, k_Tneg, and k_Lslow 
localized_indices = [7, 8, 10]; 

% Set to zero, no information for uniform model
prcc_uniform(localized_indices)  = 0; 

% plot PRCC bar chart
data = [prcc_uniform(2:8), prcc_localized(2:8), prcc_ratio(2:8)];
bar_plot = barh(data);
hold on

% plot properties
bar_plot(1).FaceColor = 'b'; 
bar_plot(2).FaceColor = 'r'; 
bar_plot(3).FaceColor = 'g'; 

xlim([-1.1, 1.1]);
xlabel('Partial Rank Correlation Coefficient (PRCC)');
ylabel('Model parameters')
title('Global Sensitivity Analysis: Steady-State Signaling (R_L + R_c)');
yticks(1:num_params);

names_sorted  = param_names;

yticklabels(names_sorted(2:8));
legend({'Uniform Model', 'Localized Model', 'Steady state signaling ratio (Localized/Uniform)'},'Location', 'southeast');
legend box off
grid off 
set(gca, 'FontSize', 30);
set(gcf,'Color','white')
box off
