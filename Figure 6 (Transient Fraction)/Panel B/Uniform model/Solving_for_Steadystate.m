function [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg)

% Non-tetraspanin model
%kL = 0;

% Equation [12]
R0 = kp / ((kL*kL_neg)/((kc*kc_neg)/(kc_neg+ki)-(kL_neg+kc)) +kL+kbi); % GOOD

% Equation [13]
RL0 = -(kL / (((kc*kc_neg)/(kc_neg+ki)) - (kL_neg+kc))) *R0;

% Equation [14]
Rc0 = (kc / (kc_neg+ki)) * RL0;


