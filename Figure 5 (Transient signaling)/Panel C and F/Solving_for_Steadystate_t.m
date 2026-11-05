function [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg)

% Tetraspanin model
%kLT = 0;
%kL = 0;

% Equation [5] 
R0 = kp / ((kbi+kT+kL - ((kT_neg*kT)/(kT_neg+kLT)) + ((kL_neg*(((kLT*kT)/(kT_neg+kLT)) +kL))/(((kc_neg*kc)/(kc_neg+ki)) - kL_neg -kc))));

% Equation [6]
T0 = (kT*R0)/(kT_neg+kLT);

%Equation [7]
RL0 = (-R0*((kLT*kT)/(kT_neg+kLT) + kL)) / ((kc_neg*kc)/(kc_neg+ki) - kL_neg - kc);

% Equation [8]
Rc0 = (kc*RL0) / (kc_neg+ki);











