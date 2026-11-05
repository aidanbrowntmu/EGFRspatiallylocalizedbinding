clc
clear 
close all

% Steady state panel H: inset 2
%-------------------------------------------%

% define the localized equations---------------------

R_func = @(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki) kp ./(kbi + kT + kL - ((kT_neg.*kT)./(kT_neg+kLT)) + ((kL_neg .*(((kLT.*kT) ./ (kT_neg+kLT)) +kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))));
T_func = @(R,kT,kT_neg,kLT) (kT.*R) ./ (kT_neg+kLT);
RL_func = @(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg) -R.*(((kLT.*kT) ./ (kT_neg+kLT)) +kL)./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kp = 1;
kT = (1/30); 
kT_neg = 0.12; 
kc =  0.1;
kc_neg = 0.05; 
ki = 0.05; 
kLT = logspace(-5,1,50); 
kL_neg = 0.05; 
kbi = 1*(3*10^-4);
kL = 0.1*kLT;

R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kLT,(R+T+RL+Rc),'b--','linewidth',10);
hold on


%plot properties

set(gca,'fontsize',40);
xlabel('k_{+L} (s^{-1})');
ylabel('Total EGF receptors');
set(gca','fontsize',30)
set(gcf,'Color','white')
axis tight
xticks([10^-4 10^-2 10^0])
box off

%legend box off

% Plot colours
% Default uniform params will be blue and default localized will be red
% kT = Orange
% kT_neg = gray
% kc = purple
% kc_neg = green
% kL_neg = cyan
% ki = magenta
% kbi = brown

