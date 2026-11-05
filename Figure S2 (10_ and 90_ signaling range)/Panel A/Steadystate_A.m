clc
clear 
close all

% Steady state panel A
%-------------------------------------------%

% define the uniform equations

R_func = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);
RL_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

kp = 1;
kc =  0.1;
kc_neg = 0.05; 
ki = 0.05; 
kL = logspace(-5,1,50); 
kL_neg = 0.05;
kbi = 1*(3*10^-4);

kc0 = kc;
kc_neg0 = kc_neg;
ki0 = ki;
kL_neg0 = kL_neg;
kbi0 = kbi;

%%%%%%%%%%

kc = 0.05;


figure(1)

R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kL,(RL+Rc),'b','linewidth',5);
hold on

clear kp kc kc_neg kL kL_neg kbi ki R RL Rc

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

semilogx(kLT,RL+Rc,'r--','linewidth',5);
hold on

%plot properties

set(gca,'fontsize',28);
xlabel('k_{+L} (s^{-1})','FontSize',30);
ylabel('Signaling EGFR (R_L+R_c)','FontSize',30);
set(gcf,'Color','white')
%axis tight
xticks([10^-4 10^-2 10^0])
yticks([0 20 40 60])
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


% additional for legend

% Line 1 (solid)
% Making line
annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'Uniform model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Line 2 (dashed)
% Making line
annotation('line', [0.15 0.2], [0.8 0.8], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.775 0.15 0.05], 'String', 'Localized model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);




% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'All k_{0}', 'Color','b',...
           'LineStyle', 'none','FontSize', 30);

% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'All k_{0}', 'Color','r',...
           'LineStyle', 'none','FontSize', 30);
