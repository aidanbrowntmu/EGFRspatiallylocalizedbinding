clc
clear 
close all

% Steady state panel B
%-------------------------------------------%

% define the Localized equations---------------------

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
kLT = logspace(-5,1,30); 
kL_neg = 0.05; 
kbi = 1*(3*10^-4); 
kL = 0.1*kLT;

kc = kc/4;

R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kLT,RL+Rc,'rd','linewidth',2,'markersize',14);
hold on
semilogx(kLT,RL,'r^','linewidth',2,'markersize',14);
semilogx(kLT,Rc,'rs','linewidth',2,'markersize',14);

kc = 0.1;

kc = 0.1*4;

R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kLT,RL+Rc,'bd','linewidth',2,'markersize',14);
hold on
semilogx(kLT,RL,'b^','linewidth',2,'markersize',14);
semilogx(kLT,Rc,'bs','linewidth',2,'markersize',14);

set(gca','fontsize',28)
xlabel('k_{+L} (s^{-1})','fontsize',30);
ylabel('Signaling EGFR (R_L+R_c)','fontsize',30);
%legend('k_{+c} = k_{+c}/4 (R_L+R_c)','k_{+c} = k_{+c}/4 (R_L)','k_{+c} = k_{+c}/4 (R_c)','k_{+c} = k_{+c}*4 (R_L+R_c)','k_{+c} = k_{+c}*4 (R_L)','k_{+c} = k_{+c}*4 (R_c)',Fontsize = 12,location = 'northwest');
set(gcf,'Color','white')
yticks([0 20 40 60 80 100])
xticks([10^-4 10^-2 10^0])
%yticklabel({'0' '20' '40' '60' '80' '100'})
ylim([0 100])
box off 

% Plot colours
% Default uniform params will be blue and default localized will be red
% kT = Orange
% kT_neg = gray
% kc = purple
% kc_neg = green
% kL_neg = cyan
% ki = magenta
% kbi = brown

%legend box off


% additional legend details

% Line 1 (solid)
% Making line
annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.82 0.15 0.05], 'String', 'R_L+R_c', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FitBoxToText','on','FontSize', 30);

% Line 2 (dashed)
% Making line
annotation('line', [0.15 0.2], [0.8 0.8], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.775 0.15 0.05], 'String', 'R_L', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Line 2 (dashed)
% Making line
annotation('line', [0.15 0.2], [0.75 0.75], 'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.725 0.15 0.05], 'String', 'R_c', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);


% for curves

% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+c} = 0.25k_{+c0}', 'Color','r',...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+c} = 4k_{+c0}', 'Color','b',...
           'LineStyle', 'none','FontSize', 30);
