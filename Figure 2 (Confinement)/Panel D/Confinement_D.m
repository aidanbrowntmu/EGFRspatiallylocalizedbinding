close all;
clear all;
clc;

% Confinement figure 1, panel 5
%--------------------------------------------%

% define Uniform equations

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

figure(1)

clear kLT kL;

kL = 0.01;

%

clear kc
kc = kc0*logspace(-2,2,51);
R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kc/kc0,(Rc)./(R+RL+Rc),'-','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on

clear kc
kc = kc0;

%

clear kc_neg
kc_neg = kc_neg0*logspace(-2,2,51);
R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kc_neg/kc_neg0,(Rc)./(R+RL+Rc),'color',[0 0.5 0],'linewidth',5);

clear kc_neg
kc_neg = kc_neg0;

%

clear ki
ki = ki0*logspace(-2,2,51);
R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(ki/ki0,(Rc)./(R+RL+Rc),'m','linewidth',5);
hold on

clear ki
ki = ki0;

%

clear kL_neg
kL_neg = kL_neg0*logspace(-2,2,51);
R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kL_neg/kL_neg0,(Rc)./(R+RL+Rc),'color',[0, 0.8078, 0.8235],'linewidth',5);

clear kL_neg
kL_neg = kL_neg0;

%

clear kbi
kbi = kbi0*logspace(-2,2,51);
R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kbi/kbi0,(Rc)./(R+RL+Rc),'-','color',[0.6 0.4 0.2],'linewidth',5);
hold on

clear kbi
kbi = kbi0;

% plot properties
set(gca','fontsize',28)
xlabel('k/k_0','FontSize',30);
ylabel('Confined EGFR fraction','FontSize',30);
set(gcf,'Color','white')
axis tight
ylim([0 1])
xticks([10^-2 10^-1 10^0 10^1 10^2])

 legend('k_{+c}','k_{-c}',...
    'k_i','k_{-L}','k_{bi}','Fontsize', 30, 'Orientation','horizontal',location = 'north');
% 
 legend('boxoff')
box off
%%%%%%%%%%

% Plot colours
% Default uniform params will be blue and default localized will be red
% kT = Orange
% kT_neg = gray
% kc = purple
% kc_neg = green
% kL_neg = cyan
% ki = magenta
% kbi = brown



% annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+c}', 'Color',[0.5210 0.0860 0.8190],...
%            'LineStyle', 'none','FontSize', 25);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{-c} = k_{+c0}', 'Color',[0 0.5 0],...
%            'LineStyle', 'none','FontSize', 25);
% 
% annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{-L}', 'Color',[0, 0.8078, 0.8235],...
%            'LineStyle', 'none','FontSize', 25);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{i}', 'Color','m',...
%            'LineStyle', 'none','FontSize', 25);
% 
% annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{bi}', 'Color',[0.6 0.4 0.2],...
%            'LineStyle', 'none','FontSize', 25);
% 

% Line 1 (solid)
% Making line
annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'Uniform model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);
