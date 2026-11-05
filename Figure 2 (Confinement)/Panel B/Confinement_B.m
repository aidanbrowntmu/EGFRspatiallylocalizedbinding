close all;
clear all;
clc;

% Confinement figure 1, panel 2
%--------------------------------------------%

% define the Localized model equations
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

kT0 = kT;
kT_neg0 = kT_neg;
kc0 = kc;
kc_neg0 = kc_neg;
ki0 = ki;
kL_neg0 = kL_neg;
kbi0 = kbi;


figure(1)

kc = kc0/4;
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kLT,(T+Rc)./(R+T+RL+Rc),'rd','linewidth',2,'markersize',14);
hold on
semilogx(kLT,(T)./(R+T+RL+Rc),'r^','linewidth',2,'markersize',14);
semilogx(kLT,(Rc)./(R+T+RL+Rc),'rsq','linewidth',2,'markersize',14);


x = kLT;
y1 = (T+Rc)./(R+T+RL+Rc);
y2 = (T)./(R+T+RL+Rc);
y3 = (Rc)./(R+T+RL+Rc);

kc = kc0;
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kLT,(T+Rc)./(R+T+RL+Rc),'bd','linewidth',2,'markersize',14);
hold on
semilogx(kLT,(T)./(R+T+RL+Rc),'b^','linewidth',2,'markersize',14);
semilogx(kLT,(Rc)./(R+T+RL+Rc),'bsq','linewidth',2,'markersize',14);


kc = 2*kc0;
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kLT,(T+Rc)./(R+T+RL+Rc),'d','color',[0 0.5 0],'linewidth',2,'markersize',14);
hold on
semilogx(kLT,(T)./(R+T+RL+Rc),'^','color',[0 0.5 0],'linewidth',2,'markersize',14);
semilogx(kLT,(Rc)./(R+T+RL+Rc),'sq','color',[0 0.5 0],'linewidth',2,'markersize',14);

% plot properties
set(gca','fontsize',28)
xlabel('k_{+L} (s^{-1})','FontSize',30);
ylabel('Confined EGFR fraction','FontSize',30);
set(gcf,'Color','white')
axis tight
ylim([0 1])
xticks([10^-4 10^-2 10^0])


%legend('k_{+c} = k_{+c0}/4 confined (T+R_c) fraction','k_{+c} = k_{+c0}/4 T fraction','k_{+c} = k_{+c0}/4 R_c fraction'...
%    ,'k_{+c} = k_{+c0} confined (T+R_c) fraction','k_{+c} = k_{+c0} T fraction','k_{+c} = k_{+c0} R_c fraction',...
%    'k_{+c} = k_{+c0}*4 confined (T+R_c) fraction','k_{+c} = k_{+c0}*4 T fraction','k_{+c} = k_{+c0}*4 R_c fraction',location ='northwest',Fontsize=12);

kc = kc0;

% Plot colours
% Default uniform params will be blue and default localized will be red
% kT = Orange
% kT_neg = gray
% kc = purple
% kc_neg = dark green [0 0.5 0]
% kL_neg = dark cyan or turqoise [0, 0.8078, 0.8235]
% ki = magenta
% kbi = brown

%legend('boxoff') 
box off



% additional legend details

% Line 1 (solid)
% Making line
annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.82 0.15 0.05], 'String', 'Confined T+R_c fraction', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FitBoxToText','on','FontSize', 30);

% Line 2 (dashed)
% Making line
annotation('line', [0.15 0.2], [0.8 0.8], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.775 0.15 0.05], 'String', 'T fraction', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Line 2 (dashed)
% Making line
annotation('line', [0.15 0.2], [0.75 0.75], 'LineWidth', 3, 'LineStyle', ':', 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.725 0.15 0.05], 'String', 'R_c fraction', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);


% for curves

% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+c} = 0.25k_{+c0}', 'Color','r',...
           'LineStyle', 'none','FontSize', 30);

% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'All k_{0}', 'Color','b',...
           'LineStyle', 'none','FontSize', 30);

% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+c} = 4k_{+c0}', 'Color',[0, 0.5, 0],...
           'LineStyle', 'none','FontSize', 30);
