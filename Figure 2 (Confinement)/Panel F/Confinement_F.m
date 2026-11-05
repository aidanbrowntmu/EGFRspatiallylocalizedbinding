close all;
clear all;
clc;

% Confinement figure 1, panel F: kLT = 10^-1
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
kLT = logspace(-5,1,50); 
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

clear kLT kL;

kLT = 0.1;
kL = 0.1*kLT;

clear kT
kT = kT0*logspace(-2,2,51);

R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kT/kT0,(T+Rc)./(R+T+RL+Rc),'--','color',[1 0.5 0],'linewidth',5);
hold on

clear kT
kT = kT0;

%

clear kT_neg
kT_neg = kT_neg0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kT_neg/kT_neg0,(T+Rc)./(R+T+RL+Rc),'--','color',[0.5 0.5 0.5],'linewidth',5);

clear kT_neg
kT_neg = kT_neg0;

%

clear kc
kc = kc0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kc/kc0,(T+Rc)./(R+T+RL+Rc),'--','color',[0.5210 0.0860 0.8190],'linewidth',5);

clear kc
kc = kc0;

%

clear kc_neg
kc_neg = kc_neg0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kc_neg/kc_neg0,(T+Rc)./(R+T+RL+Rc),'--','color',[0 0.5 0],'linewidth',5);

clear kc_neg
kc_neg = kc_neg0;

%

clear ki
ki = ki0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(ki/ki0,(T+Rc)./(R+T+RL+Rc),'m--','linewidth',5);

clear ki
ki = ki0;

%

clear kL_neg
kL_neg = kL_neg0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kL_neg/kL_neg0,(T+Rc)./(R+T+RL+Rc),'--','color',[0, 0.8078, 0.8235],'linewidth',5);

clear kL_neg
kL_neg = kL_neg0;

%

clear kbi
kbi = kbi0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kbi/kbi0,(T+Rc)./(R+T+RL+Rc),'--','color',[0.6 0.4 0.2],'linewidth',5);

clear kbi
kbi = kbi0;

set(gca','fontsize',28)
xlabel('k/k_0','FontSize',30);
ylabel('Confined EGFR fraction','FontSize',30);
set(gcf,'Color','white')
axis tight
ylim([0 1])
xticks([10^-2 10^-1 10^0 10^1 10^2])

%legend('k_{+T}','k_{-T}','k_{+c}','k_{-c}',...
%    'k_{i}','k_{-L}','k_{bi}',orientation = 'horizontal',Fontsize = 18,location = 'north');
%
%legend('boxoff') 
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

