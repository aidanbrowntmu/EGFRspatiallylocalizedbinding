close all;
clear all;
clc;

% Confinement figure 1, panel H: kL = 10^-1
%--------------------------------------------%

% define the tetraspanin model equations
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

clear kc
kc = kc0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kc/kc0,(T+Rc)./(R+T+RL+Rc),'--','color',[0.5210 0.0860 0.8190],'linewidth',5);

confinedkc = (T+Rc)./(R+T+RL+Rc);

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

confinedkc_neg = (T+Rc)./(R+T+RL+Rc);

clear kc_neg
kc_neg = kc_neg0;

%

clear ki
ki = ki0*logspace(-2,2,51);
R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
T = T_func(R,kT,kT_neg,kLT);
RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(ki/ki0,(T+Rc)./(R+T+RL+Rc),'m--','linewidth',3);

confinedki = (T+Rc)./(R+T+RL+Rc);

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

confinedkL_neg = (T+Rc)./(R+T+RL+Rc);

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

confinedkbi = (T+Rc)./(R+T+RL+Rc);

clear kbi
kbi = kbi0;

xlabel('k/k_0');
ylabel('Confined EGFR fraction');
set(gca','fontsize',20)
set(gcf,'Color','white')
axis tight
ylim([0 1])

%legend('k_T with medium k_{LT}','k_{-T} with mediun k_{LT}','k_c with medium k_{LT}','k_{-c} with medium k_{LT}',...
%    'k_i with medium k_{LT}','k_{-L} with medium k_{LT}','k_{bi} with medium k_{LT}');




%---------------------------%

% Confinement figure 1, panel H
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

kL = 0.1;

%

clear kc
kc = kc0*logspace(-2,2,51);
R2 = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL2 = RL_func(R2,kL,kL_neg,ki,kc,kc_neg);
Rc2 = Rc_func(kc,kc_neg,ki,RL2);

semilogx(kc/kc0,(Rc2)./(R2+RL2+Rc2),'-','color',[0.5210 0.0860 0.8190],'linewidth',3);

confinedkc2 = (Rc2)./(R2+RL2+Rc2);

hold on

clear kc
kc = kc0;

%

clear kc_neg
kc_neg = kc_neg0*logspace(-2,2,51);
R2 = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL2 = RL_func(R2,kL,kL_neg,ki,kc,kc_neg);
Rc2 = Rc_func(kc,kc_neg,ki,RL2);

semilogx(kc_neg/kc_neg0,(Rc2)./(R2+RL2+Rc2),'-','color',[0 0.5 0],'linewidth',5);

confinedkc_neg2 = (Rc2)./(R2+RL2+Rc2);

clear kc_neg
kc_neg = kc_neg0;

%

clear ki
ki = ki0*logspace(-2,2,51);
R2 = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL2 = RL_func(R2,kL,kL_neg,ki,kc,kc_neg);
Rc2 = Rc_func(kc,kc_neg,ki,RL2);

semilogx(ki/ki0,(Rc2)./(R2+RL2+Rc2),'m','linewidth',5);

confinedki2 = (Rc2)./(R2+RL2+Rc2);

hold on

clear ki
ki = ki0;

%

clear kL_neg
kL_neg = kL_neg0*logspace(-2,2,51);
R2 = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL2 = RL_func(R2,kL,kL_neg,ki,kc,kc_neg);
Rc2 = Rc_func(kc,kc_neg,ki,RL2);

semilogx(kL_neg/kL_neg0,(Rc2)./(R2+RL2+Rc2),'-','color',[0, 0.8078, 0.8235],'linewidth',5);

confinedkL_neg2 = (Rc2)./(R2+RL2+Rc2);

clear kL_neg
kL_neg = kL_neg0;

%

clear kbi
kbi = kbi0*logspace(-2,2,51);
R2 = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL2 = RL_func(R2,kL,kL_neg,ki,kc,kc_neg);
Rc2 = Rc_func(kc,kc_neg,ki,RL2);

semilogx(kbi/kbi0,(Rc2)./(R2+RL2+Rc2),'-','color',[0.6 0.4 0.2],'linewidth',5);

confinedkbi2 = (Rc2)./(R2+RL2+Rc2);

hold on

clear kbi
kbi = kbi0;

% plot properties

xlabel('k/k_0');
ylabel('Confined EGFR fraction');
set(gca','fontsize',20)
set(gcf,'Color','white')
axis tight
ylim([0 1])

%legend('k_c with medium k_L','k_{-c} with medium k_L',...
%    'k_i with medium k_L','k_{-L} with medium k_L','k_{bi} with medium k_L',location='northwest');

%%%%%%%%%%



% make ratio plots

figure(99)

loglog(logspace(-2,2,51), confinedkc./confinedkc2,'o-','color',[0.5210 0.0860 0.8190],'linewidth',2,'markersize',15)
hold on
loglog(logspace(-2,2,51), confinedkc_neg./confinedkc_neg2,'go-','color',[0 0.5 0],'linewidth',2,'markersize',15)
loglog(logspace(-2,2,51), confinedki./confinedki2,'mo-','linewidth',2,'markersize',15)
loglog(logspace(-2,2,51), confinedkL_neg./confinedkL_neg2,'o-','color',[0, 0.8078, 0.8235],'linewidth',2,'markersize',15)
loglog(logspace(-2,2,51), confinedkbi./confinedkbi2,'o-','color',[0.6 0.4 0.2],'linewidth',2,'markersize',15)

set(gca','fontsize',28)
xlabel('k/k_0','FontSize',30);
ylabel('Confined EGFR Fraction Ratio','FontSize',30);
%legend('k_{+c}','k_{-c}','k_i','k_{-L}','k_{bi}',Fontsize = 16,location = 'west')
set(gcf,'Color','white')
axis tight
ylim([10^0*0.5 10^2*1.5])
yticks([10^0 10^1 10^2])
xticks([10^-2 10^-1 10^0 10^1 10^2])

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


annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.82 0.15 0.05], 'String', 'Localized/Uniform', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FitBoxToText','on','FontSize', 30);

% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'o', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', '-', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);
