clc
clear 
close all

% Steady state panel C
%-------------------------------------------%

% define  model equations:

R_func = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);
RL_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

% Do for uniform model first 

% do kc max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kc_0 = 0.1;
kc_vals = kc_0*logspace(-2,2,51);
for i=1:length(kc_vals)
    kp = 1;
    kc =  kc_vals(i);
    kc_neg = 0.05; 
    ki = 0.05;  
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4); 
    kL = logspace(-5,5,501);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kc = maxsignal;


% do kcneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kcneg_0 = 0.05;
kcneg_vals = kcneg_0*logspace(-2,2,51);
for i=1:length(kcneg_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_vals(i);
    ki = 0.05; 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4);
    kL = logspace(-5,5,501);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    
end

maxsignal_kcneg = maxsignal;


% do ki max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

ki_0 = 0.05;
ki_vals = ki_0*logspace(-2,2,51);
for i=1:length(ki_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_vals(i); 
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4);
    kL = logspace(-5,5,501);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_ki = maxsignal;

% do kbi max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kbi_0 = 3*10^-4;
kbi_vals = kbi_0*logspace(-2,2,51);
for i=1:length(kbi_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0; 
    kL_neg = 0.05; 
    kbi = kbi_vals(i);
    kL =logspace(-5,5,501);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kbi = maxsignal;

% do kLneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kLneg_0 = 0.05;
kLneg_vals = kLneg_0*logspace(-2,2,51);
for i=1:length(kLneg_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;  
    kL_neg = kLneg_vals(i); 
    kbi = kbi_0; %1/s
    kL = logspace(-5,5,501);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kLneg = maxsignal;


figure(1)
loglog(kc_vals/kc_0,maxsignal_kc,'-','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
loglog(kcneg_vals/kcneg_0,maxsignal_kcneg,'-','color',[0 0.5 0],'linewidth',5);
loglog(ki_vals/ki_0,maxsignal_ki,'m','linewidth',5);
loglog(kbi_vals/kbi_0,maxsignal_kbi,'-','color',[0.5 0.3 0],'linewidth',5);
loglog(kLneg_vals/kLneg_0,maxsignal_kLneg,'-','color',[0, 0.8078, 0.8235],'linewidth',5);
set(gca','fontsize',14)
xlabel('k/k0','fontsize',16);
ylabel('max(R_L + R_c) at high k_{+L}','fontsize',16);


%----------------------%

clear R RL Rc kp kc kc_neg kL kL_neg ki kbi

% Do for localized model

% define localized model equations
% 
R_func = @(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki) kp ./(kbi + kT + kL - ((kT_neg.*kT)./(kT_neg+kLT)) + ((kL_neg .*(((kLT.*kT) ./ (kT_neg+kLT)) +kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))));
T_func = @(R,kT,kT_neg,kLT) (kT.*R) ./ (kT_neg+kLT);
RL_func = @(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg) -R.*(((kLT.*kT) ./ (kT_neg+kLT)) +kL)./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);



% do kc max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kc_0 = 0.1;
kc_vals = kc_0*logspace(-2,2,51);
for i=1:length(kc_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12;
    kc =  kc_vals(i); 
    kc_neg = 0.05; 
    ki = 0.05;  
    kLT = logspace(-5,5,501); 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4);
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kc = maxsignal;


% do kcneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kcneg_0 = 0.05;
kcneg_vals = kcneg_0*logspace(-2,2,51);
for i=1:length(kcneg_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_vals(i);
    ki = 0.05; 
    kLT = logspace(-5,5,501); 
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kcneg = maxsignal;


% do ki max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

ki_0 = 0.05;
ki_vals = ki_0*logspace(-2,2,51);
for i=1:length(ki_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12;
    kc =  kc_0; 
    kc_neg = kcneg_0;
    ki = ki_vals(i); 
    kLT = logspace(-5,5,501);
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_ki = maxsignal;

% do kbi max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kbi_0 = 3*10^-4;
kbi_vals = kbi_0*logspace(-2,2,51);
for i=1:length(kbi_vals)
    kp = 1;
    kT = (1/30);
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0; 
    kLT = logspace(-5,5,501); 
    kL_neg = 0.05;
    kbi = kbi_vals(i);
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kbi = maxsignal;

% do kLneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kLneg_0 = 0.05;
kLneg_vals = kLneg_0*logspace(-2,2,51);
for i=1:length(kLneg_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0; 
    kLT = logspace(-5,5,501); 
    kL_neg = kLneg_vals(i); 
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kLneg = maxsignal;

% do kT max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kT_0 = 1/30;
kT_vals = kT_0*logspace(-2,2,51);
for i=1:length(kT_vals)
    kp = 1;
    kT = kT_vals(i);
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0; 
    kLT = logspace(-5,5,501); 
    kL_neg = kLneg_0;
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kT = maxsignal;

% do kTneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kTneg_0 = 0.12;
kTneg_vals = kTneg_0*logspace(-2,2,51);
for i=1:length(kTneg_vals)
    kp = 1;
    kT = kT_0; 
    kT_neg = kTneg_vals(i); 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;
    kLT = logspace(-5,5,501); 
    kL_neg = kLneg_0; 
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
end

maxsignal_kTneg = maxsignal;

figure(1)
set(gca','fontsize',28)
loglog(kc_vals/kc_0,maxsignal_kc,'+','color',[0.5210 0.0860 0.8190],'linewidth',2,'markersize',15);
hold on
loglog(kcneg_vals/kcneg_0,maxsignal_kcneg,'+','color',[0 0.5 0],'linewidth',2,'markersize',15);
loglog(ki_vals/ki_0,maxsignal_ki,'m+','linewidth',2,'markersize',15);
loglog(kbi_vals/kbi_0,maxsignal_kbi,'+','color',[0.5 0.3 0],'Linewidth',2,'markersize',15);
loglog(kLneg_vals/kLneg_0,maxsignal_kLneg,'+','color',[0, 0.8078, 0.8235],'linewidth',2,'markersize',15);
loglog(kT_vals/kT_0,maxsignal_kT,'+','color',[1 0.5 0],'linewidth',2,'markersize',15);
loglog(kTneg_vals/kTneg_0,maxsignal_kTneg,'+',"Color",[0.5 0.5 0.5],'linewidth',2,'markersize',15);
xlabel('k/k_0','fontsize',30);
ylabel('max(R_L + R_c) at high k_{+L}','fontsize',30);
set(gcf,'Color','white')
yticks([10^1 10^2 10^3])
%axis tight
box off
%legend('uniform:k_{+c}','uniform:k_{-c}','uniform:k_i','uniform:k_{bi}','uniform:k_{-L}','localized: k_{+c}','localized: k_{-c}','localized: k_i','localized: k_{bi}','localized: k_{-L}','localized: k_{+T}','localized: k_{-T}','fontsize',10,Fontsize=8,location = 'west');

%legend box off
%%

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
annotation('line', [0.15 0.2], [0.25 0.25], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.2 0.15 0.05], 'String', 'Uniform model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Line 2 (dashed)
% Making text
annotation('textbox', [0.15 0.18 0.15 0.05], 'String', '+', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Making text
annotation('textbox', [0.205 0.18 0.15 0.05], 'String', 'Localized model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% Making text
annotation('textbox', [0.205 0.18 0.15 0.05], 'String', 'overlap', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 30);

% for curves

% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+c}', 'Color',[0.5210 0.0860 0.8190],...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{-c}', 'color',[0 0.5 0],...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{i}', 'Color','m',...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{bi},', 'Color',[0.5 0.3 0],...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{-L},', 'color',[0, 0.8078, 0.8235],...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{+T},', 'Color',[1 0.5, 0],...
           'LineStyle', 'none','FontSize', 30);


% Make additional text for legend
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'k_{-T}', 'Color',[0.5 0.5 0.5],...
           'LineStyle', 'none','FontSize', 30);
