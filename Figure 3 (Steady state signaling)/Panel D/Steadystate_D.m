clc
clear 
close all

% Steady state panel E

% define the Localized model equations

R_func = @(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki) kp ./(kbi + kT + kL - ((kT_neg.*kT)./(kT_neg+kLT)) + ((kL_neg .*(((kLT.*kT) ./ (kT_neg+kLT)) +kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))));
T_func = @(R,kT,kT_neg,kLT) (kT.*R) ./ (kT_neg+kLT);
RL_func = @(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg) -R.*(((kLT.*kT) ./ (kT_neg+kLT)) +kL)./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);


%--------------------------%

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
    kLT = logspace(-4,4,51); 
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
end

kLT_half_kc = kLT_half;


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
    kLT = logspace(-4,4,51); 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
end

kLT_half_kcneg = kLT_half;


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
    kLT = logspace(-4,4,51); 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
end

kLT_half_ki = kLT_half;


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
    kLT_kbi = logspace(-10,4,51);
    kL_neg = 0.05; 
    kbi = kbi_vals(i);
    kL = 0.1*kLT_kbi;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT_kbi,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT_kbi);
    RL = RL_func(R,kT,kT_neg, kLT_kbi,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT_kbi(index-1) + fraction*(kLT_kbi(index) - kLT_kbi(index-1));
end

kLT_half_kbi = kLT_half;


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
    kLT = logspace(-4,4,51); 
    kL_neg = kLneg_vals(i); 
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
end

kLT_half_kLneg = kLT_half;


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
    kLT = logspace(-8,4,51);
    kL_neg = kLneg_0;
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
end

kLT_half_kT = kLT_half;


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
    kLT = logspace(-8,4,51); 
    kL_neg = kLneg_0; 
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
end

kLT_half_kTneg = kLT_half;

% plot properties

figure(1)
loglog(kc_vals/kc_0,kLT_half_kc,'--','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
loglog(kcneg_vals/kcneg_0,kLT_half_kcneg,'--','color',[0 0.5 0],'linewidth',5);
loglog(ki_vals/ki_0,kLT_half_ki,'m--','linewidth',4);
loglog(kbi_vals/kbi_0,kLT_half_kbi,'--',"Color", [0.5 0.3 0],'linewidth',5);
loglog(kLneg_vals/kLneg_0,kLT_half_kLneg,'--','color',[0, 0.8078, 0.8235],'linewidth',5);
loglog(kT_vals/kT_0,kLT_half_kT,'--','color',[1 0.5 0],'linewidth',5);
loglog(kTneg_vals/kTneg_0,kLT_half_kTneg,'--',"Color",[0.5 0.5 0.5],'linewidth',5);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at half of max (s^{-1})','fontsize',30);
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{bi}','Uniform : k_{-L}','fontsize',14);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','k_{+T}','k_{-T}','fontsize',12,location='north');
set(gcf,'Color','white')
%axis tight
%ylim([10^0 10^2])
xlim([10^-2 10^2])
ylim([10^-5*0.5 10^-1])
xticks([10^-2 10^-1 10^0 10^1 10^2])
yticks([10^-5 10^-4 10^-3 10^-2 10^-1])
box off
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


% % For individual legend curves 
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{+c}', 'Color',[0.5210 0.0860 0.8190]',...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-c}', 'Color','g',...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{i}', 'Color','m',...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{bi}', 'Color',[0.5 0.3 0],...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-L}', 'Color','c',...
%            'LineStyle', 'none','FontSize', 20);
% 
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{+T}', 'Color',[1 0.5 0],...
           'LineStyle', 'none','FontSize', 30);
       
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-T}', 'Color',[0.5 0.5 0.5],...
           'LineStyle', 'none','FontSize', 30);