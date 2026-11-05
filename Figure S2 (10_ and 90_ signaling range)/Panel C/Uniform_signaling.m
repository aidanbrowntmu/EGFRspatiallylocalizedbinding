clc
clear 
close all

% Steady state panel D
%-------------------------------------------%

% define Uniform  model equations:

R_func = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);
RL_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

% Do for Uniform  model

% do kc half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kc_0 = 0.1;
kc_vals = kc_0*logspace(-2,2,101);
for i=1:length(kc_vals)
    kp = 1;
    kc =  kc_vals(i);
    kc_neg = 0.05; 
    ki = 0.05;
    kL_neg = 0.05;
    kbi = (3*10^-4); 
    kL = logspace(-8,0,1001);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kL_10(i) = kL(index10-1) + fraction10*(kL(index10) - kL(index10-1));
    kL_90(i) = kL(index90-1) + fraction90*(kL(index90) - kL(index90-1));
end

kL_half_kc = kL_half;
kL_10_kc = kL_10;
kL_90_kc = kL_90;

% do kcneg half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kcneg_0 = 0.05;
kcneg_vals = kcneg_0*logspace(-2,2,101);
for i=1:length(kcneg_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_vals(i);
    ki = 0.05;  
    kL_neg = 0.05;
    kbi = 1*(3*10^-4); 
    kL = logspace(-8,0,1001);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
       maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kL_10(i) = kL(index10-1) + fraction10*(kL(index10) - kL(index10-1));
    kL_90(i) = kL(index90-1) + fraction90*(kL(index90) - kL(index90-1));
end

kL_half_kcneg = kL_half;
kL_10_kcneg = kL_10;
kL_90_kcneg = kL_90;

% do ki max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

ki_0 = 0.05;
ki_vals = ki_0*logspace(-2,2,101);
for i=1:length(ki_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_vals(i); 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4);
    kL = logspace(-8,0,1001);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
       maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kL_10(i) = kL(index10-1) + fraction10*(kL(index10) - kL(index10-1));
    kL_90(i) = kL(index90-1) + fraction90*(kL(index90) - kL(index90-1));
end

kL_half_ki = kL_half;
kL_10_ki = kL_10;
kL_90_ki = kL_90;
% % do kbi max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kbi_0 = 3*10^-4;
kbi_vals = kbi_0*logspace(-2,2,101); 
for i=1:length(kbi_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;
    kL_neg = 0.05; 
    kbi = kbi_vals(i); 
    kL = logspace(-10,0,1001);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kL_10(i) = kL(index10-1) + fraction10*(kL(index10) - kL(index10-1));
    kL_90(i) = kL(index90-1) + fraction90*(kL(index90) - kL(index90-1));
end

kL_half_kbi = kL_half;
kL_10_kbi = kL_10;
kL_90_kbi = kL_90;

% do kLneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kLneg_0 = 0.05;
kLneg_vals = kLneg_0*logspace(-2,2,101);
for i=1:length(kLneg_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0; 
    kL_neg = kLneg_vals(i);
    kbi = kbi_0; 
    kL = logspace(-8,0,1001);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kL_10(i) = kL(index10-1) + fraction10*(kL(index10) - kL(index10-1));
    kL_90(i) = kL(index90-1) + fraction90*(kL(index90) - kL(index90-1));
end

kL_half_kLneg = kL_half;
kL_10_kLneg = kL_10;
kL_90_kLneg = kL_90;
% plot properties

figure(1)
semilogx(kc_vals/kc_0,(kL_90_kc./kL_10_kc),'-','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
semilogx(kcneg_vals/kcneg_0,(kL_90_kcneg./kL_10_kcneg),'-','color',[0 0.5 0],'linewidth',5);
semilogx(ki_vals/ki_0,(kL_90_ki./kL_10_ki),'m-','linewidth',5);
semilogx(kbi_vals/kbi_0,(kL_90_kbi./kL_10_kbi),'-','color',[0.5 0.3 0.1],'linewidth',5);
semilogx(kLneg_vals/kLneg_0,(kL_90_kLneg./kL_10_kLneg),'-','color',[0, 0.8078, 0.8235],'linewidth',5);

xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at 90% of max / k_{+L} at 10% of max','fontsize',30);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',14,location = 'southeast');
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{-L}','fontsize',14);
set(gca','fontsize',28)
set(gcf,'Color','white')
%axis tight
xlim([10^-2 10^2])
ylim([50 100])
%xticks([10^-4 10^-2 10^0 10^2 10^4])
%yticks([0 20 40 60 80 100])
title('90 %/10%')
%axis tight
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

%%

% For individual legend curves 

% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{+c}', 'Color',[0.5210 0.0860 0.8190]',...
           'LineStyle', 'none','FontSize', 30);

% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-c}', 'Color',[0 0.5 0],...
           'LineStyle', 'none','FontSize', 30);
       
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{i}', 'Color','m',...
           'LineStyle', 'none','FontSize', 30);
       
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{bi}', 'Color',[0.5 0.3 0],...
           'LineStyle', 'none','FontSize', 30);
       
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-L}', 'Color',[0, 0.8078, 0.8235],...
           'LineStyle', 'none','FontSize', 30);
     
       
       
%% plot 10% of ligand

figure(2)
semilogx(kc_vals/kc_0,(kL_10_kc),'-','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
semilogx(kcneg_vals/kcneg_0,(kL_10_kcneg),'-','color',[0 0.5 0],'linewidth',5);
semilogx(ki_vals/ki_0,(kL_10_ki),'m-','linewidth',5);
semilogx(kbi_vals/kbi_0,(kL_10_kbi),'-','color',[0.5 0.3 0.1],'linewidth',5);
semilogx(kLneg_vals/kLneg_0,(kL_10_kLneg),'-','color',[0, 0.8078, 0.8235],'linewidth',5);

xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at 90% of max / k_{+L} at 10% of max','fontsize',30);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',14,location = 'southeast');
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{-L}','fontsize',14);

set(gca','fontsize',28)
set(gcf,'Color','white')
axis tight
%xlim([10^-2 10^2])
%ylim([10^-5*0.5 10^-1])
%xticks([10^-2 10^-1 10^0 10^1 10^2])
%yticks([10^-5 10^-4 10^-3 10^-2 10^-1])
title('10%')

box off


% plot 90% 

figure(3)
semilogx(kc_vals/kc_0,(kL_90_kc),'-','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
semilogx(kcneg_vals/kcneg_0,(kL_90_kcneg),'-','color',[0 0.5 0],'linewidth',5);
semilogx(ki_vals/ki_0,(kL_90_ki),'m-','linewidth',5);
semilogx(kbi_vals/kbi_0,(kL_90_kbi),'-','color',[0.5 0.3 0.1],'linewidth',5);
semilogx(kLneg_vals/kLneg_0,(kL_90_kLneg),'-','color',[0, 0.8078, 0.8235],'linewidth',5);

xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at 90% of max / k_{+L} at 10% of max','fontsize',30);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',14,location = 'southeast');
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{-L}','fontsize',14);

set(gca','fontsize',28)
set(gcf,'Color','white')
axis tight
%xlim([10^-2 10^2])
%ylim([10^-5*0.5 10^-1])
%xticks([10^-2 10^-1 10^0 10^1 10^2])
%yticks([10^-5 10^-4 10^-3 10^-2 10^-1])
title('90%')

box off
