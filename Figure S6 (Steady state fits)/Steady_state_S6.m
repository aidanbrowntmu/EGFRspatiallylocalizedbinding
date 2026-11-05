clc
clear 
close all

% Steady state fits

% define the Uniform model equations 

R_func = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);
RL_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

% Do for Uniform model

% do kc half

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
    kL = logspace(-4,4,101);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kL, valuelist);
    x0(i) = params.kL_0;
end

kL_half_kc = x0;


% do kcneg half

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
    kL = logspace(-4,4,101);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kL, valuelist);
    x0(i) = params.kL_0;
end

kL_half_kcneg = x0;


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
    kL = logspace(-4,4,101);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
     % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kL, valuelist);
    x0(i) = params.kL_0;
end

kL_half_ki = x0;


% % do kbi max and half
% 
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
    kL = logspace(-10,4,101);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kL, valuelist);
    x0(i) = params.kL_0;
end

kL_half_kbi = x0;


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
    kbi = kbi_0; 
    kL = logspace(-4,4,501);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kL, valuelist);
    x0(i) = params.kL_0;
end

kL_half_kLneg = x0;


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
    kLT = logspace(-4,4,101);
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_kc = x0;


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
    kLT = logspace(-4,4,101);
    kL_neg = 0.05;
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_kcneg = x0;


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
    kLT = logspace(-4,4,101); 
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4);
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_ki = x0;


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
    kLT = logspace(-10,4,101); 
    kL_neg = 0.05; 
    kbi = kbi_vals(i);
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_kbi = x0;


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
    kLT = logspace(-4,4,101); 
    kL_neg = kLneg_vals(i); 
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_kLneg = x0;


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
    kLT = logspace(-8,4,101);
    kL_neg = kLneg_0; 
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_kT = x0;


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
    kLT = logspace(-8,4,101);
    kL_neg = kLneg_0; 
    kbi = kbi_0;
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    % NEW: Curve fitting sigmoid
    valuelist = RL+Rc;
    [params, fitted_y] = SteadyState_Sigmoid_fit(kLT, valuelist);
    x0(i) = params.kL_0;
end

kLT_half_kTneg = x0;


% figure 1: Localized model

figure(1)
loglog(kc_vals/kc_0,kLT_half_kc,'--','color',[0.5210 0.0860 0.8190],'linewidth',4,'markersize',15);
hold on
loglog(kcneg_vals/kcneg_0,kLT_half_kcneg,'--','color',[0 0.5 0],'linewidth',4,'markersize',15);
loglog(ki_vals/ki_0,kLT_half_ki,'m--','linewidth',4,'markersize',15);
loglog(kbi_vals/kbi_0,kLT_half_kbi,'--','color',[0.5 0.3 0],'linewidth',4,'markersize',15);
loglog(kLneg_vals/kLneg_0,kLT_half_kLneg,'--','color',[0, 0.8078, 0.8235],'linewidth',4,'markersize',15);
loglog(kT_vals/kT_0,kLT_half_kT,'--','color',[1 0.5 0],'linewidth',4,'markersize',15);
loglog(kTneg_vals/kTneg_0,kLT_half_kTneg,'--','color',[0.5 0.5 0.5],'linewidth',4,'markersize',15);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at half sigmoid height (s^{-1})','fontsize',30);
%legend('Uniform: k_c','Uniform: k_{-c}','Uniform: k_i','Uniform: k_{bi}','Uniform: k_{-L}','fontsize',14);
%legend('k_c','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',20,location='northwest');
set(gcf,'Color','white')
axis tight
%ylim([10^0 10^1])
%xlim([10^-2 10^2])
xticks([10^-2 10^-1 10^0 10^1 10^2])
%legend box off
box off

% figure 2: Uniform

figure(2)
loglog(kc_vals/kc_0,kL_half_kc,'-','color',[0.5210 0.0860 0.8190],'linewidth',4,'markersize',15);
hold on
loglog(kcneg_vals/kcneg_0,kL_half_kcneg,'-','color',[0 0.5 0],'linewidth',4,'markersize',15);
loglog(ki_vals/ki_0,kL_half_ki,'m-','linewidth',4,'markersize',15);
loglog(kbi_vals/kbi_0,kL_half_kbi,'-','color',[0.5 0.3 0],'linewidth',4,'markersize',15);
loglog(kLneg_vals/kLneg_0,kL_half_kLneg,'-','color',[0, 0.8078, 0.8235],'linewidth',4,'markersize',15);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at half sigmoid height (s^{-1})','fontsize',30);
%legend('Uniform: k_c','Uniform: k_{-c}','Uniform: k_i','Uniform: k_{bi}','Uniform: k_{-L}','fontsize',14);
%legend('k_c','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',20,location='northwest');
set(gcf,'Color','white')
axis tight
%ylim([10^0 10^1])
%xlim([10^-2 10^2])
xticks([10^-2 10^-1 10^0 10^1 10^2])
%legend box off
box off

% plot properties
% Ratio plot
figure(3)
loglog(kc_vals/kc_0,kLT_half_kc./kL_half_kc,'o-','color',[0.5210 0.0860 0.8190],'linewidth',4,'markersize',15);
hold on
loglog(kcneg_vals/kcneg_0,kLT_half_kcneg./kL_half_kcneg,'o-','color',[0 0.5 0],'linewidth',4,'markersize',15);
loglog(ki_vals/ki_0,kLT_half_ki./kL_half_ki,'mo-','linewidth',4,'markersize',15);
loglog(kbi_vals/kbi_0,kLT_half_kbi./kL_half_kbi,'o-','color',[0.5 0.3 0],'linewidth',4,'markersize',15);
loglog(kLneg_vals/kLneg_0,kLT_half_kLneg./kL_half_kLneg,'o-','color',[0, 0.8078, 0.8235],'linewidth',4,'markersize',15);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at half sigmoid height ratio','fontsize',30);
%legend('Uniform: k_c','Uniform: k_{-c}','Uniform: k_i','Uniform: k_{bi}','Uniform: k_{-L}','fontsize',14);
%legend('k_c','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',20,location='northwest');
set(gcf,'Color','white')
axis tight
ylim([10^0 10^1])
xlim([10^-2 10^2])
xticks([10^-2 10^-1 10^0 10^1 10^2])
%legend box off
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