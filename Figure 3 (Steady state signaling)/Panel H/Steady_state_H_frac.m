clc
clear 
close all

% Steady state panel C
%-------------------------------------------%

% define Uniform  model equations:

R_func = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);
RL_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

% Do for Uniform  model

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
    kL = logspace(-4,4,51);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

        % find index for max
    ind = find(RL+Rc==max(RL+Rc));
    maxsignal(i) = max(RL+Rc)./(R(ind)+RL(ind)+Rc(ind));

    totalsignal = (RL+Rc)./(R+RL+Rc); % for this set of parameters 
    
    %maxsignal(i) = max(RL+Rc);
    index = find(totalsignal>(0.5*maxsignal(i)),1);
    previoussignal = totalsignal(index-1); %RL(index-1) + Rc(index-1);
    nextsignal = totalsignal(index); %RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
end

kL_half_kc = kL_half;


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
    kL = logspace(-4,4,51);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
  % find index for max
    ind = find(RL+Rc==max(RL+Rc));
    maxsignal(i) = max(RL+Rc)./(R(ind)+RL(ind)+Rc(ind));

    totalsignal = (RL+Rc)./(R+RL+Rc); % for this set of parameters 
    
    %maxsignal(i) = max(RL+Rc);
    index = find(totalsignal>(0.5*maxsignal(i)),1);
    previoussignal = totalsignal(index-1); %RL(index-1) + Rc(index-1);
    nextsignal = totalsignal(index); %RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
end

kL_half_kcneg = kL_half;


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
    kL = logspace(-4,4,51);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
  % find index for max
    ind = find(RL+Rc==max(RL+Rc));
    maxsignal(i) = max(RL+Rc)./(R(ind)+RL(ind)+Rc(ind));

    totalsignal = (RL+Rc)./(R+RL+Rc); % for this set of parameters 
    
    %maxsignal(i) = max(RL+Rc);
    index = find(totalsignal>(0.5*maxsignal(i)),1);
    previoussignal = totalsignal(index-1); %RL(index-1) + Rc(index-1);
    nextsignal = totalsignal(index); %RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
end

kL_half_ki = kL_half;

% % do kbi max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL;

kbi_0 = 3*10^-4;
kbi_vals = kbi_0*logspace(-4,4,51); 
for i=1:length(kbi_vals)
    kp = 1;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;
    kL_neg = 0.05; 
    kbi = kbi_vals(i); 
    kL_kbi = logspace(-10,4,51); % Change this to be smaller ligand concentration 
    
    R = R_func(kp,kbi,kL_kbi,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL_kbi,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

  % find index for max
    ind = find(RL+Rc==max(RL+Rc));
    maxsignal(i) = max(RL+Rc)./(R(ind)+RL(ind)+Rc(ind));

    totalsignal = (RL+Rc)./(R+RL+Rc); % for this set of parameters 
    
    %maxsignal(i) = max(RL+Rc);
    index = find(totalsignal>(0.5*maxsignal(i)),1);
    previoussignal = totalsignal(index-1); %RL(index-1) + Rc(index-1);
    nextsignal = totalsignal(index); %RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL_kbi(index-1) + fraction*(kL_kbi(index) - kL_kbi(index-1));
end

kL_half_kbi = kL_half;


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
    kL = logspace(-4,4,51);
    
    R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
    RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
  % find index for max
    ind = find(RL+Rc==max(RL+Rc));
    maxsignal(i) = max(RL+Rc)./(R(ind)+RL(ind)+Rc(ind));

    totalsignal = (RL+Rc)./(R+RL+Rc); % for this set of parameters 
    
    %maxsignal(i) = max(RL+Rc);
    index = find(totalsignal>(0.5*maxsignal(i)),1);
    previoussignal = totalsignal(index-1); %RL(index-1) + Rc(index-1);
    nextsignal = totalsignal(index); %RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kL_half(i) = kL(index-1) + fraction*(kL(index) - kL(index-1));
end

kL_half_kLneg = kL_half;

% plot properties

figure(1)
loglog(kc_vals/kc_0,kL_half_kc,'-','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
loglog(kcneg_vals/kcneg_0,kL_half_kcneg,'-','color',[0 0.5 0],'linewidth',5);
loglog(ki_vals/ki_0,kL_half_ki,'m-','linewidth',5);
loglog(kbi_vals/kbi_0,kL_half_kbi,'-','color',[0.5 0.3 0.1],'linewidth',5);
loglog(kLneg_vals/kLneg_0,kL_half_kLneg,'-','color',[0, 0.8078, 0.8235],'linewidth',5);

xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at half of max fraction (s^{-1})','fontsize',30);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','fontsize',14,location = 'southeast');
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{-L}','fontsize',14);

set(gca','fontsize',28)
set(gcf,'Color','white')
axis tight
xlim([10^-2 10^2])
ylim([10^-2*0.5 10^2])
xticks([10^-2 10^-1 10^0 10^1 10^2])
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

% % For individual legend curves 
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{+c}', 'Color',[0.5210 0.0860 0.8190]',...
%            'LineStyle', 'none','FontSize', 30);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-c}', 'Color','g',...
%            'LineStyle', 'none','FontSize', 30);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{i}', 'Color','m',...
%            'LineStyle', 'none','FontSize', 30);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{bi}', 'Color',[0.5 0.3 0],...
%            'LineStyle', 'none','FontSize', 30);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-L}', 'Color','c',...
%            'LineStyle', 'none','FontSize', 30);
% 