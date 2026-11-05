clc
clear
close all

% EGFR Manuscript
% Transition duration: time to reach 90% of new SS

% uniform model

% default parameters:

kc0 = 0.1;
kc_neg0 = 0.05;
ki0 = 0.05;
kL_neg0 = 0.05;
kbi0 = 3*10^-4;
kp0 = 1;

kc = kc0;
kc_neg = kc_neg0;
ki = ki0;
kL_neg = kL_neg0;
kbi = kbi0;
kp = kp0;

Delta_t = 1e-4;
threshold = 0.1;


kc_list = logspace(log10(kc0/100),log10(kc0*100),51);
for i=1:length(kc_list)
    i
    kc= kc_list(i);

    kL = 0;
    [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kL_new = 0.01;
    [R_new, RL_new, Rc_new] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kL = kL_new;
    [duration(i)] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
kc = kc0;

% figure(1)
% loglog(kc_list/kc0,peakheight,'-','color',[0.5210 0.0860 0.8190],'LineWidth',5);
% hold on 

figure(2)
loglog(kc_list/kc0,duration,'-','color',[0.5210 0.0860 0.8190],'LineWidth',5);
hold on

disp('done kc')

%----------------------------------------------------------------------------------%

% Add in the other model parameters in the uniform model

kcneg_list = logspace(log10(kc_neg0/100),log10(kc_neg0*100),51);
for i=1:length(kcneg_list)
    i
    kc_neg= kcneg_list(i);

    kL = 0;
    [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kL_new = 0.01;
    [R_new, RL_new, Rc_new] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kL = kL_new;
    [duration(i)] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);
end

kc_neg = kc_neg0;

% figure(1)
% loglog(kcneg_list/kc_neg0,peakheight,'g-','LineWidth',5);
% hold on 

figure(2)
loglog(kcneg_list/kc_neg0,duration,'g-','LineWidth',5);
hold on

disp('done kc_neg')




ki_list = logspace(log10(ki0/100),log10(ki0*100),51);
for i=1:length(ki_list)
    i
    ki= ki_list(i);

    kL = 0;
    [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kL_new = 0.01;
    [R_new, RL_new, Rc_new] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kL = kL_new;
    [duration(i)] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
ki = ki0;

% figure(1)
% loglog(ki_list/ki0,peakheight,'m-','LineWidth',5);
% hold on 

figure(2)
loglog(ki_list/ki0,duration,'m-','LineWidth',5);
hold on

disp('done ki')

kLneg_list = logspace(log10(kL_neg0/100),log10(kL_neg0*100),101);
for i=1:length(kLneg_list)-1
    i
    kL_neg= kLneg_list(i);

    kL = 0;
    [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kL_new = 0.01;
    [R_new, RL_new, Rc_new] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kL = kL_new;
    [duration(i)] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
kL_neg = kL_neg0;

% figure(1)
% loglog(kLneg_list/kL_neg0,peakheight,'c-','LineWidth',5);
% hold on 

figure(2)
loglog(kLneg_list(1:100)/kL_neg0,duration,'c-','LineWidth',5);
hold on

disp('done kL_neg')

kbi_list = logspace(log10(kbi0/100),log10(kbi0*100),101);
for i=1:length(kbi_list)-1
    i
    kbi= kbi_list(i);

    kL = 0;
    [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kL_new = 0.01;
    [R_new, RL_new, Rc_new] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kL = kL_new;
    [duration(i)] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
kbi = kbi0;

% figure(1)
% loglog(kbi_list(1:100)/kbi0,peakheight,'-','color',[0.5 0.3 0],'LineWidth',5);
% hold on 

figure(2)
loglog(kbi_list(1:100)/kbi0,duration,'-','color',[0.5 0.3 0],'LineWidth',5);
hold on

disp('done kbi')

%additional for legend

figure(1)

% additional for legend (applied to all fig 5)

% Line 1 (solid)
% Making line
annotation('line', [0.15 0.2], [0.225 0.225], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.2 0.15 0.05], 'String', 'Uniform model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 20);

% Line 2 (dashed)
% Making line
annotation('line', [0.15 0.2], [0.175 0.175], 'LineWidth', 3, 'LineStyle', '--', 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.15 0.15 0.05], 'String', 'Localized model', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 20);


% % plot properties
% xlabel('k/k_0','fontsize',16);
% ylabel('Peak height','fontsize',20);
% legend('k_{+c}','k_{-c}','k_i','k_{-L}','k_{bi}','fontsize',14,'Orientation','horizontal',location = 'north');
% axis tight
% xlim([10^-2*3 10^2])
% set(gca','fontsize',20)
% set(gcf,'Color','white')
% legend box off
% box off



figure(2)

% plot properties
xlabel('k/k_0','fontsize',16);
ylabel('Transition duration','fontsize',20);
%legend('k_c','k_{-c}','k_i','k_{-L}','k_{bi}','k_{-T}','k_T','fontsize',12,'Orientation','horizontal',location = 'south');
axis tight
xlim([10^-2 10^2])
ylim([10^0 10^2*2])
legend('k_{+c}','k_{-c}','k_i','k_{-L}','k_{bi}','fontsize',30,'Orientation','horizontal',location = 'south');
legend box off
set(gca','fontsize',20)
set(gcf,'Color','white')
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
