clc
clear 
close all

% EGFR Manuscript
% Transition duration: time to reach 90% of new SS

% Fig6. Panel C
%-----------------------------------------------------%
% localized model temporal response

% default parameters:

kc0 = 0.1;
kc_neg0 = 0.05;
ki0 = 0.05;
kL_neg0 = 0.05;
kbi0 = 3*10^-4;
kT0 = 1/30;
kT_neg0 = 0.12;
kp0 = 1;

kc = kc0;
kc_neg = kc_neg0;
ki = ki0;
kL_neg = kL_neg0;
kbi = kbi0;
kT = kT0;
kT_neg = kT_neg0;
kp = kp0;

Delta_t = 1e-3;
threshold = 0.1;


kc_list = logspace(log10(kc0/100),log10(kc0*100),51);
for i=1:length(kc_list)
    i
    kc= kc_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
        
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [duration(i)] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
kc = kc0;

trans_dur_kcL = duration;


%----------------------------------------------------------------------------------------------------------------------------

% Add in the other remaning model parameters

kc_neg_list = logspace(log10(kc_neg0/100),log10(kc_neg0*100),51);
for i=1:length(kc_neg_list)
    i
    kc_neg= kc_neg_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [duration(i)] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);
end

kc_neg = kc_neg0;

trans_dur_kc_negL = duration;



ki_list = logspace(log10(ki0/100),log10(ki0*100),51);
for i=1:length(ki_list)
    i
    ki= ki_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [duration(i)] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
ki = ki0;

trans_dur_kiL = duration;


kL_neg_list = logspace(log10(kL_neg0/100),log10(kL_neg0*100),101);
for i=1:length(kL_neg_list)-1
    i
    kL_neg= kL_neg_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [duration(i)] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
kL_neg = kL_neg0;

trans_dur_kL_negL = duration;



kbi_list = logspace(log10(kbi0/100),log10(kbi0*100),101);
for i=1:length(kbi_list)-1
    i
    kbi= kbi_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [duration(i)] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);
end

%reset
kbi = kbi0;

trans_dur_kbiL = duration;


clear duration


% UNIFORM MODEL--------------------------------------------------

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

Delta_t = 1e-3;
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

trans_dur_kcU = duration;

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

trans_dur_kc_negU = duration;


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

trans_dur_kiU = duration;


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

trans_dur_kL_negU = duration;

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

trans_dur_kbiU = duration;

% plotting ratio of localized to uniform transition duration
figure(1)
semilogx(kc_list/kc0,trans_dur_kcL./trans_dur_kcU,'-o','color',[0.5210 0.0860 0.8190],'LineWidth',2,'markersize',10);
hold on
semilogx(kc_neg_list/kc_neg0,trans_dur_kc_negL./trans_dur_kc_negU,'o-','color',[0 0.5 0],'LineWidth',2,'markersize',10);
semilogx(ki_list/ki0,trans_dur_kiL./trans_dur_kiU,'mo-','LineWidth',2,'markersize',10);
semilogx(kL_neg_list(1:100)/kL_neg0,trans_dur_kL_negL./trans_dur_kL_negU,'o-','color',[0, 0.8078, 0.8235],'LineWidth',2,'markersize',10);
semilogx(kbi_list(1:100)/kbi0,trans_dur_kbiL./trans_dur_kbiU,'-o','color',[0.5 0.3 0],'LineWidth',2,'markersize',10);






% plot properties
set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('Transition duration ratio','fontsize',30);
%legend('k_{+c}','k_{-c}','k_i','k_{-L}','k_{bi}','fontsize',30,'Orientation','horizontal',location = 'south');
set(gcf,'Color','white')
ylim([0.5 2])
%legend box off
xticks([10^-2 10^-1 10^0 10^1 10^2])
box off


% Making line
annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.82 0.15 0.05], 'String', 'Localized/Uniform', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FitBoxToText','on','FontSize', 30);
