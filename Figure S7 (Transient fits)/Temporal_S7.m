clc
clear
close all

% Transient signaling fits (Figure S7)

% localized model temporal response

% localized MODEL-------------------------------------------

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
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kc(i) = params.k;
    % figure(2)
    % plot(t_rec,valuelist,'b.')
    % hold on
    % plot(t_rec,fitted_y,'r.')
end

%reset
kc = kc0;

kc_T_L = 1./k_kc;

figure(1)
loglog(kc_list/kc0,kc_T_L,'--','color',[0.5210 0.0860 0.8190],'LineWidth',4);
hold on 

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
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kc_neg(i) = params.k;
end

kc_neg = kc_neg0;

kc_neg_T_L = 1./k_kc_neg;


figure(1)
loglog(kc_neg_list/kc_neg0,kc_neg_T_L,'--',"Color", [0 0.5 0],'LineWidth',4);
hold on 


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
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_ki(i) = params.k;
end

%reset
ki = ki0;

ki_T_L = 1./k_ki;

figure(1)
loglog(ki_list/ki0,ki_T_L,'m--','LineWidth',4);
hold on 

kL_neg_list = logspace(log10(kL_neg0/100),log10(kL_neg0*100),51);
for i=1:length(kL_neg_list)
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
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kL_neg(i) = params.k;
end

%reset
kL_neg = kL_neg0;

kL_neg_T_L = 1./k_kL_neg;

figure(1)
loglog(kL_neg_list/kL_neg0,kL_neg_T_L,'--',"Color", [0, 0.8078, 0.8235],'LineWidth',4);
hold on 


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
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kbi(i) = params.k;
end

%reset
kbi = kbi0;

kbi_T_L = 1./k_kbi;

figure(1)
loglog(kbi_list(1:100)/kbi0,kbi_T_L,'--','color',[0.5 0.3 0],'LineWidth',4);
hold on 


clear peakheight duration

kT_neg_list = logspace(log10(kT_neg0/100),log10(kT_neg0*100),51);
for i=1:length(kT_neg_list)
    i
    kT_neg= kT_neg_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kT_neg(i) = params.k;
end

%reset
kT_neg = kT_neg0;

kT_neg_T_L = 1./k_kT_neg;

figure(1)
loglog(kT_neg_list/kT_neg0,kT_neg_T_L,'--',"Color",[0.5 0.5 0.5]','LineWidth',4);
hold on 

kT_list = logspace(log10(kT0/100),log10(kT0*100),51);
for i=1:length(kT_list)
    i
    kT= kT_list(i);

    kLT = 0;
    kL = 0.1*kLT;
    [R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kLT_new = 0.01;
    kL_new = 0.1*kLT_new;
    [R_new, T_new, RL_new, Rc_new] = Solving_for_Steadystate_t(kp,kLT_new,kT,kT_neg,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    
    kLT = kLT_new;
    kL = kL_new;
    [~, ~, valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,R0,T0,RL0,Rc0,...
        R_new, T_new, RL_new, Rc_new,Delta_t,threshold);

    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kT(i) = params.k;
end

%reset
kT = kT0;

kT_T_L = 1./k_kT;

figure(1)
loglog(kT_list/kT0,kT_T_L,'--',"Color",[1 0.5 0],'LineWidth',4);
hold on 

figure(1)
set(gca','fontsize',28)

% plot properties
xlabel('k/k_0','fontsize',30);
ylabel('Peak decay timescale (s)','fontsize',30);
%legend('k_c','k_{-c}','k_i','k_{-L}','k_{bi}','k_{-T}','k_T','fontsize',12,'Orientation','horizontal',location = 'south');
axis tight
%xlim([10^-2 10^2])
%ylim([10^1 10^4])
xticks([10^-2 10^-1 10^0 10^1 10^2])
set(gcf,'Color','white')
box off


%clear all
% Uniform MODEL---------------------------------------------------------------------


% Temporal response panel B and E
% Uniform model temporal response

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
    [~, ~, valuelist, t_rec] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);


    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kc(i) = params.k;
end

%reset
kc = kc0;

kc_T_U = 1./k_kc;

figure(2)
loglog(kc_list/kc0,kc_T_U,'-',"Color",[0.5210 0.0860 0.8190],'LineWidth',4);
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
    [~, ~, valuelist, t_rec] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);


    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kc_neg(i) = params.k;
end

kc_neg = kc_neg0;

kc_neg_T_U = 1./k_kc_neg;

figure(2)
loglog(kcneg_list/kc_neg0,kc_neg_T_U,'-','color',[0 0.5 0],'LineWidth',4,'markersize',15)
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
    [~, ~, valuelist, t_rec] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);


    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_ki(i) = params.k;
end

%reset
ki = ki0;

ki_T_U = 1./k_ki;

loglog(ki_list/ki0,ki_T_U,'m-','LineWidth',4,'markersize',15)

disp('done ki')

kLneg_list = logspace(log10(kL_neg0/100),log10(kL_neg0*100),51);
for i=1:length(kLneg_list)
    i
    kL_neg= kLneg_list(i);

    kL = 0;
    [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
    
    kL_new = 0.01;
    [R_new, RL_new, Rc_new] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL_new,kL_neg);
    
    kL = kL_new;
    [~, ~, valuelist, t_rec] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);


    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kL_neg(i) = params.k;
end

%reset
kL_neg = kL_neg0;

kL_neg_T_U = 1./k_kL_neg;

loglog(kLneg_list/kL_neg0,kL_neg_T_U,'-',"Color", [0 0.8078 0.8235],'LineWidth',4,'markersize',15)

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
    [~, ~, valuelist, t_rec] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,R0,RL0,Rc0,...
        R_new, RL_new, Rc_new,Delta_t,threshold);


    % NEW: Curve fitting to extract time to decay
    [params, fitted_y] = Transient_Exponential_decay_fit(t_rec, valuelist);
    k_kbi(i) = params.k;
end

%reset
kbi = kbi0;

kbi_T_U = 1./k_kbi;

loglog(kbi_list(1:100)/kbi0,kbi_T_U,'-','color',[0.5 0.3 0],'LineWidth',4,'markersize',15)

disp('done kbi')

% Plotting peak height (Panel C)

%figure(1)
% loglog(kc_list/kc0,peakheight_kc_L./peakheight_kc_U,'o-','color',[0.5210 0.0860 0.8190],'LineWidth',2,'markersize',15)
% hold on
% loglog(kcneg_list/kc_neg0,peakheight_kc_neg_L./peakheight_kc_neg_U,'o-','color',[0 0.5 0],'LineWidth',2,'markersize',15)
% loglog(ki_list/ki0,peakheight_ki_L./peakheight_ki_U,'mo-','LineWidth',2,'markersize',15)
% loglog(kL_neg_list/kL_neg0,peakheight_kL_neg_L./peakheight_kL_neg_U,'o-',"Color", [0 0.8078 0.8235],'LineWidth',2,'markersize',15)
% loglog(kbi_list(1:100)/kbi0,peakheight_kbi_L./peakheight_kbi_U,'o-','color',[0.5 0.3 0],'LineWidth',2,'markersize',15)

% % plot properties
% set(gca','fontsize',28)
% 
% xlabel('k/k_0','fontsize',30);
% ylabel('Peak height ratio','fontsize',30);
% %axis tight
% xlim([10^-2 10^2])
% ylim([10^-1 10^1])
% %yticks([10^-1 10^0/2 10^0])
% set(gcf,'Color','white')
% xticks([10^-2 10^-1 10^0 10^1 10^2])
% box off


% Plotting peak duration (Panel F)

figure(2)
% loglog(kc_list/kc0,peakduration_kc_L./peakduration_kc_U,'o-','color',[0.5210 0.0860 0.8190],'LineWidth',2,'markersize',15)
% hold on
% loglog(kcneg_list/kc_neg0,peakduration_kc_neg_L./peakduration_kc_neg_U,'o-','color',[0 0.5 0],'LineWidth',2,'markersize',15)
% loglog(ki_list/ki0,peakduration_ki_L./peakduration_ki_U,'mo-','LineWidth',2,'markersize',15)
% loglog(kL_neg_list/kL_neg0,peakduration_kL_neg_L./peakduration_kL_neg_U,'o-',"Color", [0 0.8078 0.8235],'LineWidth',2,'markersize',15)
% loglog(kbi_list(1:100)/kbi0,peakduration_kbi_L./peakduration_kbi_U,'o-','color',[0.5 0.3 0],'LineWidth',2,'markersize',15)

% plot properties
set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('Peak decay timescale (s)','fontsize',30);
%axis tight
%xlim([10^-2 10^2])
%ylim([10^0 10^1])
set(gcf,'Color','white')
xticks([10^-2 10^-1 10^0 10^1 10^2])
box off


%ratio plot
figure(3)
loglog(kc_list/kc0,kc_T_L./kc_T_U,'o-','color',[0.5210 0.0860 0.8190],'LineWidth',2,'markersize',15)
hold on
loglog(kcneg_list/kc_neg0,kc_neg_T_L./kc_neg_T_U,'o-','color',[0 0.5 0],'LineWidth',2,'markersize',15)
loglog(ki_list/ki0,ki_T_L./ki_T_U,'mo-','LineWidth',2,'markersize',15)
loglog(kL_neg_list/kL_neg0,kL_neg_T_L./kL_neg_T_U,'o-',"Color", [0 0.8078 0.8235],'LineWidth',2,'markersize',15)
loglog(kbi_list(1:100)/kbi0,kbi_T_L./kbi_T_U,'o-','color',[0.5 0.3 0],'LineWidth',2,'markersize',15)

% plot properties
set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('Timescale ratio','fontsize',30);
%axis tight
%xlim([10^-2 10^2])
%ylim([10^0 10^1])
set(gcf,'Color','white')
xticks([10^-2 10^-1 10^0 10^1 10^2])
ylim([10^0 10^1])
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

% Making line
annotation('line', [0.15 0.2], [0.85 0.85], 'LineWidth', 3, 'Color', 'k');
% Making text
annotation('textbox', [0.205 0.82 0.15 0.05], 'String', 'Localized/Uniform', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FitBoxToText','on','FontSize', 30);
