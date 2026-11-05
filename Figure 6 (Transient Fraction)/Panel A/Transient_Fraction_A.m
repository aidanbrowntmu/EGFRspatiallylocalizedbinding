clc
clear
close all

%Transient Signaling fraction: Fig. 6, Panel A
%-----------------------------------------------------------------%
% Localized MODEL-------------------------------------------

% default parameters:

kc0 = 0.1;
kc_neg0 = 0.05;
ki0 = 0.05;
kL_neg0 = 0.05;
kbi0 = 3*10^-4;
kT0 = 1/30;
kT_neg0 = 0.12;

%for indexx2 = 1:4 % should be same length as curves (changed this to differ from Uniform index, can change back but might cause some errors be careful)
kc_varying_list = logspace(-3,1,10);
%list2 = logspace(-4,0,5); %kc_neg_list = logspace(-3,0,10);
%curves = [0.05/2 0.05 0.05*4 0.05*10]; %[1/30/2 1/30 2/30 4/30 10/30];
kT_list = [0.12/10 0.12];

for indexy = 1:length(kT_list)

% Define the Transient Differential Equations for Localized model

% Equation [1]: dR/dt
f = @(R,RL,T,Rc,kp,kbi,kT,kT_neg,kL,kL_neg) kp - kbi * R - kT * R + kT_neg * T - kL * R + kL_neg*RL;

% Equation [2]: dT/dt
f1 = @(R,T,kT,kT_neg,kLT) kT * R - kT_neg *T - kLT *T;

% Equation [3]: dRL/dt
f2 = @ (kLT,T,kL,R,kc_neg,Rc,kL_neg,RL,kc)  kLT * T + kL * R + kc_neg *Rc - kL_neg*RL -kc*RL;

% Equation [4]: dRc/dt
f3 = @ (kc,RL,kc_neg,Rc,ki) kc * RL - kc_neg *Rc - ki *Rc;

% % Define the (NEW) parameters in the model

kp = 1;
kT = (1/30);
kT_neg =  kT_list(indexy);
kc =  0.1; 
kc_neg =   0.05; 
ki = 0.05; 
% 

% Time step
Delta_t = 5; 
maxtime = 1000;

kLT = [0]; %logspace(-4,2,1001); %logspace(-4,0,7);% kLT holder   %logspace(-3,3,2501);%[0.01 0.1 1 10 100 1000];%linspace(0,1000,502);

% End loop
stop = maxtime/Delta_t;

% Define time when to increase kLT rates:
t1 = 0; % Set time for first kLT increase
t2 = 5000;
t3 = 8000;
kLT_new = 10^-2; % Amount kLT rate is increased



% Additional parameters for NEW model
kL_neg = 0.05; 
kL = (0.1.*kLT);
kbi = (3*10^-4); 

% Define the Initial Conditions

[R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);

R = R0;
T = T0;
RL = RL0;
Rc = Rc0;

% Initial time
t = -200;

% Record keeping
count = 1;

t_rec(count) = t;
R_rec(count) = R;
T_rec(count) = T;
RL_rec(count) = RL;
Rc_rec(count) = Rc;
kLT_rec(count) = kLT;
RL_peak1 = 0;
RL_peak2 = 0;

% While loop

 while(t<=maxtime) 
    if t<t1
        kLT = kLT;
    elseif (t>t1) %(t1 <= t) && (t < t2) % Stimulating to new ligand concentration
        kLT = kLT_new;
%     elseif (t2 <=t) && (t < t3)  % Stimulating to new ligand concentration
%         kLT = kLT_new*5;
    end

           % set variables

    R_past = R;
    T_past = T;
    RL_past = RL;
    Rc_past = Rc;

    % Calculate derivatives

    dR = f(R_past,RL_past,T_past,Rc_past,kp,kbi,kT,kT_neg,kL,kL_neg);
    dT = f1(R_past,T_past,kT,kT_neg,kLT);
    dRL = f2(kLT,T_past,kL,R_past,kc_neg,Rc_past,kL_neg,RL_past,kc);
    dRc = f3(kc,RL_past,kc_neg,Rc_past,ki);

    % now update to new values (new values incremented)
    
    R = R_past + dR*Delta_t;
    T = T_past + dT*Delta_t;
    RL = RL_past + dRL*Delta_t;
    Rc = Rc_past + dRc*Delta_t;
        
    t = t + Delta_t;
        %imm_frac2 = (T_rec+Rc_rec)./(R_rec + T_rec+RL_rec+Rc_rec);
    
    count = count + 1;
    t_rec(count) = t;
    R_rec(count) = R;
    T_rec(count) = T;
    RL_rec(count) = RL;
    Rc_rec(count) = Rc;
    kLT_rec(count) = kLT;


 end

 valuelist = RL_rec + Rc_rec;

 % if indexy == 1
 %    writetable(table(t_rec(:), valuelist(:), 'VariableNames', {'t_rec', 'Signaling EGFR fraction'}), 'kT_neg0.txt', 'Delimiter', '\t');
 % else
 %    writetable(table(t_rec(:), valuelist(:), 'VariableNames', {'t_rec', 'Signaling EGFR fraction'}), '0.1kT_neg0.txt', 'Delimiter', '\t');
 % end
 % 
 % writetable(table(t_rec(:), kLT_rec(:), 'VariableNames', {'t_rec', 'Ligand concentration'}), 'Ligand_conc.txt', 'Delimiter', '\t');


 % plot activated vs time

 % ---- Large subplot (top 2/3) ----
subplot('Position', [0.12 0.38 0.80 0.57]);   
 plot(t_rec,(RL_rec+Rc_rec)./(R_rec+T_rec+RL_rec+Rc_rec),'Linewidth',14)
 hold on
% ylabel('R_L+R_c')
 legend('k_{-T} = 0.1k_{-T0}','k_{-T} = k_{-T0}',location = 'northwest')
 xlim([-100 1500])

 set(gcf,'color','w');
 axis tight
 set(gca, "FontSize",60)
 %xlim([0 0.2])
 %xticks([0 1200 2400 3600])
 ylim([0 0.5])
 legend box off
 box off

end


figure(1)
subplot('Position', [0.12 0.10 0.80 0.22]);  
plot(t_rec, kLT_rec,'Linewidth',14);
%xlabel('Time (s)')
%ylabel({'Ligand binding';'rate (s^{-1})'})
ylim([0 max(kLT_rec)*2])
xlim([-200 1000])
set(gca, "FontSize",60)
set(gcf,'color','w');
%axis tight
yticks([0 0.01])
yticklabels({'0','0.01'})
box off



% plot fraction of EGFR in system

figure(2)
plot(t_rec,R_rec./(R_rec+T_rec+RL_rec+Rc_rec),'linewidth',24)
hold on
plot(t_rec,T_rec./(R_rec+T_rec+RL_rec+Rc_rec),'linewidth',24)
plot(t_rec,RL_rec./(R_rec+T_rec+RL_rec+Rc_rec),'linewidth',24)
plot(t_rec,Rc_rec./(R_rec+T_rec+RL_rec+Rc_rec),'linewidth',24)

xlabel('Time (s)')
ylabel('Concentration')
legend('R/Total EGFR','T/Total EGFR','R_L/Total EGFR','R_c/Total EGFR')
set(gca, "FontSize",80)
set(gcf,'color','w');
legend box off
box off
axis tight
ylim([0 1])


% Make additional text for legend for figure 2
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'R/Total EGFR', 'color',[0 0, 1],...
           'LineStyle', 'none','FontSize', 50);

% Make additional text for legend for figure 2
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'T/Total EGFR,', 'Color',[1 0, 0],...
           'LineStyle', 'none','FontSize', 50);

% Make additional text for legend for figure 2
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'R_L/Total EGFR,', 'color',[0.98, 0.78, 0.05],...
           'LineStyle', 'none','FontSize', 50);

% Make additional text for legend for figure 2
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'R/Total EGFR,', 'color',[0.5210 0.0860 0.8190],...
           'LineStyle', 'none','FontSize', 50);

% Make additional text for legend for figure 2
annotation('textbox', [0.205 0.5 0.15 0.05], 'String', 'overlap', 'color',[0 0 0],...
           'LineStyle', 'none','FontSize', 60);
