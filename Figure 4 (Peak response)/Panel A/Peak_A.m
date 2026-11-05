clc
clear
close all

% Temporal response panel A
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
kc_varying_list = logspace(-3,1,100);
%list2 = logspace(-4,0,5); %kc_neg_list = logspace(-3,0,10);
%curves = [0.05/2 0.05 0.05*4 0.05*10]; %[1/30/2 1/30 2/30 4/30 10/30];
kT_list = [0.12*0.001 0.12*0.01];

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
Delta_t = 0.001; 
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

 % plot activated vs time

 % ---- Large subplot (top 2/3) ----
 figure(1)
 subplot('Position', [0.12 0.38 0.80 0.57]);   
 plot(t_rec,RL_rec+Rc_rec,'--','Linewidth',4)
 hold on
 ylabel('R_L + R_c')
 xlim([-200 1000])
 xticks([-200 0 200 400 600 800 1000])
 xticklabels({'-200','0','200','400','600','800','1000'})
 %ylim([0 15000])

 set(gca, "FontSize",10)
 set(gcf,'color','w');
 %axis tight
 %xlim([0 3600])
 %xticks([0 1200 2400 3600])

 legend box off
 box off

%end

% % Do the exact same but for the uniform model-----------------------------
% 
% %clear all
% 
% for indexy = 1:1
% 
% % define the uniform equations
% 
% % Equation [9]
% f = @(R,RL,kp,kbi,kL,kL_neg) kp +(RL *kL_neg) - R*(kL+kbi); % GOOD
% % Equation [10]
% f1 = @ (kL,R,RL,kc_neg,Rc,kL_neg,kc) R*kL -RL*(kL_neg+kc) + Rc*kc_neg;
% % Equation [11]
% f2 = @ (kc,RL,kc_neg,Rc,ki) RL*kc - Rc*(kc_neg+ki);
% 
% % % Define the (NEW) parameters in the model
% 
% kp = 1;
% kc =  0.1;
% kc_neg = 0.05; 
% ki = 0.05; 
% kL_neg = 0.05;
% kbi = 1*(3*10^-4);
% 
% kc0 = kc;
% kc_neg0 = kc_neg;
% ki0 = ki;
% kL_neg0 = kL_neg;
% kbi0 = kbi;
% 
% % 
% 
% % Time step
% Delta_t = 0.001; 
% maxtime = 3600;
% 
% kL = [0]; %logspace(-4,2,1001); %logspace(-4,0,7);% kLT holder   %logspace(-3,3,2501);%[0.01 0.1 1 10 100 1000];%linspace(0,1000,502);
% 
% % End loop
% stop = maxtime/Delta_t;
% 
% % Define time when to increase kLT rates:
% t1 = 0; % Set time for first kLT increase
% t2 = 5000;
% t3 = 8000;
% kL_new = 10^-2; % Amount kLT rate is increased
% 
% % Additional parameters for NEW model
% kL_neg = 0.05; 
% kbi = (3*10^-4); 
% 
% % Define the Initial Conditions
% 
% [R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);
% 
% R = R0;
% RL = RL0;
% Rc = Rc0;
% 
% % Initial time
% t = -200;
% 
% % Record keeping
% count = 1;
% 
% t_rec(count) = t;
% R_rec(count) = R;
% RL_rec(count) = RL;
% Rc_rec(count) = Rc;
% kL_rec(count) = kL;
% RL_peak1 = 0;
% RL_peak2 = 0;
% 
% % While loop
% 
%  while(t<=maxtime) 
%     if t<t1
%         kL = kL;
%     elseif (t>t1) %(t1 <= t) && (t < t2) % Stimulating to new ligand concentration
%         kL = kL_new;
% %     elseif (t2 <=t) && (t < t3)  % Stimulating to new ligand concentration
% %         kLT = kLT_new*5;
%     end
% 
%            % set variables
% 
%     R_past = R;
%     RL_past = RL;
%     Rc_past = Rc;
% 
%     % Calculate derivatives
% 
%     dR = f(R,RL,kp,kbi,kL,kL_neg);
%     dRL = f1(kL,R,RL,kc_neg,Rc,kL_neg,kc);
%     dRc = f2(kc,RL,kc_neg,Rc,ki);
% 
%     % now update to new values (new values incremented)
% 
%     R = R_past + dR*Delta_t;
%     RL = RL_past + dRL*Delta_t;
%     Rc = Rc_past + dRc*Delta_t;
% 
%     t = t + Delta_t;
%         %imm_frac2 = (T_rec+Rc_rec)./(R_rec + T_rec+RL_rec+Rc_rec);
% 
%     count = count + 1;
%     t_rec(count) = t;
%     R_rec(count) = R;
%     RL_rec(count) = RL;
%     Rc_rec(count) = Rc;
%     kL_rec(count) = kL;
% 

 %end

 valuelist = RL_rec + Rc_rec;

 % plot activated vs time
% 
%  % ---- Large subplot (top 2/3) ----
%  figure(1)
%  subplot('Position', [0.12 0.38 0.80 0.57]);   
%  plot(t_rec,RL_rec+Rc_rec,'Linewidth',4)
%  hold on
%  ylabel('Signaling EGFR (R_L+R_c')
%  xlabel('Time (min)')
%  xlim([-200 1000])
%  xticks([-200 0 200 400 600 800 1000])
%  xticklabels({'-200','0','200','400','600','800','1000'})
% 
%  set(gca, "FontSize",30)
%  set(gcf,'color','w');
%  axis tight
%  %xlim([0 3600])
%  %xticks([0 1200 2400 3600])
%  %xticklabels({'0','20','40','60'})
% 
%  legend box off
%  box off
% legend off
end

figure(1)
subplot('Position', [0.12 0.10 0.80 0.22]);  
plot(t_rec, kLT_rec,'Linewidth',14);
xlabel('Time (s)')
ylabel('Ligand binding rate (s^{-1})')
ylim([0 max(kLT_rec)*1.5])
xlim([-200 1000])
set(gca, "FontSize",50)
set(gcf,'color','w');
%axis tight
yticks([0 0.01])
yticklabels({'0','0.01'})
xticks([-200 0 200 400 600 800 1000])
xticklabels({'-200','0','200','400','600','800','1000'})
box off
