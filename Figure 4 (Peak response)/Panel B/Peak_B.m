clc
clear
close all


% Temporal response panel E
% Uniform model temporal response: Time series

% Uniform MODEL FIRST
%-----------------------------------------------------------------%
tic
% default parameters:

kc0 = 0.1;
kc_neg0 = 0.05;
ki0 = 0.05;
kL_neg0 = 0.05;
kbi0 = 3*10^-4;

%for indexx = 1:1 % should be same length as curves
kc_param_list = logspace(-3,1,50); % x-axis varying param
%kL_neg_list = logspace(-4,0,8);
 curves = [0.05 0.05*10]; % different curves on the plot

for indexy = 1:length(curves)

% Define the Transient Differential Equations for Regular model

% Equation [9]
f = @(R,RL,kp,kbi,kL,kL_neg) kp +(RL *kL_neg) - R*(kL+kbi); % GOOD
% Equation [10]
f1 = @ (kL,R,RL,kc_neg,Rc,kL_neg,kc) R*kL -RL*(kL_neg+kc) + Rc*kc_neg;
% Equation [11]
f2 = @ (kc,RL,kc_neg,Rc,ki) RL*kc - Rc*(kc_neg+ki);

% % Define the (NEW) parameters in the model

kp = 1;
kc =  0.1; 
kc_neg = 0.05; 
ki = 0.05;

kL_holder = [0]; 
kL_neg =  curves(indexy);
kbi =  3*10^-4;

% Time step
Delta_t = 0.01; 
maxtime = 3600;


% End loop
stop = maxtime/Delta_t;

% Define time when to increase kLT rates:
t1 = 0; % Set time for first kLT increase

kL_new = 10^-2; % Amount kLT rate is increased

for j = 1:length(kL_holder) %stop+2
    kL = kL_holder(j);

% Define the Initial Conditions

[R0, RL0, Rc0] = Solving_for_Steadystate(kp,kc,kc_neg,ki,kbi,kL,kL_neg);

R = R0;
RL = RL0;
Rc = Rc0;

% Initial time
t = 0;

% Record keeping
count = 1;

t_rec(count) = t;
R_rec(count) = R;
RL_rec(count) = RL;
Rc_rec(count) = Rc;
RL_peak1 = 0;
RL_peak2 = 0;

% While loop

while(t<maxtime)

    if t<t1
        kL = kL;
    elseif (t>t1) %(t1 <= t) && (t < t2)
        kL = kL_new;
%     elseif (t2 <=t) && (t < t3)

    end

      % set variables

    R_past = R;
    RL_past = RL;
    Rc_past = Rc;

    % Calculate derivatives

    dR = f(R_past,RL_past,kp,kbi,kL,kL_neg);
    dRL = f1(kL,R_past,RL_past,kc_neg,Rc_past,kL_neg,kc);
    dRc = f2(kc,RL_past,kc_neg,Rc_past,ki);

    % now update to new values (new values incremented)

    R = R_past + dR*Delta_t;
    RL = RL_past + dRL*Delta_t;
    Rc = Rc_past + dRc*Delta_t;

%     t = t + Delta_t;
% 
    count = count + 1;
    t_rec(count) = t;
    R_rec(count) = R;
    RL_rec(count) = RL;
    Rc_rec(count) = Rc;
    kL_rec(count) = kL;
    t = t + Delta_t;

end

 valuelist = RL_rec + Rc_rec;
 total_uniform = R_rec + RL_rec + Rc_rec;

 % plot activated vs time

 figure(1)
 plot(t_rec,(RL_rec+Rc_rec)./max(valuelist),'k','Linewidth',4)
%plot(t_rec,(RL_rec+Rc_rec),'Linewidth',4)

%(max(RL_rec)+max(Rc_rec)
 hold on
 xlabel('Time (min)','FontSize',30)
 ylabel('Signaling receptors (R_L+R_c)','FontSize',30)
 %legend('k_{-L} = 0.05 (EGF)','k_{-L} = 0.05*10 (EREG)','k_{-L} = 0.05*20')

 set(gca, "FontSize",30)
 set(gcf,'color','w');
 % xticks([0 1200 2400 3600 4800 6000 7200])
 % xticklabels({'0','20','40','60','80','100','120'})
 %  %ylim([0 1.1])
 % xlim([0 7200])
 title('EGF')

end
end

X = t_rec;
Y = (RL_rec+Rc_rec)./max(valuelist);

clear R_rec; RL_rec; Rc_rec; valuelist;

% Now do for the tetraspanin model

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
kc_varying_list = logspace(-3,1,50);
%list2 = logspace(-4,0,5); %kc_neg_list = logspace(-3,0,10);
%curves = [0.05/2 0.05 0.05*4 0.05*10]; %[1/30/2 1/30 2/30 4/30 10/30];
%kT_list = [0.001 0.01];

for indexy = 1:length(curves)

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
kT_neg = 0.12;
kc =  0.1; 
kc_neg =   0.05; 
ki = 0.05; 
% 

% Time step
Delta_t = 0.01; 
maxtime = 3600;

kLT_holder = [0]; %logspace(-4,2,1001); %logspace(-4,0,7);% kLT holder   %logspace(-3,3,2501);%[0.01 0.1 1 10 100 1000];%linspace(0,1000,502);

% End loop
stop = maxtime/Delta_t;

% Define time when to increase kLT rates:
t1 = 0; % Set time for first kLT increase
t2 = 5000;
t3 = 8000;
kLT_new = 10^-2; % Amount kLT rate is increased

for j = 1:length(kLT_holder) %stop+2
    kLT = kLT_holder(j);


% Additional parameters for NEW model
kL_neg = curves(indexy); 
kL = (0.1.*kLT);
kbi = (3*10^-4); 

% Define the Initial Conditions

[R0, T0, RL0, Rc0] = Solving_for_Steadystate_t(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg);

R = R0;
T = T0;
RL = RL0;
Rc = Rc0;

% Initial time
t = 0;

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
        kL = 0.1*kLT;
    elseif (t>t1) %(t1 <= t) && (t < t2) % Stimulating to new ligand concentration
        kLT = kLT_new;
        kL = 0.1*kLT;
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
 total_localized = R_rec + T_rec + RL_rec + Rc_rec;

 % plot activated vs time

 figure(1)
 plot(t_rec,(RL_rec+Rc_rec)./max(valuelist),'k--','Linewidth',4)
 hold on
 ylabel('(R_L+R_c)/max','FontSize',30)
 set(gca, "FontSize",30)
 set(gcf,'color','w');
 ylim([0 1])
 %xlim([0 3600])
 xticks([0 1200 2400 3600])
 xticklabels({'0','20','40','60'})
 yticks([0 0.2 0.5 0.8 1.0])
 yticklabels({'0','0.2','0.5','0.8','1'})
 %legend box off
 legend('Uniform','Localized')
 %legend('ligand = 10^{-2}','ligand = 10^{-1}','ligand = 10^{0}')
 legend box off
 box off
 title('EGF, ligand binding rate = 10^{-2}')

end
end

x = t_rec;
y = (RL_rec+Rc_rec)./max(valuelist);

%(max(RL_rec)+max(Rc_rec)

%-------------------------------%

  % import the experimental data: Freed paper

 %EGF
x = [0.07643474383744878, 314.9684706681671, 1208.0320176649186, 3586.0506592496654];
y = [-0.002659421814911895, 0.9940402165802514, 0.26982608529282337, 0.08936805249966517];

% EREG
x2 = [-4.553111236715894, 301.0694439003643, 1203.3711826560934, 3592.989014555853];
y2 = [0.0048431269886020385, 0.9955294212412905, 0.8988078456643441, 0.5976881872004511];

figure(1)
semilogy(x,y,'bo','Linewidth',4,'MarkerSize',20)
hold on

plot(x2,y2,'ro','Linewidth',4,'MarkerSize',20)

 % import the experimental data: Curran paper

x3 = [-0.30158668103857517, 33.14437624613941, 64.48928529541533, 123.22831787236181, 182.27899001971483, 302.09937839634074, 600.4088175009634, 1799.6883604295936];
y3 = [-0.000543809852242884, 0.6601188423481975, 0.6891264556861288, 1.0010345650847547, 0.910709074992705, 0.9311748945539433, 0.8099981430883094, 0.3429317982863358];

figure(1)
plot(x3,y3,'b+','Linewidth',4,'MarkerSize',20)

%legend('k_{-L} = 0.05 (EGF)','k_{-L} = 0.05*10 (EREG)','k_{-L} = 0.05*20','Experimental: EGF (Freed)','Experimental: EREG (Freed)','Experimental: EGF (Curran)',location='northeast')

%legend box off
box off

% additional for legend

% Line 1 (solid)

% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'o', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 100);

% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'Experimental data (Freed)', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 60);

% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', '+', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 100);
       
% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'Experimental data (Curran)', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 60);

% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'EGF (k_{-L} = 0.05)','color','b', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 60);
% Making text
annotation('textbox', [0.205 0.825 0.15 0.05], 'String', 'EREG (k_{-L} = 0.05*10)','color','r', ...
           'LineStyle', 'none', 'VerticalAlignment', 'middle','FontSize', 60);
       
       
