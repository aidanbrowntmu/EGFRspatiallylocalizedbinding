clc
clear 
close all

% Steady state panel G inset 1
%-------------------------------------------%

% define the uniform equations

R_func = @(kp,kbi,kL,kL_neg,kc,kc_neg,ki) kp ./(((kL.*kL_neg) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))+kbi+kL);
RL_func = @(R,kL,kL_neg,ki,kc,kc_neg) (R.*-kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

kp = 1;
kc =  0.1;
kc_neg = 0.05; 
ki = 0.05; 
kL = logspace(-5,1,50); 
kL_neg = 0.05;
kbi = 1*(3*10^-4);

kc0 = kc;
kc_neg0 = kc_neg;
ki0 = ki;
kL_neg0 = kL_neg;
kbi0 = kbi;

%%%%%%%%%%

kc = 0.05;


figure(1)

R = R_func(kp,kbi,kL,kL_neg,kc,kc_neg,ki);
RL = RL_func(R,kL,kL_neg,ki,kc,kc_neg);
Rc = Rc_func(kc,kc_neg,ki,RL);

semilogx(kL,(R+RL+Rc),'b','linewidth',10);
hold on

X3 = kL;
Y3 = (R+RL+Rc);
Localized_data3 = [X3(:), Y3(:)];
% Save as a tab-delimited text file with headers
header = {'k_{L}', '(R+RL+Rc)'};
writetable(array2table(Localized_data3, 'VariableNames', header), 'Total_EGFR_pop_inset_H.txt', 'Delimiter', '\t');



%plot properties

set(gca,'fontsize',40);
xlabel('k_{+L} (s^{-1})');
ylabel('Total EGF receptors');
set(gca','fontsize',30)
set(gcf,'Color','white')
axis tight
xticks([10^-4 10^-2 10^0])
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


