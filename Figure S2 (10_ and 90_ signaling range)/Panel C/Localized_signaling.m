clc
clear 
close all

% Steady state panel E (test for ligand concentration ratio: 10 and 90%)

% define the Localized model equations

R_func = @(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki) kp ./(kbi + kT + kL - ((kT_neg.*kT)./(kT_neg+kLT)) + ((kL_neg .*(((kLT.*kT) ./ (kT_neg+kLT)) +kL) ./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc))));
T_func = @(R,kT,kT_neg,kLT) (kT.*R) ./ (kT_neg+kLT);
RL_func = @(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg) -R.*(((kLT.*kT) ./ (kT_neg+kLT)) +kL)./ (((kc.*kc_neg) ./ (kc_neg+ki)) -kL_neg - kc);
Rc_func = @(kc,kc_neg,ki,RL) (kc.*RL) ./(kc_neg+ki);

color = ["bo" , "ro", "ko"];
color2 = ["b-", "r-", "k-"];
color3 = ["b--", "r--", "k--"];
color4 = ["bo", "ro", "ko"];
color5 = ["bx", "rx", "kx"];

%--------------------------%

% do kc max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90;

kc_0 = 0.1;
kc_vals = kc_0*logspace(-2,2,101); %[kc_0 kc_0/10 kc_0/100]; %kc_0*logspace(-2,2,51);
for i=1:length(kc_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12; 
    kc =  kc_vals(i);
    kc_neg = 0.05; 
    ki = 0.05; 
    kLT = logspace(-8,0,1001); 
    kL_neg = 0.05; 
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;

    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

%     figure(99)
%     loglog(kLT,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)

    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);

    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);

    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));

    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);

    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);

    kLT_10(i) = kLT(index10-1) + fraction10*(kLT(index10) - kLT(index10-1));
    kLT_90(i) = kLT(index90-1) + fraction90*(kLT(index90) - kLT(index90-1));

    ratio(i) = kLT_90(i)./kLT_10(i);
    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);

end


% for i = 1:length(kc_vals)
%     figure(99)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
% end

% figure(99)
% legend('k_c = k_{c0}', 'k_c = k_{c0}/10', 'k_c = k_{c0}/100','location','southeast')
% legend box off
% box off
% axis tight
%title(ratio)

kLT_half_kc = kLT_half;
kLT_10_kc = kLT_10;
kLT_90_kc = kLT_90;


% do kcneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90 ratio;

kcneg_0 = 0.05;
kcneg_vals = kcneg_0*logspace(-2,2,101); %[kcneg_0 kcneg_0*10 kcneg_0*100];  %*logspace(-2,2,51);
for i=1:length(kcneg_vals)
    kp = 1;
    kT = (1/30);
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_vals(i);
    ki = 0.05;
    kLT = logspace(-8,0,1001); 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;

    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

%     figure(98)
%     loglog(kLT,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)

    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);

    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);

    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));

    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);

    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);

    kLT_10(i) = kLT(index10-1) + fraction10*(kLT(index10) - kLT(index10-1));
    kLT_90(i) = kLT(index90-1) + fraction90*(kLT(index90) - kLT(index90-1));

    ratio(i) = kLT_90(i)./kLT_10(i);

    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);

end

% for i = 1:length(kcneg_vals)
%     figure(98)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
% end 
% 
% figure(98)
% legend('k_{-c} = k_{-c0}', 'k_{-c} = k_{-c0}*50', 'k_{-c} = k_{-c0}*500','location','southeast')
% legend box off
% box off
% axis tight
% title(ratio)

kLT_half_kcneg = kLT_half;
kLT_10_kcneg = kLT_10;
kLT_90_kcneg = kLT_90;

% do ki max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90 ratio;

ki_0 = 0.05;
ki_vals = ki_0*logspace(-2,2,101); %[ki_0 ki_0/50 ki_0/500]; %ki_0*logspace(-2,2,51);
for i=1:length(ki_vals)
    kp = 1;
    kT = (1/30);
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_vals(i);
    kLT = logspace(-8,0,1001); 
    kL_neg = 0.05;
    kbi = 1*(3*10^-4); 
    kL = 0.1*kLT;

    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
% 
%     figure(97)
%     loglog(kLT,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)

    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);

    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);

    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));

    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);

    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);

    kLT_10(i) = kLT(index10-1) + fraction10*(kLT(index10) - kLT(index10-1));
    kLT_90(i) = kLT(index90-1) + fraction90*(kLT(index90) - kLT(index90-1));
    ratio(i) = kLT_90(i)./kLT_10(i);
    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);

end

% for i = 1:length(ki_vals)
%     figure(97)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
% end 

% figure(97)
% legend('k_i = k_{i0}', 'k_i = k_{i0}/10', 'k_i = k_{i0}/100','location','southeast')
% legend box off
% box off
% axis tight
% title(ratio)

kLT_half_ki = kLT_half;
kLT_10_ki = kLT_10;
kLT_90_ki = kLT_90;

% do kbi max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90 ratio;

kbi_0 = 3*10^-4;
kbi_vals = kbi_0*logspace(-2,2,101); %[kbi_0 kbi_0*10 kbi_0*100]; %kbi_0*logspace(-2,2,51);
for i=1:length(kbi_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12;
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;
    kLT_kbi = logspace(-10,0,1001);
    kL_neg = 0.05; 
    kbi = kbi_vals(i);
    kL = 0.1*kLT_kbi;

    R = R_func(kp,kbi,kT,kL,kT_neg,kLT_kbi,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT_kbi);
    RL = RL_func(R,kT,kT_neg, kLT_kbi,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

%     figure(96)
%     loglog(kLT_kbi,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)

    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);

    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);

    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT_kbi(index-1) + fraction*(kLT_kbi(index) - kLT_kbi(index-1));

    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);

    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);

    kLT_10(i) = kLT_kbi(index10-1) + fraction10*(kLT_kbi(index10) - kLT_kbi(index10-1));
    kLT_90(i) = kLT_kbi(index90-1) + fraction90*(kLT_kbi(index90) - kLT_kbi(index90-1));
    ratio(i) = kLT_90(i)./kLT_10(i);
    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);

end

% for i = 1:length(kbi_vals)
%     figure(96)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
% end
% 
% figure(96)
% legend('k_{bi} = k_{bi0}', 'k_{bi} = k_{bi0}*10', 'k_{bi} = k_{bi0}*100','location','southeast')
% legend box off
% box off
% axis tight
% title(ratio)

kLT_half_kbi = kLT_half;
kLT_10_kbi = kLT_10;
kLT_90_kbi = kLT_90;

% do kLneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90 ratio;

kLneg_0 = 0.05;
kLneg_vals = kLneg_0*logspace(-2,2,101); %[kLneg_0 kLneg_0/10 kLneg_0/100]; %kLneg_0*logspace(-2,2,51);
for i=1:length(kLneg_vals)
    kp = 1;
    kT = (1/30); 
    kT_neg = 0.12; 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;  
    kLT = logspace(-8,0,1001); 
    kL_neg = kLneg_vals(i); 
    kbi = kbi_0; 
    kL = 0.1*kLT;

    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);

%     figure(95)
%     loglog(kLT,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)

    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);

    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);

    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));

    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);

    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);

    kLT_10(i) = kLT(index10-1) + fraction10*(kLT(index10) - kLT(index10-1));
    kLT_90(i) = kLT(index90-1) + fraction90*(kLT(index90) - kLT(index90-1));
    ratio(i) = kLT_90(i)./kLT_10(i);
    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);
%     figure(1)
%     %semilogx(kLT,R+T+RL+Rc,color2(1),'linewidth',4)
%     semilogx(kLT,T.*kLT,color2(i),'linewidth',2)
%     hold on
%     semilogx(kLT,R.*kL,color3(i),'linewidth',2)
%     xlabel('ligand concentration')
%     ylabel('Ligand binding flux')
%     box off
%     legend box off
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)
end

% for i = 1:length(kLneg_vals)
%     figure(95)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
%     figure(1)
%     plot(kLT_10(i),0,color4(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),0,color5(i),'Linewidth',2,'markersize',20)
% end
% 
% figure(1)
% legend('T*kLT (bind through tetraspanins)','R*kL (bind outside tetraspanins)','location','northwest')
% 
% figure(95)
% legend('k_{-L} = k_{-L0}', 'k_{-L} = k_{-L0}*10', 'k_{-L} = k_{-L0}*100','location','southeast')
% legend box off
% box off
% axis tight
% title(ratio)

kLT_half_kLneg = kLT_half;
kLT_10_kLneg = kLT_10;
kLT_90_kLneg = kLT_90;

% do kT max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90 ratio;

kT_0 = 1/30;
kT_vals = kT_0*logspace(-2,2,101); %[kT_0 kT_0/50 kT_0/500]; %kT_0*logspace(-2,2,51);
for i=1:length(kT_vals)
    kp = 1;
    kT = kT_vals(i); 
    kT_neg = 0.12; 
    kc =  kc_0; 
    kc_neg = kcneg_0;
    ki = ki_0;  
    kLT = logspace(-8,0,1001);
    kL_neg = kLneg_0;
    kbi = kbi_0; 
    kL = 0.1*kLT;
    
    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kLT_10(i) = kLT(index10-1) + fraction10*(kLT(index10) - kLT(index10-1));
    kLT_90(i) = kLT(index90-1) + fraction90*(kLT(index90) - kLT(index90-1));
    ratio(i) = kLT_90(i)./kLT_10(i);

%     figure(94)
%     loglog(kLT,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)

    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);


end

% for i = 1:length(kT_vals)
%     figure(94)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
% end
% 
% figure(94)
% 
% legend('k_T = k_{T0}', 'k_T = k_{T0}/10', 'k_T = k_{T0}/100','location','southeast')
% legend box off
% box off
% axis tight
% title(ratio)

kLT_half_kT = kLT_half;
kLT_10_kT = kLT_10;
kLT_90_kT = kLT_90;

% do kTneg max and half

clear kp kT kT_neg kc kc_neg ki kLT kL_neg kbi kL kLT_half kLT_10 kLT_90 ratio;

kTneg_0 = 0.12;
kTneg_vals = kTneg_0*logspace(-2,2,101); %[kTneg_0 kTneg_0/10 kTneg_0/100]; %kTneg_0*logspace(-2,2,51);
for i=1:length(kTneg_vals)
    kp = 1;
    kT = kT_0;
    kT_neg = kTneg_vals(i); 
    kc =  kc_0;
    kc_neg = kcneg_0;
    ki = ki_0;
    kLT = logspace(-8,0,1001); 
    kL_neg = kLneg_0; 
    kbi = kbi_0; 
    kL = 0.1*kLT;

    R = R_func(kp,kbi,kT,kL,kT_neg,kLT,kL_neg,kc,kc_neg,ki);
    T = T_func(R,kT,kT_neg,kLT);
    RL = RL_func(R,kT,kT_neg, kLT,kL,kL_neg,ki,kc,kc_neg);
    Rc = Rc_func(kc,kc_neg,ki,RL);
    
    maxsignal(i) = max(RL+Rc);
    index = find((RL+Rc)>(0.5*maxsignal(i)),1);

%     figure(93)
%     loglog(kLT,RL+Rc,'Linewidth',4)
%     hold on
%     xlabel('ligand concentration')
%     ylabel('R_L + R_c')
%     set(gcf,'Color','white')
%     set(gca','fontsize',28)
%     
    
    % define 10 and 90 percents
    index10 = find((RL+Rc)>(0.1*maxsignal(i)),1);
    index90 = find((RL+Rc)>(0.9*maxsignal(i)),1);
    
    previoussignal = RL(index-1) + Rc(index-1);
    nextsignal = RL(index) + Rc(index);
    fraction = (0.5*maxsignal(i) - previoussignal)/(nextsignal - previoussignal);
    kLT_half(i) = kLT(index-1) + fraction*(kLT(index) - kLT(index-1));
    
    % define before and after for 10 and 90
    previoussignal10 = RL(index10-1) + Rc(index10-1);
    nextsignal10 = RL(index10) + Rc(index10);

    previoussignal90 = RL(index90-1) + Rc(index90-1);
    nextsignal90 = RL(index90) + Rc(index90);
    
    % need to define fractions for 10 and 90 percents 
    fraction10 = (0.1*maxsignal(i) - previoussignal10)/(nextsignal10 - previoussignal10);
    fraction90 = (0.9*maxsignal(i) - previoussignal90)/(nextsignal90 - previoussignal90);
    
    kLT_10(i) = kLT(index10-1) + fraction10*(kLT(index10) - kLT(index10-1));
    kLT_90(i) = kLT(index90-1) + fraction90*(kLT(index90) - kLT(index90-1));
    ratio(i) = kLT_90(i)./kLT_10(i);
    % record
    valuelist10(i) = RL(index10)+Rc(index10);
    valuelist90(i) = RL(index90)+Rc(index90);

end

% for i = 1:length(kTneg_vals)
%     figure(93)
%     plot(kLT_10(i),valuelist10(i),color(i),'Linewidth',2,'markersize',20)
%     plot(kLT_90(i),valuelist90(i),color(i),'Linewidth',2,'markersize',20)
% end
% 
% figure(93)
% legend('k_{-T} = k_{-T0}', 'k_{-T} = k_{-T0}/10', 'k_{-T} = k_{-T0}/100','location','southeast')
% legend box off
% box off
% axis tight
% title(ratio)

kLT_half_kTneg = kLT_half;
kLT_10_kTneg = kLT_10;
kLT_90_kTneg = kLT_90;

% plot properties

% figure(1)
% loglog(kc_vals/kc_0,kLT_half_kc,'--','color',[0.5210 0.0860 0.8190],'linewidth',5);
% hold on
% loglog(kcneg_vals/kcneg_0,kLT_half_kcneg,'--','color',[0 0.5 0],'linewidth',5);
% loglog(ki_vals/ki_0,kLT_half_ki,'m--','linewidth',4);
% loglog(kbi_vals/kbi_0,kLT_half_kbi,'--',"Color", [0.5 0.3 0],'linewidth',5);
% loglog(kLneg_vals/kLneg_0,kLT_half_kLneg,'--','color',[0, 0.8078, 0.8235],'linewidth',5);
% loglog(kT_vals/kT_0,kLT_half_kT,'--','color',[1 0.5 0],'linewidth',5);
% loglog(kTneg_vals/kTneg_0,kLT_half_kTneg,'--',"Color",[0.5 0.5 0.5],'linewidth',5);
% 
% set(gca','fontsize',28)
% xlabel('k/k_0','fontsize',30);
% ylabel('k_{+L} at half of max (s^{-1})','fontsize',30);
% %legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{bi}','Uniform : k_{-L}','fontsize',14);
% %legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','k_{+T}','k_{-T}','fontsize',12,location='north');
% set(gcf,'Color','white')
% %axis tight
% %ylim([10^0 10^2])
% xlim([10^-2 10^2])
% ylim([10^-5*0.5 10^-1])
% xticks([10^-2 10^-1 10^0 10^1 10^2])
% yticks([10^-5 10^-4 10^-3 10^-2 10^-1])
% box off

%% plot the localized ratio

figure(1)
semilogx(kc_vals/kc_0,(kLT_90_kc./kLT_10_kc),'--','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
semilogx(kcneg_vals/kcneg_0,(kLT_90_kcneg./kLT_10_kcneg),'--','color',[0 0.5 0],'linewidth',5);
semilogx(ki_vals/ki_0,(kLT_90_ki./kLT_10_ki),'m--','linewidth',4);
semilogx(kbi_vals/kbi_0,(kLT_90_kbi./kLT_10_kbi),'--',"Color", [0.5 0.3 0],'linewidth',5);
semilogx(kLneg_vals/kLneg_0,(kLT_90_kLneg./kLT_10_kLneg),'--','color',[0, 0.8078, 0.8235],'linewidth',5);
semilogx(kT_vals/kT_0,(kLT_90_kT./kLT_10_kT),'--','color',[1 0.5 0],'linewidth',5);
semilogx(kTneg_vals/kTneg_0,(kLT_90_kTneg./kLT_10_kTneg),'--',"Color",[0.5 0.5 0.5],'linewidth',5);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at 90% of max / k_{+L} at 10% of max','fontsize',30);
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{bi}','Uniform : k_{-L}','fontsize',14);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','k_{+T}','k_{-T}','fontsize',12,location='north');
set(gcf,'Color','white')
%axis tight
%ylim([10^0 10^2])
xlim([10^-2 10^2])
%ylim([0 120])
xticks([10^-2 10^0 10^2])
%yticks([0 50 100 120])
title('90 %/10%')
box off


% plot the 10%

figure(2)
loglog(kc_vals/kc_0,(kLT_10_kc),'--','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
loglog(kcneg_vals/kcneg_0,(kLT_10_kcneg),'--','color',[0 0.5 0],'linewidth',5);
loglog(ki_vals/ki_0,(kLT_10_ki),'m--','linewidth',4);
loglog(kbi_vals/kbi_0,(kLT_10_kbi),'--',"Color", [0.5 0.3 0],'linewidth',5);
loglog(kLneg_vals/kLneg_0,(kLT_10_kLneg),'--','color',[0, 0.8078, 0.8235],'linewidth',5);
loglog(kT_vals/kT_0,(kLT_10_kT),'--','color',[1 0.5 0],'linewidth',5);
loglog(kTneg_vals/kTneg_0,(kLT_10_kTneg),'--',"Color",[0.5 0.5 0.5],'linewidth',5);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at 90% of max / k_{+L} at 10% of max','fontsize',30);
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{bi}','Uniform : k_{-L}','fontsize',14);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','k_{+T}','k_{-T}','fontsize',12,location='north');
set(gcf,'Color','white')
%axis tight
%ylim([10^0 10^2])
%xlim([10^-2 10^2])
%ylim([10^-5*0.5 10^-1])
%xticks([10^-2 10^-1 10^0 10^1 10^2])
%yticks([10^-5 10^-4 10^-3 10^-2 10^-1])
title('10 %')
box off


% plot the 90%


figure(3)
loglog(kc_vals/kc_0,(kLT_90_kc),'--','color',[0.5210 0.0860 0.8190],'linewidth',5);
hold on
loglog(kcneg_vals/kcneg_0,(kLT_90_kcneg),'--','color',[0 0.5 0],'linewidth',5);
loglog(ki_vals/ki_0,(kLT_90_ki),'m--','linewidth',4);
loglog(kbi_vals/kbi_0,(kLT_90_kbi),'--',"Color", [0.5 0.3 0],'linewidth',5);
loglog(kLneg_vals/kLneg_0,(kLT_90_kLneg),'--','color',[0, 0.8078, 0.8235],'linewidth',5);
loglog(kT_vals/kT_0,(kLT_90_kT),'--','color',[1 0.5 0],'linewidth',5);
loglog(kTneg_vals/kTneg_0,(kLT_90_kTneg),'--',"Color",[0.5 0.5 0.5],'linewidth',5);

set(gca','fontsize',28)
xlabel('k/k_0','fontsize',30);
ylabel('k_{+L} at 90% of max / k_{+L} at 10% of max','fontsize',30);
%legend('Uniform : k_c','Uniform : k_{-c}','Uniform : k_i','Uniform : k_{bi}','Uniform : k_{-L}','fontsize',14);
%legend('k_{+c}','k_{-c}','k_i','k_{bi}','k_{-L}','k_{+T}','k_{-T}','fontsize',12,location='north');
set(gcf,'Color','white')
%axis tight
%ylim([10^0 10^2])
%xlim([10^-2 10^2])
%ylim([10^-5*0.5 10^-1])
%xticks([10^-2 10^-1 10^0 10^1 10^2])
%yticks([10^-5 10^-4 10^-3 10^-2 10^-1])
title('90 %')
box off


%% 

% Plot colours
% Default uniform params will be blue and default localized will be red
% kT = Orange
% kT_neg = gray
% kc = purple
% kc_neg = green
% kL_neg = cyan
% ki = magenta
% kbi = brown


% % For individual legend curves 
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{+c}', 'Color',[0.5210 0.0860 0.8190]',...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-c}', 'Color','g',...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{i}', 'Color','m',...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{bi}', 'Color',[0.5 0.3 0],...
%            'LineStyle', 'none','FontSize', 20);
% 
% % Make additional text for legend
% annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-L}', 'Color','c',...
%            'LineStyle', 'none','FontSize', 20);
% 
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{+T}', 'Color',[1 0.5 0],...
           'LineStyle', 'none','FontSize', 30);
       
% Make additional text for legend
annotation('textbox', [0.205 0.7 0.15 0.05], 'String', 'k_{-T}', 'Color',[0.5 0.5 0.5],...
           'LineStyle', 'none','FontSize', 30);