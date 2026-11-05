function [peakheight,duration valuelist, t_rec] = EvolveEGFR(kp,kLT,kT,kT_neg,kc,kc_neg,ki,kbi,kL,kL_neg,...
    R0,T0,RL0,Rc0,R_new, T_new, RL_new, Rc_new,Delta_t,threshold)

% Equation [1]: dR/dt
f = @(R,RL,T,Rc,kp,kbi,kT,kT_neg,kL,kL_neg) kp - kbi * R - kT * R + kT_neg * T - kL * R + kL_neg*RL;
% Equation [2]: dT/dt
f1 = @(R,T,kT,kT_neg,kLT) kT * R - kT_neg *T - kLT *T;
% Equation [3]: dRL/dt
f2 = @ (kLT,T,kL,R,kc_neg,Rc,kL_neg,RL,kc)  kLT * T + kL * R + kc_neg *Rc - kL_neg*RL -kc*RL;
% Equation [4]: dRc/dt
f3 = @ (kc,RL,kc_neg,Rc,ki) kc * RL - kc_neg *Rc - ki *Rc;

R = R0;
T = T0;
RL = RL0;
Rc = Rc0;

foundstoppinglevel = 0;
foundpeak = 0;
t = 0;
counter =1;
while(foundstoppinglevel == 0)
    if(t < (0.1*Delta_t))
        last = (RL + Rc) ;
    else
        lastlast = last;
        last = (RL + Rc);
    end

    dR = f(R,RL,T,Rc,kp,kbi,kT,kT_neg,kL,kL_neg);
    dT = f1(R,T,kT,kT_neg,kLT);
    dRL = f2(kLT,T,kL,R,kc_neg,Rc,kL_neg,RL,kc);
    dRc = f3(kc,RL,kc_neg,Rc,ki);

    R = R + dR*Delta_t;
    T = T + dT*Delta_t;
    RL = RL + dRL*Delta_t;
    Rc = Rc + dRc*Delta_t;
    t_rec(counter) = t;
    t = t + Delta_t;
        
    if(t > (1.1*Delta_t))
        now = (RL + Rc);
        if((last > lastlast) && (last > now))
            foundpeak = 1;
            peakval = now;
            peaktime = t;
            peakheight = peakval - (RL0 + Rc0);
            stoppingval = threshold*(peakval - (RL_new + Rc_new)) + (RL_new + Rc_new);
        end
    end

    if(foundpeak == 1)
        now = RL + Rc;
        if(now < stoppingval)
            foundstoppinglevel = 1;
            duration = t - peaktime;
        end
    end

    valuelist(counter) = RL + Rc;
    counter = counter +1;
end

end