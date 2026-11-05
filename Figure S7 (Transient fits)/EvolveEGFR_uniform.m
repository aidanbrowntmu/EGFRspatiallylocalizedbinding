function [peakheight,duration, valuelist, t_rec] = EvolveEGFR_uniform(kp,kc,kc_neg,ki,kbi,kL,kL_neg,...
    R0,RL0,Rc0,R_new, RL_new, Rc_new,Delta_t,threshold)

% Equation [9]
f = @(R,RL,kp,kbi,kL,kL_neg) kp +(RL *kL_neg) - R*(kL+kbi); % GOOD
% Equation [10]
f1 = @ (kL,R,RL,kc_neg,Rc,kL_neg,kc) R*kL -RL*(kL_neg+kc) + Rc*kc_neg;
% Equation [11]
f2 = @ (kc,RL,kc_neg,Rc,ki) RL*kc - Rc*(kc_neg+ki);


R = R0;
RL = RL0;
Rc = Rc0;

foundstoppinglevel = 0;
foundpeak = 0;
t = 0;
counter =1; 
while(foundstoppinglevel == 0)
    if(t < (0.1*Delta_t))
        last = RL + Rc;
    else
        lastlast = last;
        last = RL + Rc;
    end

    dR = f(R,RL,kp,kbi,kL,kL_neg);
    dRL = f1(kL,R,RL,kc_neg,Rc,kL_neg,kc);
    dRc = f2(kc,RL,kc_neg,Rc,ki);

    R = R + dR*Delta_t;
    RL = RL + dRL*Delta_t;
    Rc = Rc + dRc*Delta_t;
    t_rec(counter) = t;
    t = t + Delta_t;

    if(t > (1.1*Delta_t))
        now = RL + Rc;
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

    valuelist(counter) = RL+Rc;
    counter = counter +1;

end
end