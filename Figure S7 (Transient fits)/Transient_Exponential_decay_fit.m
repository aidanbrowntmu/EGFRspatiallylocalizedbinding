function [params, fitted_y] = Transient_Exponential_decay_fit(t, y)

    % define the decay timescale
    t = transpose(t);       % Make sure data is in column
    % define the Steady state values (RL+Rc) 
    y  = transpose(y);      % Make sure data is in column

    % define peak concentration and time
    [peak_IC, max_index] = max(y);
    t_peak = t(max_index);
    
    % ignore time points prior to peak (meaningless for fitting)
    negative = t <= t_peak;     % data prior to peak
    positive = t >= t_peak;     % data following peak (fitted portion)

    % Set our decay parameters to data ONLY after peak
    t_decay = t(positive);
    y_decay = y(positive);

    % Initial condition for new steady state B (new steady state following peak decay)
    B = min(y_decay);
    
    % Initial condition for k 
    dt = t_decay - t_peak; % change in t in exponential

    % check if data is within peak decay range (concentration level and time)

    decay_index = (y_decay - B) > 0 & dt > 0;

    % solve exponential decay function (y(X) = B + (A-B) * exp(-k*(t-t0)))
    % for k and use as IC:
    k = -mean(log((y_decay(decay_index) - B) / (peak_IC - B)) ./ dt(decay_index));
    
    % check that the IC for k is a positive number
    if isnan(k) || k <= 0
        k = 0.1; 
    end

    % Vector of initial conditions
    parameters = [B, peak_IC, k]; 

    % Define equation of line fit: y(X) = B + (A-B) * exp(-k*(t-t0))
    Function = @(p, t_val) p(1) + (p(2) - p(1)) .* exp(-p(3) .* (t_val - t_peak));

    % Define cost function that minimizes fit using square of residuals
    costFunction = @(p) sum((y_decay - Function(p, t_decay)).^2);

    % Use fminsearch to iterate until yfit is minimized
    options = optimset('Display', 'off', 'MaxFunEvals', 5000, 'MaxIter', 5000);
    p_opt = fminsearch(costFunction, parameters, options);

    % Display output parameters
    %params.B = p_opt(1);        % New steady state concentration
    %params.A = p_opt(2);        % Max peak value
    params.k = p_opt(3);         % Decay constant
    %params.t_peak = t_peak;     % time to reach max peak
    %params.half_life = log(2) / p_opt(3);   

    % store fit of y
    fitted_y(positive) = Function(p_opt, t_decay);

    % Output the peak decay timescale and the fit of y from the function

end