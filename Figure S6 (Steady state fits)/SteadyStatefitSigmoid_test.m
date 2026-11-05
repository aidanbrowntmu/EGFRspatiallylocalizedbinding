function [params, fitted_y] = SteadyState_Sigmoid_fit(kL, y)

    % define the ligand concentration
    kL = transpose(kL);     % Make sure data is in column
    % define the Steady state values (RL+Rc) 
    y  = transpose(y);      % Make sure data is in column

    % ignore any negative ligand concentrations (meaningless)
    negative = kL<0;
    positive = kL>0;

    % define range of kL and y for fit
    kL_fit = kL(positive);
    y_fit  = y(positive);

    % log transform ligand concentration
    x = log10(kL_fit);

    % Initial conditions based on max and min of sigmoid
    ymin = min(y_fit);
    ymax = max(y_fit);
    
    % Initial condition for ligand required for half max (signaling onset)
    y_half = (ymin + ymax) / 2; % take the average as the mid point
    [~, half_index] = min(abs(y_fit - y_half));

    % store the estimate for ligand required for half max
    X0 = x(half_index); 
    
    % Initial condition for slope, s
    s = 1.0; 

    % Vector of initial conditions
    parameters = [ymin, ymax, X0, s];

    % Define equation of line fit: y(X) = ymin + (ymax - ymin)/(1 + exp(-(X - X0)/s))

    Function = @(p, X_val) p(1) + (p(2) - p(1)) ./ (1 + exp(-(X_val - p(3)) ./ p(4)));

    % Define cost function that minimizes fit using square of residuals
    CostFunction  = @(p) sum((y_fit - Function(p, x)).^2);

    % Use fminsearch to iterate until yfit is minimized
    options = optimset('Display', 'off', 'MaxFunEvals', 10000, 'MaxIter', 10000);
    p_opt = fminsearch(CostFunction, parameters, options);

    % store fit of y
    fitted_y(positive) = Function(p_opt, x);


    % Display output parameters
    %params.y_min    = p_opt(1);         % Sigmoid minimum
    %params.y_max    = p_opt(2);         % Sigmoid maximum
    params.log10_x0 = p_opt(3);         
    params.kL_0     = 10^(p_opt(3));    % log transform ligand required for half max
    %params.s        = p_opt(4);         % slope

    % Output the ligand requried for half max and the fit of y from the function
end