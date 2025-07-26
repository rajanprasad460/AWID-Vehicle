function dy = carDynamics_ODE(t, y, tau_func, params, reset)

    persistent lastTime

    % Reset persistent variables if requested
    if nargin > 4 && reset
        lastTime = [];
        return
    end

    % If called with no inputs, return logs
    if nargin == 0
        dy = [];  % First output still required
        return
    end

    % Normal operation during ODE solving
    if isempty(lastTime)
        lastTime = tic;
    end

    [dy, ~, ~] = carDynamics(t, y, tau_func, params);

    if toc(lastTime) > 1
        fprintf('Time: %.2f s | Elapsed Time: %.2f s\n', t, toc(lastTime));
        lastTime = tic;
    end

end
