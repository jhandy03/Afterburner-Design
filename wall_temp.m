clear; clc;

% --- Geometry ---
ri = 0.071 / 2;         % inner radius [m]
ro = ri + 0.005;        % outer radius [m], 5mm wall
N  = 10;                % radial nodes

% --- Material: 316 SS ---
rho = 8000;             % density [kg/m3]
cp  = 500;              % specific heat [J/kgK]
% k is temperature-dependent - evaluated each step (from Gemini, good idea)
% k(T) = 0.015*T + 8.5  [W/mK], valid ~300-1100K for 316SS

% ---- Conservative Gas-Side Assumptions ---
h_gas     = 800;        % worst-case forced convection coefficient [W/m2K]
T_gas_dry = 1000;       % dry operation gas temp [K]
T_gas_wet = 2200;       % wet operation gas temp [K]
dTdt_threshold = 80;    % outer wall heating rate that triggers wet assumption [K/s]

% --- Stability ---
r_nodes = linspace(ri, ro, N);
dr      = r_nodes(2) - r_nodes(1);

% Use worst-case (highest) k for stability calculation
% Higher k = higher alpha = more restrictive timestep requirement
k_max   = 0.015 * 2200 + 8.5;      % ~41.5 W/mK at max temp
alpha_max = k_max / (rho * cp);
dt_max  = dr^2 / (2 * alpha_max);
fprintf('Max stable timestep: %.4f s (%.1f ms)\n', dt_max, dt_max*1000);

dt = dt_max * 0.8;      % run at 80% of stability limit for margin

% --- Initialize ---
T_wall        = ones(1, N) * 300;   % start at ambient [K]
T_outer_prev  = 300;
t             = 0;

% In real system: replace this loop with your control system timer callback
% and read T_outer from your DAQ each iteration

%--- Simulation Loop (replace with real-time callback) ---
t_end         = 15;     % simulate 30 seconds
n_steps       = floor(t_end / dt);

% Preallocate logging
log_t       = zeros(1, n_steps);
log_Tinner  = zeros(1, n_steps);
log_Touter  = zeros(1, n_steps);
log_Tgas    = zeros(1, n_steps);

for step = 1:n_steps

    % === IN REAL SYSTEM: REPLACE THIS BLOCK WITH DAQ READ ===
    % Simulated outer wall temperature - dry then wet transient at t=10s
    if t < 10
        T_outer = 600 + 10*sin(0.5*t);     % dry operation, some noise
    else
        T_outer = 600 + (t-10)*80;          % outer wall climbing post-ignition
        T_outer = min(T_outer, 1050);        % clamp to physical limit
    end
    % ==========================================================

    % --- Infer operating mode from outer wall rate of change ---
    dTdt = (T_outer - T_outer_prev) / dt;

    if dTdt > dTdt_threshold
        T_gas = T_gas_wet;      % wall heating fast - assume wet, worst case
    else
        T_gas = T_gas_dry;
    end

    % --- Update temperature-dependent properties ---
    T_avg   = mean(T_wall);
    k_val   = 0.015 * T_avg + 8.5;         % 316SS k(T) [W/mK]
    alpha   = k_val / (rho * cp);

    % --- Explicit finite difference: interior nodes ---
    T_new = T_wall;

    for i = 2:N-1
        d2T = (1 + dr/(2*r_nodes(i))) * T_wall(i+1) ...
            - 2 * T_wall(i) ...
            + (1 - dr/(2*r_nodes(i))) * T_wall(i-1);
        T_new(i) = T_wall(i) + alpha * dt/dr^2 * d2T;
    end

% --- Inner boundary: convective flux from gas ---
    T_new(1) = T_wall(1) + 2*alpha*dt/dr * ...
        ( h_gas/k_val * (T_gas - T_wall(1)) ...
        - (T_wall(1) - T_wall(2))/dr );

    % --- Outer boundary: pin to measurement ---
    T_new(N) = T_outer;

    % --- Sanity check ---
    if any(isnan(T_new)) || any(T_new < 200) || any(T_new > 3500)
        warning('Step %d: unphysical temperatures - holding previous state', step);
        T_new = T_wall;     % hold last good state, do not update
    end

    T_wall       = T_new;
    T_inner      = T_wall(1);
    T_outer_prev = T_outer;
    t            = t + dt;

    % --- Log ---
    log_t(step)      = t;
    log_Tinner(step) = T_inner;
    log_Touter(step) = T_outer;
    log_Tgas(step)   = T_gas;

    % --- Safety checks (material limits for 316 SS) ---
    if T_inner > 1100
        warning('Step %d | t=%.2fs: CRITICAL - T_inner %.0fK exceeds structural limit', ...
                step, t, T_inner);
    end
    if T_inner > 1645
        error('Step %d | t=%.2fs: DANGER - T_inner %.0fK exceeds melt point. SHUTDOWN.', ...
              step, t, T_inner);
    end

end

% --- Plot results ---
figure;
plot(log_t, log_Tinner, 'r-', 'LineWidth', 2); hold on;
plot(log_t, log_Touter, 'b-', 'LineWidth', 2);
plot(log_t, log_Tgas,   'k--', 'LineWidth', 1);
yline(1100, 'r:', 'Structural Limit', 'LineWidth', 1.5);
yline(1645, 'm:', 'Melt Point',       'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Temperature [K]');
title('Afterburner Wall Temperature - Inner vs Outer');
legend('T inner (estimated)', 'T outer (measured)', 'T gas (assumed)', ...
       'Location', 'northwest');
grid on;