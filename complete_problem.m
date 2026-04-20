clc
clear all
close all

%% ------------------ INITIAL PARAMETERS ------------------
% Initial geometry of the tubes
L = 500; % mm assumed rod length
D_init = 30; % mm
d_init = 12; % mm
D_lb_init = 27; % mm
d_lb_init = 13.5; % mm
K = 1;
E = 69;
S = (pi * D_init^2 - pi * d_init^2) / 4; % cross-sectional area
I = (pi / 64) * (D_init^4 - d_init^4); % mm^4 moment of inertia of the section
r = sqrt(I / S);             % mm radius of gyration
Leff = K * L;                % mm effective length

% Critical buckling load (Euler)
F_cr = (pi^2 * E * I) / (Leff^2);   % N

% Critical buckling stress
sigma_cr = F_cr / S;         % MPa

% Material (Magnesium)
pho = 0.0027; % density in g/mm^3
E = 69; % GPa
sigma_max = 30; % MPa

% Initial coordinates of the nodes
x2_init = -500; y2_init = 500;
x3_init = 500;  y3_init = 0;
x4_init = 0;    y4_init = 50-6000;

nodes_init = [
    0, 0;
    x2_init, y2_init;
    x3_init, y3_init;
    x4_init, y4_init
];

% Elements (rods)
elements = [
    1, 2;
    2, 4;
    1, 3;
    1, 4;
    3, 4
];

% Loads on the nodes (Fx, Fy)
forces = [
    0,   0;
    0, 400;
    0, 400;
    0, -800
];

%% ------------------ OPTIMIZATION ------------------

% Initial values of the variables to optimize
x0 = [D_init, d_init, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init];

% Limits
D_lb_min = 10; D_lb_max = 40;
d_lb_min = 9; d_lb_max = 100;
D_min = 10; D_max = 40;
d_min = 9;  d_max = D_max - 0.5;
x2_min = -600; x2_max = -400;
y2_min = 400;  y2_max = 600;
x3_min = 400;  x3_max = 600;
y3_min = -100; y3_max = 100;
x4_min = -100; x4_max = 100;
y4_min = 400;  y4_max = 600;

lb = [D_min, d_min, x2_min, y2_min, x3_min, y3_min, x4_min, y4_min];
ub = [D_max, d_max, x2_max, y2_max, x3_max, y3_max, x4_max, y4_max];

% Functions
mass_function = @(x) compute_mass(x, elements, pho);
stress_constraint = @(x) compute_stress(x, elements, forces, sigma_max);
mc_function = @(x) compute_mass_center(x, elements, pho);

[~,~,sigma_massimo, forze_interne] =  compute_stress(x0, elements, forces, sigma_max);

% Initial values used for scaling in the objective function
y_mc_init = mc_function(x0);
m_tot_init = mass_function(x0);

% Objective function definition
obj_function = @(x) mass_function(x); %+ mc_function(x)/500;

% Optimization with fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_opt, f_opt, ~, ~, lambda_fmincon] = fmincon(obj_function, x0, [], [], [], [], lb, ub, stress_constraint, options);

% Optimization with global algorithm
options = optimoptions('ga','Display','iter','MaxGenerations',100);
[x_ga, fval_ga, ~, ~, lambda_ga] = ga(obj_function, length(x0), [], [], [], [], lb, ub, stress_constraint, options);

% Optimization with fmincon using ga as initial guess
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_fga, f_opt_fga] = fmincon(obj_function, x_ga, [], [], [], [], lb, ub, stress_constraint, options);

% test feasibility of ga point
c_ga = stress_constraint(x_ga);

c_fmincon = stress_constraint(x_opt);
%% --------- MULTISTART TEST WITH INITIALIZATION GRID ---------

D_vals_grid = linspace(15, 40, 6);   % Possible values for D (e.g. 6 points between 15 and 40 mm)
d_vals_grid = linspace(9, 38.5, 6);  % Possible values for d (e.g. 6 points between 9 and 38.5 mm)

results = [];
counter = 1;

for i = 1:length(D_vals_grid)
    for j = 1:length(d_vals_grid)
        D0 = D_vals_grid(i);
        d0 = d_vals_grid(j);
        
        % Satisfies the constraint D - d >= 1.5
        if D0 - d0 < 1.5
            continue
        end

        % Build x0 with new initial values
        x0_grid = x0;
        x0_grid(1) = D0;
        x0_grid(2) = d0;

        % Run fmincon
        try
            [x_temp, f_temp] = fmincon(obj_function, x0_grid, [], [], [], [], lb, ub, stress_constraint, options);

            results(counter).D0 = D0;
            results(counter).d0 = d0;
            results(counter).fval = f_temp;
            results(counter).x = x_temp;
            counter = counter + 1;
        catch ME
            % Skip in case of numerical errors or unsatisfied constraints
            disp(['Error with D = ', num2str(D0), ', d = ', num2str(d0), ': ', ME.message]);
        end
    end
end

% Visualize results in 3D
D0_vec = [results.D0];
d0_vec = [results.d0];
fval_vec = [results.fval];

figure(1);
scatter3(D0_vec, d0_vec, fval_vec, 100, fval_vec, 'filled');
xlabel('Initial D');
ylabel('Initial d');
zlabel('Objective function value');
title('fmincon convergence from different initial points');
colorbar;
grid on;

% 10 DIFFERENT INITIAL GUESSES

% Number of design variables
n_vars = 8;

% Number of initial guesses
n_guesses = 10;
% Preallocate 10x8 matrix
x0_all = zeros(n_guesses, n_vars);

% Random generation
for i = 1:n_guesses
    x0_all(i, :) = lb + rand(1, n_vars) .* (ub - lb);
end

% Display the results
disp('Initial guesses x0:');
disp(x0_all);

% TEST WITH DIFFERENT INITIAL GUESSES FOR FMINCON
x_opt_fmincon = zeros(n_guesses,8);
f_opt_fmincon = zeros(n_guesses,1);

for i = 1:n_guesses
    [x_opt_fmincon(i,:), f_opt_fmincon(i,:)] = fmincon(obj_function, x0_all(i,:), [], [], [], [], lb, ub, stress_constraint, options);
end

% TEST WITH DIFFERENT GENETIC ALGORITHM RUNS
x_opt_ga = zeros(n_guesses,8);
f_opt_ga = zeros(n_guesses,1);

for i = 1:n_guesses
    options = optimoptions('ga','Display','iter','MaxGenerations',100);
    [x_opt_ga(i,:), f_opt_ga(i,:)] = ga(obj_function, 8, [], [], [], [], lb, ub, stress_constraint, options);
end

% TEST WITH GA FOLLOWED BY FMINCON

x_opt_ga_fmincon = zeros(n_guesses,8);
f_opt_ga_fmincon = zeros(n_guesses,1);

for i = 1:n_guesses
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    [x_opt_ga_fmincon(i,:), f_opt_ga_fmincon(i,:)] = fmincon(obj_function, x_opt_ga(i,:), [], [], [], [], lb, ub, stress_constraint, options);
end



%% ------------------ RESULTS ------------------

D_opt = x_opt(1);
d_opt = x_opt(2);
x2_opt = x_opt(3); y2_opt = x_opt(4);
x3_opt = x_opt(5); y3_opt = x_opt(6);
x4_opt = x_opt(7); y4_opt = x_opt(8);

disp(['Objective function fmincon: ', num2str(f_opt)])
disp(['Objective function genetic algorithm: ', num2str(fval_ga)])
disp(['Objective function fmincon and genetic algorithm: ', num2str(f_opt_fga)])
disp(['Optimized outer diameter: ', num2str(D_opt), ' mm']);
disp(['Optimized inner diameter: ', num2str(d_opt), ' mm']);
disp(['Node 2: (', num2str(x2_opt), ', ', num2str(y2_opt), ')']);
disp(['Node 3: (', num2str(x3_opt), ', ', num2str(y3_opt), ')']);
disp(['Node 4: (', num2str(x4_opt), ', ', num2str(y4_opt), ')']);

% Plot optimized frame
figure(2);
nodes_opt = [0, 0; x2_opt, y2_opt; x3_opt, y3_opt; x4_opt, y4_opt];
hold on;
for i = 1:size(elements,1)
    n1 = elements(i,1);
    n2 = elements(i,2);
    plot([nodes_opt(n1,1), nodes_opt(n2,1)], ...
         [nodes_opt(n1,2), nodes_opt(n2,2)], ...
         'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end
grid on;
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
hold off;


%% ----------- 2D GRAPHICAL ANALYSIS ------------

% Here I need to plot the trend of the objective function with respect to 
% individual variables, it shouldn't be difficult

D_vals = D_min:1:D_max;
d_vals = d_min:1:d_max;
y2_vals = y2_min:1:y2_max;
y3_vals = y3_min:1:y3_max;
y4_vals = y4_min:1:y4_max;
x2_vals = x2_min:1:x2_max;
x3_vals = x3_min:1:x3_max;
x4_vals = x4_min:1:x4_max;

% plot with respect to the inner diameter
% first set the maximum value for the outer diameter to avoid negative mass

D_d = D_max;
obj = zeros(length(d_vals),1);
for i = 1:length(d_vals)
    x_d = [D_d, d_vals(i), x2_init, y2_init, x3_init, y3_init, x4_init, y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_d);
end
figure(3)
plot(d_vals,obj, 'LineWidth',2.5);

% plot with respect to the outer diameter

obj = zeros(length(D_vals),1);
for i = 1:length(D_vals)
    x_D = [D_vals(i), d_min, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_D);
end
figure(4)
plot(D_vals,obj,'LineWidth',2.5);

% plot with respect to y2

obj = zeros(length(y2_vals),1);
for i = 1:length(y2_vals)
    x_y2 = [D_max, d_max, x2_init, y2_vals(i), x3_init, y3_init, x4_init, y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_y2);
end
figure(5)
plot(y2_vals,obj,'LineWidth',2.5);

% plot with respect to x2

obj = zeros(length(x2_vals),1);
for i = 1:length(x2_vals)
    x_x2 = [D_max, d_max, x2_vals(i), y2_init, x3_init, y3_init, x4_init, y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_x2);
end
figure(6)
plot(x2_vals,obj, 'LineWidth',2.5);

% plot with respect to y3

obj = zeros(length(y3_vals),1);
for i = 1:length(y3_vals)
    x_y3 = [D_max, d_max, x2_init, y2_init, x3_init, y3_vals(i), x4_init, y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_y3);
end
figure(7)
plot(y3_vals,obj,'LineWidth',2.5);

% plot with respect to x3

obj = zeros(length(x3_vals),1);
for i = 1:length(x3_vals)
    x_x3 = [D_max, d_max, x2_init, y2_init, x3_vals(i), y3_init, x4_init, y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_x3);
end
figure(8)
plot(x3_vals,obj,'LineWidth',2.5);

% plot with respect to y4

obj = zeros(length(y4_vals),1);
for i = 1:length(y4_vals)
    x_y4 = [D_max, d_max, x2_init, y2_init, x3_init, y3_init, x4_init, y4_vals(i), D_lb_init, d_lb_init];
    obj(i) = mass_function(x_y4);
end
figure(9)
plot(y4_vals,obj, 'LineWidth',2.5);

% plot with respect to x4

obj = zeros(length(x4_vals),1);
for i = 1:length(x4_vals)
    x_x4 = [D_max, d_max, x2_init, y2_init, x3_init, y3_init, x4_vals(i), y4_init, D_lb_init, d_lb_init];
    obj(i) = mass_function(x_x4);
end
figure(10)
plot(x4_vals,obj, 'LineWidth',2.5);



%% ----------- ANALYSIS OF FINITE DIFFERENCE GRADIENTS ------------

% compute the value of the partial derivative of the objective function for
% each design variable over different step sizes

hx = logspace(-20,0,100); % vector of finite difference steps
for i=1:1:length(hx)

  % Finite difference step
  hxi = hx(i);

  % Objective function
  fx = mass_function(x0);
  fx1plush = mass_function([x0(1)+hxi, x0(2), x0(3), x0(4), x0(5), x0(6), x0(7), x0(8)]);
  fx2plush = mass_function([x0(1), x0(2)+hxi, x0(3), x0(4), x0(5), x0(6), x0(7), x0(8)]);
  fx3plush = mass_function([x0(1), x0(2), x0(3)+hxi, x0(4), x0(5), x0(6), x0(7), x0(8)]);
  fx4plush = mass_function([x0(1), x0(2), x0(3), x0(4)+hxi, x0(5), x0(6), x0(7), x0(8)]);
  fx5plush = mass_function([x0(1), x0(2), x0(3), x0(4), x0(5)+hxi, x0(6), x0(7), x0(8)]);
  fx6plush = mass_function([x0(1), x0(2), x0(3), x0(4), x0(5), x0(6)+hxi, x0(7), x0(8)]);
  fx7plush = mass_function([x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7)+hxi, x0(8)]);
  fx8plush = mass_function([x0(1), x0(2), x0(3), x0(4), x0(5), x0(6), x0(7), x0(8)+hxi]);
  dfdx1(i) = (fx1plush - fx)/hxi;
  dfdx2(i) = (fx2plush - fx)/hxi;
  dfdx3(i) = (fx3plush - fx)/hxi;
  dfdx4(i) = (fx4plush - fx)/hxi;
  dfdx5(i) = (fx5plush - fx)/hxi;
  dfdx6(i) = (fx6plush - fx)/hxi;
  dfdx7(i) = (fx7plush - fx)/hxi;
  dfdx8(i) = (fx8plush - fx)/hxi;
end

figure(11)
semilogx(hx, dfdx1, 'LineWidth', 2.5);

figure(12)
semilogx(hx, dfdx2, 'LineWidth', 2.5);

figure(13)
semilogx(hx, dfdx3, 'LineWidth', 2.5);

figure(14)
semilogx(hx, dfdx4, 'LineWidth', 2.5);

figure(15)
semilogx(hx, dfdx5, 'LineWidth', 2.5);

figure(16)
semilogx(hx, dfdx6, 'LineWidth', 2.5);

figure(17)
semilogx(hx, dfdx7, 'LineWidth', 2.5);

figure(18)
semilogx(hx, dfdx8, 'LineWidth', 2.5);




%% ------------------ FUNCTIONS ------------------

function m = compute_mass(x, elements, pho)
    D = x(1); d = x(2);
     
    S = (pi * D.^2 - pi * d.^2) / 4; % Cross-sectional area

    nodes = [0, 0; x(3), x(4); x(5), x(6); x(7), x(8)];
    lengths = zeros(size(elements,1),1);
    masses = zeros(size(elements,1),1);
    for i = 1:size(elements,1)
        n1 = elements(i,1);
        n2 = elements(i,2);
        p1 = nodes(n1, :);
        p2 = nodes(n2, :);
        lengths(i) = norm(p2 - p1);
        masses(i) = abs(lengths(i)) * pho * S;
    end

    m = sum(masses); % Total mass
end

function [c, ceq, sigma_maximum, internal_forces] = compute_stress(x, elements, forces, sigma_max)
    D = x(1); d = x(2);
    D_max = 40;
    D_lb_max = 40;
    x_saddle = x(7);
    nodes = [0, 0; x(3), x(4); x(5), x(6); x(7), x(8)];
    S = (pi * D^2 - pi * d^2) / 4;
    dof = size(nodes,1) * 2;
    A = zeros(dof, size(elements,1));
    L_vect = zeros(size(elements,1),1); % Element lengths

    for i = 1:size(elements,1)
        n1 = elements(i,1);
        n2 = elements(i,2);
        x1 = nodes(n1,1); y1 = nodes(n1,2);
        x2 = nodes(n2,1); y2 = nodes(n2,2);
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        L_vect(i) = L;
        cosT = (x2 - x1) / L;
        sinT = (y2 - y1) / L;

        A(2*n1-1, i) = cosT;
        A(2*n1,   i) = sinT;
        A(2*n2-1, i) = -cosT;
        A(2*n2,   i) = -sinT;
    end

    F = reshape(forces', [], 1);
    T = A \ F; % Element forces (tensions)
    internal_forces = T;
    sigma = T ./ S;
    sigma_maximum = sigma;

    % Constraints
    stress_violation = abs(max(sigma)) - sigma_max;
    int_geometry_violation = d - D + 0.1;
    saddle_violation = -x_saddle + 50;

    % Buckling constraint (only on compressive elements)
    E = 210e9 * 10^(-6); % Young's modulus in N/mm^2
    I = (pi / 64) * (D^4 - d^4); % Moment of inertia for hollow tube
    K = 1; % End condition coefficient (pinned-pinned)
    P_cr = (pi^2 * E * I) ./ (K^2 * max(L_vect).^2); % Critical buckling load
    buckling_violation = abs(max(T)) - P_cr;

    % Strain constraint: ε <= 0.0001 (0.01%)
    strain = T ./ (E * S);
    strain_violation = max(abs(strain)) - 1e-4;

    % All inequality constraints
    c = [stress_violation; int_geometry_violation; saddle_violation; strain_violation; buckling_violation];
    ceq = []; % No equality constraints
end

function [y_cm] = compute_mass_center(x, elements, pho)

    D = x(1); d = x(2);
    S = (pi * D.^2 - pi * d.^2) / 4;

    nodes = [0, 0; x(3), x(4); x(5), x(6); x(7), x(8)];
    lengths = zeros(size(elements,1),1);
    masses = zeros(size(elements,1),1);
    y_cm_s = zeros(size(elements,1),1);

    % Compute mass of each element
    for i = 1:size(elements,1)
        n1 = elements(i,1);
        n2 = elements(i,2);
        p1 = nodes(n1, :);
        p2 = nodes(n2, :);
        lengths(i) = norm(p2 - p1);
        masses(i) = lengths(i) * pho * S;
        y_cm_s(i) = (p1(2) + p2(2)) / 2; % Midpoint y-coordinate
    end

    % Weighted average to find center of mass
    y_cm = sum(masses .* y_cm_s) / sum(masses); 
end



