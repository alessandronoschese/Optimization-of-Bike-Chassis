clc
clear all
close all

% SIMPLIFIED PROBLEM — only varying y-coordinates of end nodes
% Initial tubular geometry
D = 25; % mm (outer diameter)
d = 12; % mm (inner diameter)

% Material properties (Magnesium)
pho = 0.0027; % density in g/mm^3
E = 69;       % GPa (Young's modulus)
sigma_max = 150; % MPa (max allowable stress)

% Initial coordinates of nodes
x2_init = -500; y2_init = 500;
x3_init = 500;  y3_init = 0;
x4_init = 0;    y4_init = 500;

nodes_init = [
    0, 0;
    x2_init, y2_init;
    x3_init, y3_init;
    x4_init, y4_init
];

% Bar elements (connectivity matrix)
elements = [
    1, 2;
    2, 4;
    1, 3;
    1, 4;
    3, 4
];

% Node forces (Fx, Fy) in N
forces = [
    0,   0;
    0, 400;
    0, 400;
    0, -800
];

%% ------------------ OPTIMIZATION ------------------
% Only varying outer and inner diameters
x0 = [D, d];

% Design bounds
D_min = 15;  D_max = 35;
d_min = 8;   d_max = 28;

alpha = 50000; % penalty factor (used in objective)

lb = [D_min, d_min];
ub = [D_max, d_max];

% Constraint function: includes stress, geometry, etc.
stress_constraint = @(x) compute_stress(x, elements, forces, sigma_max, ...
    x2_init, y2_init, x3_init, y3_init, x4_init, y4_init);

% Objective function with buckling penalty
mass_buckling_obj = @(x) objective_with_buckling(x, elements, pho, E, ...
    x2_init, y2_init, x3_init, y3_init, x4_init, y4_init, alpha);

% Objective function with penalty (possibly alternative)
obj_with_pen = @(x) objective_with_penalty(x, elements, pho, E, forces,...
    x2_init, y2_init, x3_init, y3_init, x4_init, y4_init);

% Local optimization with fmincon
options = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'sqp', ...
    'ConstraintTolerance', 1e-20);

[x_opt, f_opt] = fmincon(mass_buckling_obj, x0, [], [], [], [], ...
    lb, ub, stress_constraint, options);

% Post-processing: extract final stress and internal force results
[~,~,sigma_massimo, forze_interne] =  compute_stress(x_opt, ...
    elements, forces, sigma_max, ...
    x2_init, y2_init, x3_init, y3_init, x4_init, y4_init);

% Global optimization using Genetic Algorithm (GA)
options = optimoptions('ga','Display','iter','MaxGenerations',100);
[x_ga, fval_ga] = ga(mass_buckling_obj, length(x0), [], [], [], [], ...
    lb, ub, stress_constraint, options);

% Use GA result as initial guess for refined fmincon optimization
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
[x_fga, f_opt_fga] = fmincon(mass_buckling_obj, x_ga, [], [], [], [], ...
    lb, ub, stress_constraint, options);



%% ----------- GRIGLIA VALORI INIZIALI ------------
D_x_zeros = 22:1:30;
d_x_zeros = 8:1:16;

[DX, dX] = meshgrid(D_x_zeros, d_x_zeros); % crea la griglia
Z = zeros(size(DX)); % per salvare i risultati finali

for i = 1:numel(DX)
    x0 = [DX(i), dX(i)];
    [x_opt_grid, fval] = fmincon(mass_buckling_obj, x0, [], [], [], [], lb, ub, stress_constraint, options);
    Z(i) = fval; % salva valore ottimo trovato
end

figure
scatter3(DX(:), dX(:), Z(:), 100, Z(:), 'filled') % Grafico a dispersione 3D
colorbar;
colormap(turbo)
xlabel('Starting D');
ylabel('Starting d');
zlabel('fmincon Value');
grid on;
%% ----------- GRAPHIC ANALYSIS ------------


D_vals = D_min:0.1:D_max;
d_vals = d_min:0.1:d_max;
obj = zeros(length(d_vals), length(D_vals));
obj_p = zeros(length(d_vals), length(D_vals));
constraints = zeros(length(D_vals),length(d_vals),4);
constraint_stress = zeros(length(D_vals),length(d_vals));
constraint_geom = zeros(length(D_vals),length(d_vals));
constraint_buckling = zeros(length(D_vals),length(d_vals));
constraint_def = zeros(length(D_vals),length(d_vals));
for j = 1:length(d_vals)
    for i = 1:length(D_vals)
        x = [D_vals(i), d_vals(j)];
        obj(j,i) = mass_buckling_obj(x);
        obj_p(j,i) = obj_with_pen(x);
        constraints(j,i,:) = stress_constraint(x);
        constraint_stress(j,i) = constraints(j,i,1); 
        constraint_buckling(j,i) = constraints(j,i,2); 
        constraint_geom(j,i) = constraints(j,i,3); 
        constraint_def(j,i) = constraints(j,i,4);
    end
end

% plotto la funzione obiettivo in funzione dei diametri in 3d
% --- Superficie 3D della funzione obiettivo ---
figure(4)
surf(D_vals, d_vals, obj);
xlabel('D'); ylabel('d'); zlabel('f(x)');

figure(5)
surf(D_vals, d_vals, obj_p);
xlabel('D'); ylabel('d'); zlabel('f(x)');
title('Superficie 3D della funzione obiettivo');

figure(6)
surf(D_vals, d_vals, constraint_stress);
xlabel('D'); ylabel('d'); zlabel('constraint_stress');
title('Superficie 3D del constraint di stress');



% --- Contour plot with constraints ---
figure(7);
contour_levels = [0 200 400 800 1600 3200 4000];
contour(D_vals, d_vals, obj, contour_levels, 'LineColor', 'k', 'ShowText', 'on','LineWidth', 1);
hold on;


h1 = contour(D_vals, d_vals, constraint_stress, [0 0], 'LineColor', 'g', 'LineWidth', 2);
h2 = contour(D_vals, d_vals, constraint_geom, [0 0], 'LineColor', 'b', 'LineWidth', 2);
h3 = contour(D_vals, d_vals, constraint_buckling, [0 0], 'LineColor', 'm', 'LineWidth', 2);
h4 = contour(D_vals, d_vals, constraint_def, [0 0], 'LineColor', 'r', 'LineWidth', 2);

% Labels
xlabel('D (mm)');
ylabel('d (mm)');
hold on
plot(x_opt(1),x_opt(2),'.', 'MarkerSize', 25);
hold on






%% ------------------ OPTIMIZATION WITH STEEPEST DESCENT ------------------


% 2. Initialization for the steepest descent process

% Initial design point
xq = x0;
plot(xq(1), xq(2), 'ks', 'MarkerFaceColor', 'k');
text(xq(1), xq(2), '  xq', 'FontSize', 10);

% Set the maximum number of optimization cycles
max_cycles = 20;
cycle = 0;
converged = false;

% Design variable bounds (for D and d)
D_lower = D_min; D_upper = D_max;
d_lower = d_min; d_upper = d_max;

scale = 0.001;  % Value adjustable based on your design context

% 3. Steepest Descent Loop
while ~converged && cycle < max_cycles
    cycle = cycle + 1;
    fprintf('\nCycle %d:\n', cycle);
    hi = 1e-8;
    % Evaluate the objective at the current design point
    fx = obj_with_pen(xq);
    % Perturb D
    fx_Dplus = obj_with_pen([xq(1)+hi,xq(2)]);
    % Perturb d
    fx_dplus = obj_with_pen([xq(1), xq(2)+hi]);
    % Compute finite differences (partial derivatives)
    dfdx1 = (fx_Dplus - fx) / hi;
    dfdx2 = (fx_dplus - fx) / hi;
    grad = [dfdx1, dfdx2];
    fprintf('Gradient: [%e, %e]\n', grad(1), grad(2));

    % --- Determine the steepest descent direction ---
    % In steepest descent, the search direction is -gradient.
    % Normalize the gradient to obtain a unit vector.
    norm_grad = norm(grad);
    if norm_grad == 0
        fprintf('Gradient is zero. Convergence achieved.\n');
        break;
    end
    % Unscaled descent direction:
    descent_dir = -grad / norm_grad;

    % --- Scale the search vector ---
    % Multiply the unit direction by the scaling factor to obtain a "reasonable" step length.
    sq = descent_dir * scale;
    fprintf('Scaled search direction: [%e, %e]\n', sq(1), sq(2));

    % --- Determine an appropriate upper bound for alpha ---
    % Ensure that the new design x = xq + alpha*sq stays within the bounds.
    if sq(1) > 0
        alpha_max_D = (D_upper - xq(1)) / sq(1);
    elseif sq(1) < 0
        alpha_max_D = (xq(1) - D_lower) / abs(sq(1));
    else
        alpha_max_D = Inf;
    end
    if sq(2) > 0
        alpha_max_d = (d_upper - xq(2))/ sq(2);
    elseif sq(2) < 0
        alpha_max_d = (xq(2) - d_lower)/ abs(sq(2));
    else
        alpha_max_d = Inf;
    end
    alpha_max_bound = min(alpha_max_D, alpha_max_d);
    fprintf('Alpha upper bound based on design limits: %e\n', alpha_max_bound);

    alpha_lower_bound = 0;

    % --- Perform a line search using fminbnd ---
    options = optimset('TolX', 1e-8);  % Use default options (you may adjust if desired)
    obj_alpha = @(alpha) obj_with_pen([xq(1) + alpha * sq(1), xq(2) + alpha * sq(2)]);  
    [alpha_opt, fval, exitflag] = fminbnd(obj_alpha, alpha_lower_bound, alpha_max_bound, options);
    fprintf('Optimal alpha: %e, f(alpha) = %e, exitflag = %d\n', alpha_opt, fval, exitflag);

    % --- Update the design point ---
    xnew = xq + alpha_opt * sq;
    fprintf('New design point: D = %e, d = %e\n', xnew(1), xnew(2));
    hold on
    % Visualize the step on the contour plot
    plot([xq(1), xnew(1)], [xq(2), xnew(2)], 'yo-', 'LineWidth', 2);
    hold on
    % --- Check for convergence ---
    if abs(fval - fx) < 1e-8
        converged = true;
        fprintf('Convergence reached: Change in f(alpha) < 1e-6\n');
    end

    % Update the current design point
    xq = xnew;

    pause(1);  % Optional pause to observe each cycle
end

fprintf('\nOptimization completed in %d cycles.\n', cycle);

%% ------------------ RESULTS ------------------
y2_opt = x_opt(1);
y3_opt = x_opt(2);
disp(['Objective function fmincon: ', num2str(f_opt)])
disp(['Objective function genetic algorithm: ', num2str(fval_ga)])
disp(['Objective function fmincon and genetic algorithm: ', num2str(f_opt_fga)])
disp(['Nodo 2: (', num2str(x2_init), ', ', num2str(y2_opt), ')']);
disp(['Nodo 3: (', num2str(x3_init), ', ', num2str(y3_opt), ')']);
disp(['Nodo 4: (', num2str(x4_init), ', ', num2str(y4_init), ')']);

% Plot optimized frame
figure(8);
nodes_opt = [0, 0; x2_init, y2_init; x3_init, y3_init; x4_init, y4_init];
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
title('Telaio Ottimizzato');
axis equal;
hold off;


%% ------------------ FUNCTIONS ------------------

function m = compute_mass(x, elements, pho, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init)
    D = x(1); d = x(2);
    S = (pi * D.^2 - pi * d.^2) / 4;
    nodes = [0, 0; x2_init, y2_init; x3_init, y3_init; x4_init, y4_init];
    lengths = zeros(size(elements,1),1);
    masses = zeros(size(elements,1),1);
    for i = 1:size(elements,1)
        n1 = elements(i,1);
        n2 = elements(i,2);
        p1 = nodes(n1, :);
        p2 = nodes(n2, :);
        lengths(i) = norm(p2 - p1);
        masses(i) = abs(lengths(i))*pho*S;
    end
    m = sum(masses);
end

function [c, ceq, sigma_massimo, forze_interne] = compute_stress(x, elements, forces, sigma_max, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init)
    D = x(1); d = x(2);
    nodes = [0, 0; x2_init, y2_init; x3_init, y3_init; x4_init, y4_init];
    dof = size(nodes,1) * 2;
    S = (pi * D.^2 - pi * d.^2) / 4;
    A = zeros(dof, size(elements,1));
    L_vect = zeros(size(elements,1),1); 

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
    T = A \ F;
    forze_interne = T;
    sigma = T ./ S;
    sigma_massimo = max(abs(sigma));
    D; d;

    stress_violation = max(abs(sigma)) - sigma_max;
    
    
    E = 210e9*10^(-6); 
    I = (pi / 64) * (D^4 - d^4); 
    K = 1;
    P_cr = (pi^2 * E * I) ./ (K^2 * max(L_vect).^2);
    buckling_violation = abs(max(T)) - P_cr;
    int_geometry_violation = d - D + 0.1;
    
    strain = T ./ (E * S);
    strain_violation = max(abs(strain)) - 1e-4;
    
    c = [stress_violation; buckling_violation; int_geometry_violation; strain_violation];
    ceq = [];
end

function [y_cm] = compute_mass_center(x, elements,pho, D, d)
    S = (pi * D.^2 - pi * d.^2) / 4;
    nodes = [0, 0; x(3), x(4); x(5), x(6); x(7), x(8)];
    lengths = zeros(size(elements,1),1);
    masses = zeros(size(elements,1),1);
    y_cm_s = zeros(size(elements,1),1);
    
    for i = 1:size(elements,1)
        n1 = elements(i,1);
        n2 = elements(i,2);
        p1 = nodes(n1, :);
        p2 = nodes(n2, :);
        lengths(i) = norm(p2 - p1);
        masses(i) = lengths(i)*pho*S;
        y_cm_s(i) = (p1(2) + p2(2))/2;
    end
    
    y_cm = sum(masses.*y_cm_s)/sum(masses); 
end



function min_buckling = compute_buckling_stress(x, elements, E, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init)
    D = x(1); d = x(2);
    I = (pi/64) * (D^4 - d^4);
    A = (pi/4) * (D^2 - d^2);
    K = 1; 

    nodes = [0, 0; x2_init, y2_init; x3_init, y3_init; x4_init, y4_init];
    buckling_stresses = zeros(size(elements,1), 1);

    for i = 1:size(elements,1)
        n1 = elements(i,1);
        n2 = elements(i,2);
        p1 = nodes(n1, :);
        p2 = nodes(n2, :);
        L = norm(p2 - p1);
        sigma_cr = (pi^2 * E * I) / ((K * L)^2 * A);
        buckling_stresses(i) = sigma_cr;
    end

    
    min_buckling = min(buckling_stresses);
end

function f = objective_with_buckling(x, elements, pho, E, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init, alpha)
    
    mass = compute_mass(x, elements, pho, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init);

    
    buckling_stress = compute_buckling_stress(x, elements, E, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init);
    
    
    f = mass; %- alpha * buckling_stress;
end

function f = objective_with_penalty(x, elements, pho, E, forces, ...
    x2_init, y2_init, x3_init, y3_init, x4_init, y4_init)

    D = x(1); d = x(2);
    S = (pi * D.^2 - pi * d.^2) / 4;
    nodes = [0, 0; x2_init, y2_init; x3_init, y3_init; x4_init, y4_init];
    dof = size(nodes,1) * 2;

    
    A = zeros(dof, size(elements,1));
    L_vect = zeros(size(elements,1),1);

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
    T = A \ F;

    
    strain = T ./ (E * S);
    strain_max = max(abs(strain));
    strain_limit = 1e-4;
    strain_violation = max(0, strain_max - strain_limit);

    lambda = 94;     
    alpha = 3;      


    
     penalty = lambda * (exp(alpha * strain_violation) - 1);
    
     
    lambda_geom = 1000;
    alpha_geom = 100;

    geom_violation = max(0, 0.1 - (D - d));
    penalty_geom = lambda_geom * (exp(alpha_geom * geom_violation) - 1);  
    


    % Calcolo massa (o altro obiettivo)
    mass = compute_mass(x, elements, pho, x2_init, y2_init, x3_init, y3_init, x4_init, y4_init);

    % Funzione obiettivo penalizzata
    f = mass + penalty + penalty_geom;
end
