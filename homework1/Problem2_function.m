function [all_y, all_v] = Problem2_function(N,dt,t_tot)
% Number of spheres
l = 0.1; % beam length
R0 = 0.001; % beam radius
dl = l / (N - 1); % length of each segment

% radius of sphere
Rn = dl/10;
Rm = 0.025;

% Density of fluid
p_metal = 7000;
p_fluid = 1000;
p_tot = p_metal - p_fluid;

v_fluid = 1000; % viscosity of fluid Pa-s
E = 1e9; % Youngs Modulus
g = 9.81; % gravity

EA = E*pi*R0^2; % Stretching stiffness
EI = (E*pi*R0^4)/4; % Bending stiffness

% Nodes initial configuration
nodes_I = zeros(N,2);
for i = 1:N
    nodes_I(i,1) = (i-1) * dl;
end

% Mass & Damping Matrix
M = zeros(2*N,2*N);
C = zeros(2*N,2*N);
for i = 1 : 2*N
    if i == N || i == N+1
        M(i,i) = 4/3 * pi * Rm^3 * p_metal;
        C(i,i) = 6 * pi * v_fluid * Rm;
    else
        M(i,i) = 4/3 * pi * Rn^3 * p_metal;
        C(i,i) = 6 * pi * v_fluid * Rn;
    end
end

% Weight vector
W = zeros(2*N,1);
for i = 1: N
    if i == (N+1)/2
        W(2*i) = - 4/3 * pi * Rm^3 * p_tot * g;
    else
        W(2*i) = - 4/3 * pi * Rn^3 * p_tot * g;
    end
end

% Initial DOF vector
q0 = zeros(2*N,1);
for i = 1:N
    q0 (2*i-1) = nodes_I(i,1); % x coordinate
    q0 (2*i) = nodes_I(i,2); % y coordinate
end

% New position and velocity
q = q0;
u = (q - q0) / dt;

% Number of time steps
N_step = round (t_tot / dt);
all_mid_y = zeros(N_step, 1); % y position of R2
all_mid_v = zeros(N_step, 1); % velocity of R2

all_mid_y(1) = q(N+1);
all_mid_v(1) = u(N+1);

tol = EI / l^2 * 1e-3; % tolerance

% time marching scheme
for i = 2 : N_step

    fprintf('Time = %f\n',(i-1) * dt);

    q = q0; % Guess

    err = 10 * tol;
    % Newton Raphson
    while err > tol

        % Inertia
        f = M/dt * ((q - q0) / dt - u);
        J = M/ dt^2;

        % Elastic forces
        % Linear Spring 1 between node 1 and 2
        for k = 1: N-1
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
            dJ = hessEs(xk, yk, xkp1, ykp1, dl, EA);
            f(2*k-1:2*k+2) = f(2*k-1:2*k+2) + dF;
            J(2*k-1:2*k+2, 2*k-1:2*k+2) = ...
                J(2*k-1:2*k+2, 2*k-1:2*k+2) + dJ;
        end

        % Bending spring
        for k = 2:N-1
            xkm1 = q(2*k-3);
            ykm1 = q(2*k-2);
            xk = q(2*k-1);
            yk = q(2*k);
            xkp1 = q(2*k+1);
            ykp1 = q(2*k+2);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, dl, EI);
            dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, dl, EI);
            f(2*k-3:2*k+2) = f(2*k-3:2*k+2) + dF;
            J(2*k-3:2*k+2, 2*k-3:2*k+2) = ...
                J(2*k-3:2*k+2, 2*k-3:2*k+2) + dJ;
        end

        % viscous force
        f = f + C * (q - q0) / dt;
        J = J + C / dt;

        % Weight
        f = f - W;

        % Update
        q = q - J \ f;
        err = sum(abs(f));
    end

    % Update
    u = (q - q0) / dt;
    q0 = q;

    % Store
    %     all_mid_y(i) = q(N+1);
    %     all_mid_v(i) = u(N+1);
    all_y(:,i) = q;
    all_v(:,i) = u;
end
end


