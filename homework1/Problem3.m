clc; clear; close all;
% Problem 3 Implicit method 
% Physical variables
N = 50; % Number of spheres 
l = 1; % beam length
R = 0.013; % beam outer radius 
r = 0.011; % beam inner radius
rho = 2700; % Density of aluminum 
dl = l / (N - 1); % length of each segment 

p = 2000; % external force 
d = 0.75; % distance from left side 

% find the closest node 
P = zeros(2*N,1);
P(round(N*3/4)*2) = -p;

E = 70e9; % Youngs Modulus
A = pi*(R^2-r^2); % Beam cross section 
g = 9.81; % gravity 

dt = 0.01; % time step size 
t_tot = 1; % total time

m = pi*(R^2 - r^2)*l*rho/(N-1); % mass of each node  
I = pi/4*(R^4 - r^4); % Moment of inerita 

EA = E*A; % Stretching stiffness 
EI = E*I; % Bending stiffness 

% Nodes initial configuration 
nodes_I = zeros(N,2);
for i = 1:N 
    nodes_I(i,1) = (i-1) * dl;
end 

% Mass Matrix 
M = zeros(2*N,2*N);
for i = 1 : 2*N
    M(i,i) = m;
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
ymax = zeros(N_step, 1);

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

    % Weight 
    f = f - P;

    % Update 
    q(3:(2*N-1)) = q(3:(2*N-1)) - J(3:(2*N-1),3:(2*N-1)) \ f(3:(2*N-1));
    err = sum(abs(f(3:(2*N-1))));
    end

    % Update
    u = (q - q0) / dt;
    q0 = q;
    
    % plot of sphere position
    figure(1);
    plot(q(1:2:end),q(2:2:end),'ro-');
    axis([0 1 -0.4 0])
    drawnow

    % Store 
    all_mid_y(i) = q(N+1);
    all_mid_v(i) = u(N+1);
    ymax(i) = min(q);
end 

%% Plot of time vs. max displacement
figure (2);
timearr = (1:N_step) * dt;
plot(timearr, ymax, 'k-')
xlabel('Time')
ylabel('Maximum Displacement')
grid on

% Euler beam theory 
c = min(d, 1-d);
ymax_E = -(p*c*(l^2 - c^2)^1.5)/(9*sqrt(3)*EI*l)

