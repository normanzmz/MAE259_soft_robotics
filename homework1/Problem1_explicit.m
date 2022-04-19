clc; clear; close all;
% Problem 1 Explicit method 
%% Physical variables
% Number of spheres 
N = 3;

% radius of sphere
R1 = 0.005; 
R2 = 0.025;
R3 = 0.005;

l = 0.1; % beam length
R0 = 0.001; % beam radius 
dl = l / (N - 1); % length of each segment 

% Density of fluid 
p_metal = 7000;
p_fluid = 1000;
p_tot = p_metal - p_fluid;

v_fluid = 1000; % viscosity of fluid Pa-s
E = 1e9; % Youngs Modulus
g = 9.81; % gravity 

dt = 1e-5; % time step size 
t_tot = 10; % total time

EA = E*pi*R0^2; % Stretching stiffness 
EI = (E*pi*R0^4)/4; % Bending stiffness 

% Nodes initial configuration 
nodes_I = zeros(N,2);
for i = 1:N 
    nodes_I(i,1) = (i-1) * dl;
end 

% Mass Matrix 
M = zeros(2*N, 2*N);
M(1,1) = 4/3 * pi * R1^3 * p_metal;
M(2,2) = 4/3 * pi * R1^3 * p_metal;
M(3,3) = 4/3 * pi * R2^3 * p_metal;
M(4,4) = 4/3 * pi * R2^3 * p_metal;
M(5,5) = 4/3 * pi * R3^3 * p_metal;
M(6,6) = 4/3 * pi * R3^3 * p_metal;

% Damping Matrix 
C = zeros(2*N,2*N);
C(1,1) = 6 * pi * v_fluid * R1;
C(2,2) = 6 * pi * v_fluid * R1;
C(3,3) = 6 * pi * v_fluid * R2;
C(4,4) = 6 * pi * v_fluid * R2;
C(5,5) = 6 * pi * v_fluid * R3;
C(6,6) = 6 * pi * v_fluid * R3;

% Weight vector 
W = zeros(2*N,1);
W(2) = - 4/3 * pi * R1^3 * p_tot * g;
W(4) = - 4/3 * pi * R2^3 * p_tot * g;
W(6) = - 4/3 * pi * R3^3 * p_tot * g;

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
all_mid_y = zeros(6,N_step); % y position of R2
all_mid_v = zeros(6,N_step); % velocity of R2 

all_mid_y(:,1) = q;
all_mid_v(:,1) = u;

tol = EI / l^2 * 1e-3; % tolerance 

% time marching scheme
for i = 2 : N_step
    
    fprintf('Time = %f\n',(i-1) * dt);
    f = zeros(6,1);

    % Elastic forces 
    % Linear Spring 1 between node 1 and 2
    xk = q(1);
    yk = q(2);
    xkp1 = q(3);
    ykp1 = q(4); 
    dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
    f(1:4) = f(1:4) - dF;

    % Linear Spring 2 between node 2 and 3 
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6); 
    dF = gradEs(xk, yk, xkp1, ykp1, dl, EA);
    f(3:6) = f(3:6) - dF;

    % Bending spring 
    xkm1 = q(1);
    ykm1 = q(2);
    xk = q(3);
    yk = q(4);
    xkp1 = q(5);
    ykp1 = q(6);
    curvature0 = 0;
    dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, dl, EI);
    f(1:6) = f(1:6) - dF;

    % viscous force 
    v_fluid = C * u;
    f(1:6) = f(1:6)-v_fluid;

    % Weight 
    f = f + W;
    
    % Intermidia step 1
    step1 = dt*inv(M) * f + u;

    % Update 
    q = step1*dt + q;
    
    % Update
    u = (q - q0) / dt;
    q0 = q;
    
    % Store 
    all_mid_y(:,i) = q;
    all_mid_v(:,i) = u;
end 

figure(1);
for i = 1:10000:length(all_mid_y)
    plot(all_mid_y(1:2:end,i),all_mid_y(2:2:end,i),'ro-');
    axis equal 
    drawnow
end

figure (2);
timearr = (1:N_step) * dt;
plot (timearr,all_mid_v(N+1,:), 'k-');
xlabel ('time')
ylabel ('velocity of mid-node')

