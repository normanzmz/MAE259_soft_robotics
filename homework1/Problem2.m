clc; clear; close all;
% Problem 2 Implicit method
N = 21; % Number of spheres 
dt = 0.01; % time step size 
t_tot = 50; % total time
N_step = round (t_tot / dt); % Number of time steps

[all_y, all_v] = Problem2_function(N, dt, t_tot);

figure(1)
for i = 2 : N_step
    plot(all_y(1:2:end,i),all_y(2:2:end,i),'ro-');
    axis equal
    drawnow
end

all_mid_y = all_y(N+1,:);
all_mid_v = all_v(N+1,:);

figure (2);
timearr = (1:N_step) * dt;
plot (timearr,all_mid_y, 'k-');
xlabel ('time')
ylabel ('veritical position of mid-node')

figure (3);
plot (timearr,all_mid_v, 'k-');
xlabel ('time')
ylabel ('velocity of mid-node')

%% Question 3.1: terminal velocity vs. num of nodes 
Node_number = 5;
dt = 1e-2;
all_yn = cell(Node_number,1);
all_vn = cell(Node_number,1);

for N = 5:4:21
    [all_yn{N}, all_vn{N}] = Problem2_function(N, dt, t_tot);
end
Nodes_N = 6:6:Node_number*6;
vel_ter = zeros(Node_number,1);

for i = 5:4:21
    [wid,len] = size(all_yn{i});
    vel_ter(i) = all_vn{i}(wid/2+1,len);
end
vel_ter = nonzeros(vel_ter);

figure (4);
axis equal;
plot(Nodes_N,vel_ter, 'ro-');
xlabel('Number of Nodes');
ylabel('Terminal Velocity');

%% Question 3.2: terminal velocity vs. time step size 
step_number = 5;
N = 21;
all_ys = cell(step_number,1);
all_vs = cell(step_number,1);

for i = 1:step_number
    dt = 0.01/(10*i);
    time_N(i) = dt;
    [all_ys{i}, all_vs{i}] = Problem2_function(N, dt, t_tot);
end

vel_ter = zeros(step_number,1);

for i = 1:step_number
    [wid,len] = size(all_ys{i});
    vel_ter(i) = all_vs{i}(wid/2+1,len);
end

figure (5);
axis equal;
plot(time_N,vel_ter, 'ro-');
xlabel('Step size');
ylabel('Terminal Velocity');
