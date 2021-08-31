
clc; clear all; 


% -------------------------------------------------------------------------
% Simulation Parameters
% -------------------------------------------------------------------------

% System params
global J Jinv
% Inertia Matrix
J = [17.5,   -0.8,    0.3; ...
     -0.8,   14.9,    0.4; ...
      0.3,    0.4,   20.8];
Jinv = inv(J);

wmax = .005;        % Magnitude of maximum possible disturbance
W = wmax*[-1, 1].*[1;1;1];

% Backup controller paramaters
global umax kp kd
umax = 0.5;         % Magnitude of maximum possible control input
kp = 0.6;           % Porportional gain
kd = 2.25;          % Derivative gain

dt = .01;               % Simulation Time Step
T_max = 15;             % Simulation Time
timevec = 0:dt:T_max;


% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------

% initial state (no transformation)
x0 = [sqrt(3)/2, sqrt(1)/2, 0, 0, .1, .1, .1]'; 
traj = zeros(7, size(timevec, 2));
traj(:, 1) = x0;
[zl, zu] =  z_getter([x0;x0]);
thing_sauce = acosd(zl);


% initial state (with transformation)
x0T = x0; 
x0T(5:7) = J*x0(5:7); 

% embedding state (no transformation)
E_state = zeros(14, size(timevec, 2));
E_state(:, 1) = [x0; x0];
z_holder = zeros(2, size(timevec, 2));
[zl, zu] =  z_getter(E_state(:, 1));
z_holder(:, 1) = [acosd(zu); acosd(zl)];

% embedding state (with transformation)
E_stateT = zeros(14, size(timevec, 2));
E_stateT(:, 1) = [x0T; x0T];
z_holderT = zeros(2, size(timevec, 2));
[zl, zu] =  z_getter(E_stateT(:, 1));
z_holderT(:, 1) = [acosd(zu); acosd(zl)];


% -------------------------------------------------------------------------
% Simulation
% -------------------------------------------------------------------------

for i = 1:size(timevec, 2) - 1
    xnow  = E_state(1:7, i);
    xhnow = E_state(8:14, i);
    E_state(:, i + 1) = E_state(:, i) + ...
                        dt*E_init( xnow, xhnow, W(:, 1), W(:, 2), J, Jinv); 

                    
    [zl, zu] =  z_getter(E_state(:, i+1));
    z_holder(:, i+1) = [acosd(zu); acosd(zl)];
    
    xnowT  = E_stateT(1:7, i);
    xhnowT = E_stateT(8:14, i);
    E_stateT(:, i + 1) = E_stateT(:, i) + ...
                        dt*E_tran( xnowT, xhnowT, W(:, 1), W(:, 2), J, Jinv); 

    [zl, zu] =  z_getter(E_stateT(:, i+1));
    z_holderT(:, i+1) = [acosd(zu); acosd(zl)];
    
    traj(:, i+1) = traj(:, i) + dt*F(traj(:, i), J, Jinv);
    [zl, zu] =  z_getter([traj(:, i+1);traj(:, i+1)]);
    thing_sauce(:, end+1) = acosd(zl);
end






% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

LW_backup_traj = 2;
backup_traj_color = [ 0, 199, 199]./255;
backup_traj_color2 = [1, 0, 0];


figure(1); clf; 
hold on; grid on;

scatter(0, z_holder(1, 1), 80, 'filled', 'b', 'DisplayName', 'Initial State')

plot(timevec, z_holder(1, :), 'Color', backup_traj_color, 'LineWidth', LW_backup_traj, 'HandleVisibility', 'off');
plot(timevec, z_holder(2, :), 'Color', backup_traj_color, 'LineWidth', LW_backup_traj, 'DisplayName', 'RSO w/o Transform');
reach_set = patch([timevec,fliplr(timevec)], ...
                  [z_holder(1, :), fliplr(z_holder(2, :))], ...
                  'c', 'facealpha', 0.1, 'edgealpha', 0, 'HandleVisibility', 'off');


plot(timevec, z_holderT(1, :), 'Color', backup_traj_color2, 'LineWidth', LW_backup_traj, 'HandleVisibility', 'off');
plot(timevec, z_holderT(2, :), 'Color', backup_traj_color2, 'LineWidth', LW_backup_traj, 'DisplayName', 'RSO w/ Transform');
reach_set2 = patch([timevec,fliplr(timevec)], ...
                  [z_holderT(1, :), fliplr(z_holderT(2, :))], ...
                  'r', 'facealpha', 0.1, 'edgealpha', 0, 'HandleVisibility', 'off');

plot(timevec, thing_sauce(1, :), 'Color', 'k', 'LineWidth', LW_backup_traj, 'DisplayName', 'Deterministic Trajectory');

scatter(0, z_holder(1, 1), 80, 'filled', 'b', 'HandleVisibility', 'off')

xlabel('time [s]','Interpreter','latex')
ylabel('Opening Angle $\beta$ [deg]', 'Interpreter', 'latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1, 16, 0, 100])
xticks([0, 5, 10, 15])

Leg = legend('Location', 'southwest');
%set(Leg,'visible','off')
drawnow



% -------------------------------------------------------------------------
% Functoins
% -------------------------------------------------------------------------

% System Dynamics
function out = F(x, J, Jinv)
    global kp kd umax
    out(1:4, 1) = [1/2*(-x(2)*x(5) - x(3)*x(6)  -x(4)*x(7)); ...
                   1/2*(x(1)*x(5) - x(4)*x(6) + x(3)*x(7)); ...
                   1/2*(x(4)*x(5) + x(1)*x(6)  - x(2)*x(7)); ...
                   1/2*(-x(3)*x(5) + x(2)*x(6) + x(1)*x(7))];
       
    eta = [2*(x(3)*x(4) + x(1)*x(2)); ...
           2*(x(1)*x(3) - x(2)*x(4)); ...
           0];
       
    out(5:7, 1) = -Jinv*(cross(x(5:7), J*x(5:7))) + ...
        Jinv*umax*tanh((1/umax)*( cross(x(5:7), J*x(5:7)) - kp*J*eta - kd*J*x(5:7) ));
end

% Embedding function for initial system dynamics
function out = E_init(x, xh, w, wh, J, Jinv)
    global umax kp kd
    out = [decomp_init( x, xh,  w, wh, J, Jinv, umax, kp, kd); ...
           decomp_init(xh,  x, wh,  w, J, Jinv, umax, kp, kd)];
end

% Embedding function for system dynamics with transformed first subsystem
function out = E_tran(x, xh, w, wh, J, Jinv)
    global umax kp kd
    out = [decomp_tran( x, xh,  w, wh, J, Jinv, umax, kp, kd); ...
           decomp_tran(xh,  x, wh,  w, J, Jinv, umax, kp, kd)];
end

% Saturation Function
function out = sigma(in)
    out = umax*tanh((1/umax)*in);
end


function [zl, zu] = z_getter(EE)
    zl = min([ 1 - 2*EE(2)^2 - 2*EE(3)^2, ...
        1 - 2*EE(9)^2 - 2*EE(3)^2, ...
        1 - 2*EE(2)^2 - 2*EE(10)^2, ...
        1 - 2*EE(9)^2 - 2*EE(10)^2]);
    zu = 1;
    
    if EE(2) <= 0 && 0 <= EE(9)
    else
        zu = zu + max([- 2*EE(2)^2, - 2*EE(9)^2]);
    end
    if EE(3) <= 0 && 0 <= EE(10)
    else
        zu = zu + max([- 2*EE(3)^2, - 2*EE(10)^2]);
    end
end

