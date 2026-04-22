%% ================================================================
%  Soft Robotic Gripper  –  LQR Animation
%  Fingers grip the object at equilibrium.
%  Disturbance pushes them open, LQR pulls them back to grasp.
%
%  LQR METHOD : Output-weighted  Q = C'C  (replaces Bryson's Rule)
%               R = ρ·I₂  (scalar torque penalty)
%
%  Keep fk_finger_4.m in the same folder.
%% ================================================================
clear; clc; close all;

%% ── 1. PARAMETERS ───────────────────────────────────────────────
L   = 0.15;
l   = L/3;
m   = 0.1;
Jl  = m*l^2;
b   = 0.011;
k   = 0.5;
kc  = 50;
d   = 0.12;    % gap between finger bases
w   = 0.04;    % object width

%% ── 2. FIND GRASPING EQUILIBRIUM ANGLE ──────────────────────────
base1 = -d/2;
base2 =  d/2;

reach_needed = (d - w) / 2 + 0.005;
tip_y = 3 * l * cos(asin(reach_needed / (3*l)));
th0   = asin(reach_needed / (3*l));
fprintf('Grasp equilibrium angle: th0 = %.2f deg\n', rad2deg(th0));

% Equilibrium state: all joints at th0, all velocities zero
x_eq = zeros(12,1);
x_eq(1) = th0;  x_eq(3) = th0;  x_eq(5) = th0;   % finger 1
x_eq(7) = th0;  x_eq(9) = th0;  x_eq(11)= th0;   % finger 2

%% ── 3. BUILD MATRICES ───────────────────────────────────────────
al  = -2*k/Jl;   be = -2*b/Jl;
ka  =  k/Jl;     mu =  b/Jl;
gc  =  kc*l^2/Jl;
a66 = -ka - gc;

A11 = [
  0,        1,       0,    0,    0,      0    ;
  al,       be,      ka,   mu,   0,      0    ;
  0,        0,       0,    1,    0,      0    ;
  ka,       mu,      al,   be,   ka,     mu   ;
  0,        0,       0,    0,    0,      1    ;
 -gc,       0,    ka-gc,   mu,  -ka-gc, -mu
];
A12 = zeros(6);
A12(6,1) = -gc;  A12(6,3) = -gc;  A12(6,5) = -gc;
A21 = zeros(6);
A21(6,1) = -gc;  A21(6,3) = -gc;  A21(6,5) = -gc;
A   = [A11,A12; A21,A11];

B      = zeros(12,2);
B(2,1) = 1/Jl;
B(8,2) = 1/Jl;

C        = zeros(7,12);
C(1,1)=1; C(2,3)=1; C(3,5)=1;
C(4,7)=1; C(5,9)=1; C(6,11)=1;
C(7,:) = kc*l*[1,0,1,0,1,0,1,0,1,0,1,0];
D = zeros(7,2);

%% ── 4. EQUILIBRIUM TORQUE ────────────────────────────────────────
% At equilibrium x_eq, x_dot=0:  0 = A*x_eq + B*u_eq
% => u_eq = -(B'B)^{-1} B' (A*x_eq)
u_eq = -(B'*B)\(B'*(A*x_eq));
fprintf('Equilibrium torques: tau1=%.4f  tau2=%.4f N.m\n', u_eq(1), u_eq(2));

%% ── 5. LQR ON PERTURBATION VARIABLES ────────────────────────────
%
%  CHANGE: Bryson's Rule  →  Output-Weighted LQR
%  ─────────────────────────────────────────────
%  OLD (Bryson's Rule):
%    Qblock = diag([1/theta_max^2, 1/dtheta_max^2, ...])
%    Q = blkdiag(Qblock, Qblock)
%    R = diag([1/tau_max^2, 1/tau_max^2])
%
%  NEW (Output-Weighted):
%    Q = C'*C          ← penalises all 7 sensor outputs including Fc
%    Q(ai,ai) *= 10    ← extra emphasis on angular position states
%    R = rho * eye(2)  ← single scalar torque penalty knob
%
%  Key advantage: row 7 of C encodes contact force Fc = kc·l·Σθ,
%  so the LQR cost now automatically penalises grip force deviations
%  — something Bryson's Rule completely missed.
%
fprintf('\n=== LQR Controller Design (Output-Weighted) ===\n');

% Angle (position) state indices — odd entries in state vector
ai = [1, 3, 5, 7, 9, 11];

% Build Q
Q = C' * C;                    % 12×12 output-weighted cost
Q(ai,ai) = Q(ai,ai) * 10;      % ×10 boost on angle states
                               %   → tighter positional regulation
                               %   → velocity transients tolerated more

% Build R
rho = 0.01;                    % torque penalty scalar
                               %   rho ↑  → smaller torques, slower settling
                               %   rho ↓  → larger torques, faster recovery
R   = rho * eye(2);            % equal penalty on both fingers

% Solve CARE:  A'P + PA - PBR^{-1}B'P + Q = 0
% Returns:     K = R^{-1} B' P
[K_lqr, ~, cl_poles] = lqr(A, B, Q, R);

% Diagnostics
lqr_slowest = max(real(cl_poles));
fprintf('  rho (torque penalty)     = %.4f\n',   rho);
fprintf('  LQR slowest Re(lambda)   = %.4f rad/s\n', lqr_slowest);
fprintf('  All Re(lambda) < 0?      %s\n', mat2str(all(real(cl_poles)<0)));
fprintf('  LQR gain K (2x12):\n');
disp(K_lqr);

A_cl = A - B*K_lqr;

%% ── 6. SIMULATE ─────────────────────────────────────────────────

% ✅ NOISE TOGGLE — change these to switch noise on/off
use_noise = true;          % true = noisy,  false = clean
noise_std = deg2rad(0.1);  % try: 0.1, 0.5, 1.0, 2.0 degrees

% Disturbance: push both fingers OPEN (away from object) by 8 deg
disturbance = deg2rad(-8);
x0_tilde = zeros(12,1);
x0_tilde(1) = disturbance;   % finger 1, link 1 pushed open
x0_tilde(3) = disturbance;   % finger 1, link 2
x0_tilde(7) = disturbance;   % finger 2, link 1
x0_tilde(9) = disturbance;   % finger 2, link 2

tspan  = 0:0.001:3;
N_sim  = numel(tspan);
dt     = 0.001;

if ~use_noise
    %% ── CLEAN simulation ─────────────────────────────────────────
    sys_cl  = ss(A_cl, B, C, D);
    [~, t_sim, xt] = initial(sys_cl, x0_tilde, tspan);

    % Control torques (perturbation)
    u_tilde = -(K_lqr * xt')';
    u_full  = u_tilde + repmat(u_eq', N_sim, 1);

else
    %% ── NOISY simulation ─────────────────────────────────────────
    xt      = zeros(N_sim, 12);
    u_tilde = zeros(N_sim, 2);
    xt(1,:) = x0_tilde';

    for i = 1:N_sim-1
        x_t = xt(i,:)';

        % Add measurement noise
        x_noisy = x_t + noise_std * randn(12,1);

        % LQR control on noisy measurement
        u_i          = -K_lqr * x_noisy;
        u_tilde(i,:) = u_i';

        % True (noise-free) dynamics
        x_dot        = A_cl * x_t + B * u_i;
        xt(i+1,:)    = x_t + dt * x_dot;
    end

    t_sim  = tspan';
    u_full = u_tilde + repmat(u_eq', N_sim, 1);
end

% Full state = equilibrium + perturbation
x_full = xt + repmat(x_eq', N_sim, 1);

% Extract angles
TH = x_full(:, [1,3,5,7,9,11]);

% Contact force (nonlinear calculation)
Fc_sim = zeros(N_sim, 1);
for i = 1:N_sim
    xT1 = l*sum(sin(TH(i,1:3)));
    xT2 = l*sum(sin(TH(i,4:6)));
    pen = xT1 + xT2 - (d-w);
    Fc_sim(i) = kc*max(0,pen);
end

fprintf('Simulation done. Initial grip angle: %.1f deg\n', rad2deg(th0));
fprintf('Disturbance: %.1f deg opening\n', rad2deg(abs(disturbance)));
if use_noise
    fprintf('Noise level: %.2f deg\n', rad2deg(noise_std));
else
    fprintf('Noise: OFF (clean simulation)\n');
end

%% ── 7. PRE-COMPUTE GEOMETRY ──────────────────────────────────────
N  = numel(t_sim);
P1 = zeros(2,4,N);
P2 = zeros(2,4,N);
for i = 1:N
    P1(:,:,i) = fk_finger_4(base1, TH(i,1:3)', l, +1);
    P2(:,:,i) = fk_finger_4(base2, TH(i,4:6)', l, -1);
end

%% ── 8. ANIMATION FIGURE ─────────────────────────────────────────
fps    = 25;
stride = max(1, round(1/(fps*0.001)));
frames = 1:stride:N;

fig = figure(2);
clf(fig);
set(fig, 'Name',     'Gripper LQR Animation  –  Output-Weighted', ...
         'Color',    'white', ...
         'Position', [120 80 700 740]);

ax = axes('Parent',fig,'Position',[0.09 0.09 0.87 0.87]);
hold(ax,'on');
axis(ax,'equal');
axis(ax,[-0.115  0.115  -0.020  tip_y + 0.07]);
set(ax,'FontSize',11);
box(ax,'on'); grid(ax,'on');
set(ax,'GridAlpha',0.20);
xlabel(ax,'x  (m)','FontSize',12);
ylabel(ax,'y  (m)','FontSize',12);

%% ── 9. STATIC SCENE ──────────────────────────────────────────────

% Base floor
fill(ax,[-0.115 0.115 0.115 -0.115],[-0.020 -0.020 0 0], ...
     [0.80 0.80 0.83],'EdgeColor','none');
for xh = -0.11:0.013:0.11
    plot(ax,[xh xh+0.009],[-0.018 -0.002], ...
         'Color',[0.62 0.62 0.66],'LineWidth',0.9);
end

% Object
tip_x1 = base1 + 3*l*sin(th0);
tip_x2 = base2 - 3*l*sin(th0);
w_draw = tip_x2 - tip_x1;
obj_y0 = tip_y - 0.04;
obj_ht = 0.08;
fill(ax,[tip_x1, tip_x2, tip_x2, tip_x1], ...
         [obj_y0, obj_y0, obj_y0+obj_ht, obj_y0+obj_ht], ...
     [0.93 0.78 0.48], ...
     'EdgeColor',[0.58 0.38 0.12],'LineWidth',2.2,'FaceAlpha',0.92);
for yg = obj_y0+0.012 : 0.018 : obj_y0+obj_ht-0.005
    plot(ax,[tip_x1+0.003, tip_x2-0.003],[yg yg], ...
         'Color',[0.75 0.55 0.25],'LineWidth',0.7,'LineStyle','-');
end
text(ax, 0, obj_y0+obj_ht/2, 'OBJECT', ...
     'HorizontalAlignment','center','FontSize',10, ...
     'FontWeight','bold','Color',[0.48 0.28 0.05]);

% Equilibrium finger lines
eq_pts1 = fk_finger_4(base1, [th0;th0;th0], l, +1);
eq_pts2 = fk_finger_4(base2, [th0;th0;th0], l, -1);
plot(ax, eq_pts1(1,:), eq_pts1(2,:), '--', ...
     'Color',[0.55 0.65 0.90],'LineWidth',1.2);
plot(ax, eq_pts2(1,:), eq_pts2(2,:), '--', ...
     'Color',[0.90 0.60 0.55],'LineWidth',1.2);
text(ax, eq_pts1(1,4)-0.005, eq_pts1(2,4)+0.008, 'target', ...
     'FontSize',7,'Color',[0.40 0.50 0.80],'HorizontalAlignment','right');
text(ax, eq_pts2(1,4)+0.005, eq_pts2(2,4)+0.008, 'target', ...
     'FontSize',7,'Color',[0.80 0.35 0.30]);

%% ── 10. FINGER HANDLES ───────────────────────────────────────────
c1 = [0.18 0.45 0.85];
c2 = [0.85 0.22 0.15];

p1_0 = P1(:,:,1);
p2_0 = P2(:,:,1);

h_f1 = plot(ax, p1_0(1,:), p1_0(2,:), '-o', ...
            'Color',c1,'LineWidth',6, ...
            'MarkerSize',11,'MarkerFaceColor',c1,'MarkerEdgeColor','w');
h_f2 = plot(ax, p2_0(1,:), p2_0(2,:), '-o', ...
            'Color',c2,'LineWidth',6, ...
            'MarkerSize',11,'MarkerFaceColor',c2,'MarkerEdgeColor','w');

h_tip1 = plot(ax, p1_0(1,4), p1_0(2,4), 'o', ...
              'MarkerSize',15,'MarkerFaceColor',c1,'MarkerEdgeColor','w','LineWidth',2);
h_tip2 = plot(ax, p2_0(1,4), p2_0(2,4), 'o', ...
              'MarkerSize',15,'MarkerFaceColor',c2,'MarkerEdgeColor','w','LineWidth',2);

plot(ax, base1, 0,'s','MarkerSize',16,'MarkerFaceColor',c1,'MarkerEdgeColor','w','LineWidth',2);
plot(ax, base2, 0,'s','MarkerSize',16,'MarkerFaceColor',c2,'MarkerEdgeColor','w','LineWidth',2);

legend(ax,[h_f1,h_f2],{'Finger 1 (left)','Finger 2 (right)'}, ...
       'Location','southwest','FontSize',10);

%% ── 11. TEXT OVERLAYS ────────────────────────────────────────────
h_title = title(ax, sprintf('Gripper LQR (Output-Weighted)  –  Equilibrium = %.1f°', ...
                rad2deg(th0)), 'FontSize',13,'FontWeight','bold');

h_status = text(ax, 0, -0.013, ...
    'Disturbance applied  –  fingers pushed open', ...
    'HorizontalAlignment','center','FontSize',11, ...
    'FontWeight','bold','Color',[0.80 0.35 0.05]);

h_info = text(ax, -0.112, 0.175, '', ...
              'FontSize',8,'Color',[0.20 0.20 0.35], ...
              'FontName','Courier','VerticalAlignment','top');

%% ── 12. ANIMATE ─────────────────────────────────────────────────
settled = false;

for fi = 1:numel(frames)
    i = frames(fi);

    set(h_f1,   'XData',P1(1,:,i),'YData',P1(2,:,i));
    set(h_f2,   'XData',P2(1,:,i),'YData',P2(2,:,i));
    set(h_tip1, 'XData',P1(1,4,i),'YData',P1(2,4,i));
    set(h_tip2, 'XData',P2(1,4,i),'YData',P2(2,4,i));

    set(h_title,'String', sprintf( ...
        'Gripper LQR (Output-Weighted)  –  t = %.3f s  |  F_c = %.3f N', ...
        t_sim(i), Fc_sim(i)));

    pert_deg = max(abs(rad2deg(xt(i,[1,3,5,7,9,11]))));

    if pert_deg < 0.3 && ~settled
        settled = true;
    end

    if settled
        set(h_status,'String', ...
            sprintf('Grip restored  ✓  |  F_c = %.3f N', Fc_sim(i)), ...
            'Color',[0.05 0.45 0.15]);
    else
        set(h_status,'String', ...
            sprintf('LQR Active  –  closing grip  |  error = %.2f°', pert_deg), ...
            'Color',[0.05 0.50 0.15]);
    end

    ang = rad2deg(TH(i,:));
    set(h_info,'String', sprintf( ...
        ['F1: θ1=%+.1f°  θ2=%+.1f°  θ3=%+.1f°\n' ...
         'F2: θ1=%+.1f°  θ2=%+.1f°  θ3=%+.1f°\n' ...
         'τ1=%+.3f  τ2=%+.3f N·m'], ...
        ang(1),ang(2),ang(3),ang(4),ang(5),ang(6), ...
        u_full(i,1),u_full(i,2)));

    drawnow;
    pause(1/fps);
end

set(h_title,'String','LQR Stabilisation Complete  –  Object Gripped (Output-Weighted)');
fprintf('\nAnimation done. Figure 2 stays open.\n');

%% ── 13. RESPONSE PLOTS (Figure 3) ───────────────────────────────
figure(3);
set(gcf,'Name','LQR Response (Output-Weighted)','Color','white','Position',[840 80 580 720]);

subplot(3,1,1);
hold on;
lbl = {'\theta_{11}','\theta_{12}','\theta_{13}', ...
       '\theta_{21}','\theta_{22}','\theta_{23}'};
clrs = [c1;c1*0.6+[0 0.3 0];c1*0.4+[0.3 0.3 0]; ...
        c2;c2*0.6+[0 0.2 0];c2*0.4+[0.2 0.2 0]];
for j = 1:6
    ls = '-'; if j>3, ls='--'; end
    plot(t_sim, rad2deg(TH(:,j)), ls,'Color',clrs(j,:), ...
         'LineWidth',1.5,'DisplayName',lbl{j});
end
yline(rad2deg(th0),'k:','LineWidth',1.2,'DisplayName','Equilibrium');
ylabel('Angle (deg)','FontSize',11);
title('Joint Angles  (solid=F1, dashed=F2)','FontSize',12,'FontWeight','bold');
legend('NumColumns',4,'FontSize',7,'Location','southeast');
grid on; xlim([0 3]);

subplot(3,1,2);
plot(t_sim,u_full(:,1),'-','Color',c1,'LineWidth',1.8,'DisplayName','\tau_1');
hold on;
plot(t_sim,u_full(:,2),'-','Color',c2,'LineWidth',1.8,'DisplayName','\tau_2');
yline(u_eq(1),'--','Color',c1*0.6,'LineWidth',1.0,'DisplayName','\tau_{1,eq}');
yline(u_eq(2),'--','Color',c2*0.6,'LineWidth',1.0,'DisplayName','\tau_{2,eq}');
ylabel('Torque (N·m)','FontSize',11);
title('Control Torques','FontSize',12,'FontWeight','bold');
legend('FontSize',9); grid on; xlim([0 3]);

subplot(3,1,3);
plot(t_sim,Fc_sim,'-','Color',[0.80 0.15 0.10],'LineWidth',1.8);
yline(Fc_sim(end),'k:','LineWidth',1.2);
ylabel('Force (N)','FontSize',11);
xlabel('Time (s)','FontSize',11);
title('Contact Force F_c','FontSize',12,'FontWeight','bold');
grid on; xlim([0 3]);

sgtitle('Output-Weighted LQR – Grip Recovery Response', ...
        'FontSize',13,'FontWeight','bold');

%% ── 14. SENSOR NOISE COMPARISON PLOT ────────────────────────────
% Runs all 4 noise levels using the output-weighted K_lqr

noise_levels = [0, 0.5, 1.0, 2.0];   % degrees (0 = clean)
line_colors  = [
    0.10 0.10 0.80;   % clean  = blue
    0.20 0.70 0.20;   % 0.5°   = green
    0.90 0.60 0.10;   % 1.0°   = orange
    0.85 0.10 0.10;   % 2.0°   = red
];

figure(7);
set(gcf,'Name','Noise Comparison (Output-Weighted LQR)','Color','white', ...
        'Position',[100 100 900 650]);

sp1 = subplot(3,1,1); hold on; grid on;
title('Joint Angle θ_{11} Recovery','FontSize',11,'FontWeight','bold');
ylabel('Perturbation (deg)'); xlabel('Time (s)');

sp2 = subplot(3,1,2); hold on; grid on;
title('Control Torque τ_1','FontSize',11,'FontWeight','bold');
ylabel('Torque (N·m)'); xlabel('Time (s)');

sp3 = subplot(3,1,3); hold on; grid on;
title('Contact Force Fc','FontSize',11,'FontWeight','bold');
ylabel('Force (N)'); xlabel('Time (s)');

for ni = 1:numel(noise_levels)
    ns  = deg2rad(noise_levels(ni));
    clr = line_colors(ni,:);

    if ns == 0
        % Clean — use ss/initial
        sys_cl_t = ss(A_cl, B, C, D);
        [~, ~, xt_n] = initial(sys_cl_t, x0_tilde, tspan);
        ut_n = -(K_lqr * xt_n')';
        lbl  = 'Clean (\sigma=0°)';
        lw   = 2.5;
    else
        % Noisy — Euler integration with measurement noise
        xt_n = zeros(N_sim,12);
        ut_n = zeros(N_sim,2);
        xt_n(1,:) = x0_tilde';

        for i = 1:N_sim-1
            x_t  = xt_n(i,:)';
            x_ny = x_t + ns*randn(12,1);
            u_i  = -K_lqr * x_ny;
            ut_n(i,:) = u_i';
            xt_n(i+1,:) = x_t + dt*(A_cl*x_t + B*u_i);
        end
        lbl = sprintf('\\sigma = %.1f°', noise_levels(ni));
        lw  = 1.5;
    end

    % Contact force for this noise run
    x_full_n = xt_n + repmat(x_eq',N_sim,1);
    TH_n     = x_full_n(:,[1,3,5,7,9,11]);
    Fc_n     = zeros(N_sim,1);
    for i = 1:N_sim
        xT1 = l*sum(sin(TH_n(i,1:3)));
        xT2 = l*sum(sin(TH_n(i,4:6)));
        Fc_n(i) = kc*max(0, xT1+xT2-(d-w));
    end

    plot(sp1, tspan, rad2deg(xt_n(:,1)), ...
         'Color',clr,'LineWidth',lw,'DisplayName',lbl);
    plot(sp2, tspan, ut_n(:,1)+u_eq(1), ...
         'Color',clr,'LineWidth',lw,'DisplayName',lbl);
    plot(sp3, tspan, Fc_n, ...
         'Color',clr,'LineWidth',lw,'DisplayName',lbl);
end

legend(sp1,'Location','southeast','FontSize',9);
legend(sp2,'Location','northeast','FontSize',9);
legend(sp3,'Location','southeast','FontSize',9);
xlim(sp1,[0 3]); xlim(sp2,[0 3]); xlim(sp3,[0 3]);

sgtitle('Effect of Sensor Noise on LQR Grip Recovery (Output-Weighted)', ...
        'FontSize',13,'FontWeight','bold');

fprintf('\nDone. Active controller: Output-Weighted LQR (rho = %.4f)\n', rho);
