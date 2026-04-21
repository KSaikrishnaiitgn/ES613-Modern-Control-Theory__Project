%% ================================================================
%  Soft Robotic Gripper  –  LQR Animation  (AUTO DISTURBANCE SWEEP)
%
%  Automatically cycles through disturbance values in ONE figure:
%    Left  : animation panel
%    Right : response plot (all runs overlaid, colour-coded)
%
%  LQR METHOD : Output-weighted  Q = C'C
%               R = rho * I2
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
d   = 0.12;
w   = 0.04;

%% ── 2. GRASPING EQUILIBRIUM ANGLE ───────────────────────────────
base1 = -d/2;
base2 =  d/2;

reach_needed = (d - w) / 2 + 0.005;
tip_y = 3 * l * cos(asin(reach_needed / (3*l)));
th0   = asin(reach_needed / (3*l));
fprintf('Grasp equilibrium angle: th0 = %.2f deg\n', rad2deg(th0));

x_eq = zeros(12,1);
x_eq(1)=th0; x_eq(3)=th0; x_eq(5)=th0;
x_eq(7)=th0; x_eq(9)=th0; x_eq(11)=th0;

%% ── 3. BUILD MATRICES ───────────────────────────────────────────
al  = -2*k/Jl;   be = -2*b/Jl;
ka  =  k/Jl;     mu =  b/Jl;
gc  =  kc*l^2/Jl;

A11 = [
  0,  1,  0,  0,  0,  0;
  al, be, ka, mu, 0,  0;
  0,  0,  0,  1,  0,  0;
  ka, mu, al, be, ka, mu;
  0,  0,  0,  0,  0,  1;
 -gc, 0,  ka-gc, mu, -ka-gc, -mu
];

A12 = zeros(6);
A12(6,1) = -gc; A12(6,3) = -gc; A12(6,5) = -gc;

A21 = zeros(6);
A21(6,1) = -gc; A21(6,3) = -gc; A21(6,5) = -gc;

A = [A11, A12; A21, A11];

B      = zeros(12,2);
B(2,1) = 1/Jl;
B(8,2) = 1/Jl;

C        = zeros(7,12);
C(1,1)=1; C(2,3)=1; C(3,5)=1;
C(4,7)=1; C(5,9)=1; C(6,11)=1;
C(7,:)   = kc*l*[1,0,1,0,1,0,1,0,1,0,1,0];
D        = zeros(7,2);

%% ── 4. EQUILIBRIUM TORQUE ────────────────────────────────────────
u_eq = -(B'*B)\(B'*(A*x_eq));

%% ── 5. LQR ──────────────────────────────────────────────────────
ai = [1, 3, 5, 7, 9, 11];
Q  = C' * C;
Q(ai,ai) = Q(ai,ai) * 10;
rho = 0.01;
R   = rho * eye(2);

[K_lqr, ~, cl_poles] = lqr(A, B, Q, R);
fprintf('All poles stable: %d\n', all(real(cl_poles)<0));

A_cl   = A - B*K_lqr;
sys_cl = ss(A_cl, B, C, D);
tspan  = 0:0.001:3;

%% ── 6. DISTURBANCE SWEEP SETTINGS ──────────────────────────────
%  ▶ Edit this list to add/remove disturbance levels:
dist_list_deg = [2, 5, 8, 15, 20];

%  Colour for each disturbance (matched to severity bands)
% dist_colors = [
%     0.20 0.70 0.20;   % 2  deg – green
%     0.10 0.50 0.90;   % 5  deg – blue
%     0.90 0.70 0.10;   % 8  deg – yellow
%     0.90 0.40 0.10;   % 15 deg – orange
%     0.85 0.10 0.10;   % 20 deg – red
% ];

n_dist = numel(dist_list_deg);
dist_colors = colormap(lines(n_dist));
fps    = 25;

%% ── 7. PRE-SIMULATE ALL DISTURBANCES ────────────────────────────
fprintf('\nPre-simulating %d disturbance cases...\n', n_dist);
ALL_TH   = cell(n_dist,1);
ALL_xt   = cell(n_dist,1);
ALL_Fc   = cell(n_dist,1);
settle_t = nan(n_dist,1);

for di = 1:n_dist
    dist = deg2rad(dist_list_deg(di));

    x0 = zeros(12,1);
    x0(1)=-dist; x0(3)=-dist; x0(7)=-dist; x0(9)=-dist;

    [~, t_sim, xt] = initial(sys_cl, x0, tspan);
    TH_d = xt(:,[1,3,5,7,9,11]) + th0;

    Fc_d = zeros(numel(t_sim),1);
    for i = 1:numel(t_sim)
        xT1 = l*sum(sin(TH_d(i,1:3)));
        xT2 = l*sum(sin(TH_d(i,4:6)));
        pen = xT1 + xT2 - (d-w);
        Fc_d(i) = kc*max(0,pen);
    end

    err_norm    = max(abs(xt(:,[1,3,5,7,9,11])), [], 2);
    settled_idx = find(err_norm < deg2rad(0.5), 1, 'first');
    if ~isempty(settled_idx)
        settle_t(di) = t_sim(settled_idx);
        fprintf('  %2d deg → settled at %.3f s\n', dist_list_deg(di), settle_t(di));
    else
        fprintf('  %2d deg → DID NOT SETTLE\n', dist_list_deg(di));
    end

    ALL_TH{di} = TH_d;
    ALL_xt{di} = xt;
    ALL_Fc{di} = Fc_d;
end

%% ── 8. BUILD FIGURE (single figure, two panels) ─────────────────
fig = figure(10); clf(fig);
set(fig, 'Name',     'Gripper LQR – Auto Disturbance Sweep', ...
         'Color',    'white', ...
         'Position', [60 60 1200 680]);

%% Left panel – animation
ax_anim = subplot('Position', [0.03 0.10 0.50 0.82]);
hold(ax_anim,'on'); axis equal;
axis(ax_anim, [-0.115  0.115  -0.020  tip_y+0.09]);
set(ax_anim,'FontSize',10); box on; grid on;
set(ax_anim,'GridAlpha',0.18);
xlabel(ax_anim,'x  (m)','FontSize',11);
ylabel(ax_anim,'y  (m)','FontSize',11);

%% Right panel – response plot
ax_plot = subplot('Position', [0.58 0.10 0.40 0.82]);
hold(ax_plot,'on');
xlabel(ax_plot,'Time (s)','FontSize',11);
ylabel(ax_plot,'θ₁ perturbation  (deg)','FontSize',11);
title(ax_plot, 'Recovery Response – All Disturbances', ...
      'FontSize',12,'FontWeight','bold');
grid(ax_plot,'on');
xlim(ax_plot,[0 3]); ylim(ax_plot,[-25 5]);
yline(ax_plot, 0, 'k:', 'LineWidth',1.5);

%% ── 9. STATIC SCENE (drawn once, never redrawn) ─────────────────
% Floor
fill(ax_anim, [-0.115 0.115 0.115 -0.115],[-0.020 -0.020 0 0], ...
     [0.80 0.80 0.83],'EdgeColor','none');
for xh = -0.11:0.013:0.11
    plot(ax_anim,[xh xh+0.009],[-0.018 -0.002], ...
         'Color',[0.62 0.62 0.66],'LineWidth',0.9);
end

% Object
tip_x1 = base1 + 3*l*sin(th0);
tip_x2 = base2 - 3*l*sin(th0);
obj_y0  = tip_y - 0.04;
obj_ht  = 0.08;
fill(ax_anim, [tip_x1,tip_x2,tip_x2,tip_x1], ...
     [obj_y0,obj_y0,obj_y0+obj_ht,obj_y0+obj_ht], ...
     [0.93 0.78 0.48], 'EdgeColor',[0.58 0.38 0.12], ...
     'LineWidth',2.2,'FaceAlpha',0.92);
text(ax_anim, 0, obj_y0+obj_ht/2, 'OBJECT', ...
     'HorizontalAlignment','center','FontSize',10, ...
     'FontWeight','bold','Color',[0.48 0.28 0.05]);

% Equilibrium guide lines
eq_pts1 = fk_finger_4(base1,[th0;th0;th0],l,+1);
eq_pts2 = fk_finger_4(base2,[th0;th0;th0],l,-1);
plot(ax_anim, eq_pts1(1,:),eq_pts1(2,:),'--', ...
     'Color',[0.55 0.65 0.90],'LineWidth',1.2);
plot(ax_anim, eq_pts2(1,:),eq_pts2(2,:),'--', ...
     'Color',[0.90 0.60 0.55],'LineWidth',1.2);
text(ax_anim, eq_pts1(1,4)-0.005, eq_pts1(2,4)+0.008, 'target', ...
     'FontSize',7,'Color',[0.40 0.50 0.80],'HorizontalAlignment','right');
text(ax_anim, eq_pts2(1,4)+0.005, eq_pts2(2,4)+0.008, 'target', ...
     'FontSize',7,'Color',[0.80 0.35 0.30]);

% Base mounts
plot(ax_anim, base1, 0,'s','MarkerSize',16, ...
     'MarkerFaceColor',[0.3 0.3 0.8],'MarkerEdgeColor','w','LineWidth',2);
plot(ax_anim, base2, 0,'s','MarkerSize',16, ...
     'MarkerFaceColor',[0.8 0.3 0.3],'MarkerEdgeColor','w','LineWidth',2);

%% ── 10. DYNAMIC HANDLES (reused each disturbance run) ────────────
h_f1   = plot(ax_anim,[0 0],[0 0],'-o','Color',[0 0 0],'LineWidth',6, ...
              'MarkerSize',10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','w');
h_f2   = plot(ax_anim,[0 0],[0 0],'-o','Color',[0 0 0],'LineWidth',6, ...
              'MarkerSize',10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor','w');
h_tip1 = plot(ax_anim, 0, 0,'o','MarkerSize',14,'MarkerFaceColor',[0 0 0], ...
              'MarkerEdgeColor','w','LineWidth',2);
h_tip2 = plot(ax_anim, 0, 0,'o','MarkerSize',14,'MarkerFaceColor',[0 0 0], ...
              'MarkerEdgeColor','w','LineWidth',2);

h_title  = title(ax_anim, 'Initialising...', ...
                 'FontSize',12,'FontWeight','bold');
h_status = text(ax_anim, 0, -0.013, '', ...
                'HorizontalAlignment','center','FontSize',12, ...
                'FontWeight','bold','Color',[0 0 0]);
h_info   = text(ax_anim, -0.112, tip_y+0.075, '', ...
                'FontSize',9,'Color',[0.2 0.2 0.35], ...
                'VerticalAlignment','top');

% Live line on right panel (updated each frame)
h_live = plot(ax_plot, NaN, NaN, '-', 'LineWidth',2.5, 'Color',[0 0 0]);

%% ── 11. Legend handles collected during sweep ────────────────────
leg_handles = gobjects(n_dist,1);   % filled inside main loop

%% ── 12. MAIN LOOP – CYCLE THROUGH DISTURBANCES ──────────────────
for di = 1:n_dist
    clr      = dist_colors(di,:);
    dist_deg = dist_list_deg(di);
    TH_d     = ALL_TH{di};
    xt_d     = ALL_xt{di};
    Fc_d     = ALL_Fc{di};
    N        = size(TH_d,1);

    % Update color of dynamic elements
    set(h_f1,   'Color',clr,'MarkerFaceColor',clr);
    set(h_f2,   'Color',clr,'MarkerFaceColor',clr);
    set(h_tip1, 'Color',clr,'MarkerFaceColor',clr);
    set(h_tip2, 'Color',clr,'MarkerFaceColor',clr);
    set(h_live, 'Color',clr,'XData',NaN,'YData',NaN);
    set(h_status,'Color',clr);
    set(h_title, 'Color',clr);

    %% Phase A – Flash disturbance
    P1_0 = fk_finger_4(base1, TH_d(1,1:3)', l, +1);
    P2_0 = fk_finger_4(base2, TH_d(1,4:6)', l, -1);
    set(h_f1,  'XData',P1_0(1,:),'YData',P1_0(2,:));
    set(h_f2,  'XData',P2_0(1,:),'YData',P2_0(2,:));
    set(h_tip1,'XData',P1_0(1,4),'YData',P1_0(2,4));
    set(h_tip2,'XData',P2_0(1,4),'YData',P2_0(2,4));
    set(h_title, 'String', sprintf('Case %d/%d  –  %d° Disturbance Applied ⚡', ...
                  di, n_dist, dist_deg));
    set(h_status,'String', sprintf('⚡ %d° disturbance!', dist_deg));
    set(h_info,  'String', sprintf('Case %d of %d\nDisturbance: %d°\nrho = %.4f', ...
                  di, n_dist, dist_deg, rho));
    drawnow;

    for flash = 1:3
        set(h_status,'Visible','off'); drawnow; pause(0.10);
        set(h_status,'Visible','on');  drawnow; pause(0.10);
    end
    pause(0.3);

    %% Phase B – Animate recovery
    stride = max(1, round(1/(fps*0.001)));
    frames = 1:stride:N;
    settled_shown = false;

    for fi = 1:numel(frames)
        i = frames(fi);

        P1_i = fk_finger_4(base1, TH_d(i,1:3)', l, +1);
        P2_i = fk_finger_4(base2, TH_d(i,4:6)', l, -1);
        set(h_f1,   'XData',P1_i(1,:),'YData',P1_i(2,:));
        set(h_f2,   'XData',P2_i(1,:),'YData',P2_i(2,:));
        set(h_tip1, 'XData',P1_i(1,4),'YData',P1_i(2,4));
        set(h_tip2, 'XData',P2_i(1,4),'YData',P2_i(2,4));

        % Update live line
        set(h_live, 'XData', tspan(1:i), ...
                    'YData', rad2deg(xt_d(1:i, 1)));

        pert = rad2deg(max(abs(xt_d(i, [1,3,5,7,9,11]))));

        set(h_title,'String', sprintf( ...
            'Case %d/%d  |  %d° Disturbance  |  t = %.2f s  |  error = %.1f°', ...
            di, n_dist, dist_deg, tspan(i), pert));

        set(h_info,'String', sprintf( ...
            'Case %d of %d\nDisturbance: %d°\nFc = %.4f N\nrho = %.4f', ...
            di, n_dist, dist_deg, Fc_d(i), rho));

        if pert > 1
            set(h_status,'String', ...
                sprintf('LQR ACTIVE – closing grip | error = %.1f°', pert));
        elseif ~settled_shown
            settled_shown = true;
            if ~isnan(settle_t(di))
                set(h_status,'String', ...
                    sprintf('✓ Grip Restored in %.2f s  |  Fc = %.3f N', ...
                    settle_t(di), Fc_d(i)));
            else
                set(h_status,'String','✗ LQR FAILED – could not recover!', ...
                    'Color',[0.85 0.10 0.10]);
            end
        end

        drawnow;
        pause(1/fps);
    end

    %% Freeze settled state and draw permanent trace on right panel
    % Replace live line with a permanent coloured trace
    leg_handles(di) = plot(ax_plot, tspan(1:N), rad2deg(xt_d(:,1)), '-', ...
         'Color', clr, 'LineWidth', 2.0, ...
         'DisplayName', sprintf('%d°', dist_deg));

    % Reset live line to invisible until next case
    set(h_live,'XData',NaN,'YData',NaN);

    if ~isnan(settle_t(di))
        set(h_title,'String', sprintf( ...
            '✓  Case %d/%d  –  %d° Recovered in %.2f s', ...
            di, n_dist, dist_deg, settle_t(di)), 'Color',[0.05 0.55 0.20]);
    else
        set(h_title,'String', sprintf( ...
            '✗  Case %d/%d  –  %d° FAILED', di, n_dist, dist_deg), ...
            'Color',[0.85 0.10 0.10]);
    end

    pause(0.8);   % brief pause between cases
end

%% ── 13. FINAL STATE ──────────────────────────────────────────────
set(h_title,'String', ...
    sprintf('✓  Auto Sweep Complete  –  %d Disturbance Cases (rho = %.4f)', ...
    n_dist, rho), 'Color',[0.05 0.40 0.15]);
set(h_status,'String', 'All cases finished. Edit dist_list_deg to change sweep.', ...
    'Color',[0.25 0.25 0.35]);
set(h_info,'String', sprintf( ...
    'Sweep: %s deg\nMethod: Output-Weighted (Q=C''C)\nrho = %.4f', ...
    mat2str(dist_list_deg), rho));

title(ax_plot, 'Recovery Response – All Disturbances (Overlaid)', ...
      'FontSize',12,'FontWeight','bold');
legend(ax_plot, leg_handles, 'Location','southeast','FontSize',9, ...
       'Title','Disturbance');

fprintf('\nSweep complete. All %d cases shown.\n', n_dist);
fprintf('To change sweep values, edit dist_list_deg on line ~130.\n');
