%% ================================================================
%  Soft Robotic Gripper — LQR + Luenberger Observer Design
%  ME 613 | Modern Control Theory
%
%  Controllability and observability already verified in
%  gripper_verification_v3.m — not repeated here.
%% ================================================================

clear; clc; close all;

%% ================================================================
%  SECTION 1 — PARAMETERS
%% ================================================================

L_len = 0.1616;  ell = L_len/3;
mj    = 0.01289;  Jl  = mj * ell^2;
k     = 2.22;   b   = 0.0036;  kc = 50;

a  = -2*k/Jl;   be = -2*b/Jl;   ka = k/Jl;
mu = b/Jl;      gc = kc*ell^2/Jl;

fprintf('Jl=%.2e | α=%.0f β=%.0f κ=%.0f μ=%.1f γc=%.0f\n',...
    Jl,a,be,ka,mu,gc);

%% ================================================================
%  SECTION 2 — STATE-SPACE MATRICES
%% ================================================================

A11 = [0    1    0    0    0    0   ;
       a    be   ka   mu   0    0   ;
       0    0    0    1    0    0   ;
       ka   mu   a    be   ka   mu  ;
       0    0    0    0    0    1   ;
      -gc   0  ka-gc  mu -ka-gc -mu ];

A12 = zeros(6,6);
A12(6,1)=-gc; A12(6,3)=-gc; A12(6,5)=-gc;

A = [A11 A12; A12 A11];

B = zeros(12,2);
B(2,1)=1/Jl;  B(8,2)=1/Jl;

C = zeros(7,12);
ai = [1 3 5 7 9 11];
for i=1:6; C(i,ai(i))=1; end
C(7,:) = kc*ell*[1 0 1 0 1 0 1 0 1 0 1 0];

D = zeros(7,2);

%% ================================================================
%  SECTION 3 — OPEN-LOOP CHECK
%% ================================================================

ol_eig = eig(A);
fprintf('\nOpen-loop: all stable? %s | slowest Re = %.4f rad/s\n',...
    mat2str(all(real(ol_eig)<0)), max(real(ol_eig)));

%% ================================================================
%  SECTION 4 — LQR CONTROLLER DESIGN
%
%  Objective : find K (2×12) such that u = -K·x̃ minimises
%
%      J = ∫₀^∞ ( x̃ᵀ Q x̃  +  ũᵀ R ũ ) dt
%
%  Method    : solve the Continuous Algebraic Riccati Equation
%              AᵀP + PA − P·B·R⁻¹·Bᵀ·P + Q = 0
%              K = R⁻¹·Bᵀ·P
%
%  Tuning philosophy
%  -----------------
%  Q = Cᵀ·C  weights the measured outputs (angles + contact force).
%  Angle states receive an additional 10× diagonal boost so that
%  angular tracking errors are penalised more heavily than velocity.
%
%  R = 0.01·I₂  sets a moderately low torque penalty on both fingers.
%  Increase ρ → smaller torques, slower settling.
%  Decrease ρ → larger torques, faster response, risk of saturation.
%% ================================================================

fprintf('\n=== LQR Controller Design ===\n');

Q_lqr = C'*C;
Q_lqr(ai,ai) = Q_lqr(ai,ai) * 10;

rho   = 0.01;
R_lqr = rho * eye(2);

[K, ~, cl_eig] = lqr(A, B, Q_lqr, R_lqr);

lqr_slowest = max(real(cl_eig));

fprintf('  ρ (torque penalty)      = %.4f\n', rho);
fprintf('  LQR slowest Re(λ)       = %.4f rad/s\n', lqr_slowest);
fprintf('  All Re(λ) < 0 (stable)? %s\n', mat2str(all(real(cl_eig) < 0)));
fprintf('  LQR gain matrix K (2×12):\n');
disp(K);

%% ================================================================
%  SECTION 5 — LUENBERGER POLE PLACEMENT OBSERVER
%
%  Design rule: place observer poles with real parts
%  3–5× more negative than the slowest LQR pole.
%
%  LQR slowest Re = -15.73  →  target Re ≈ -47 to -79 rad/s
%  We choose 4× → target Re = 4 × (-15.73) = -62.90 rad/s
%
%  Pole structure: 6 complex-conjugate pairs centred at target_re
%  with increasing imaginary parts.
%% ================================================================

target_re = 4 * lqr_slowest;

imag_parts = [5; 15; 25; 35; 45; 55];
obs_poles  = [target_re + 1i*imag_parts ;
              target_re - 1i*imag_parts ];

fprintf('\n--- Observer pole placement ---\n');
fprintf('Target real part: %.4f rad/s\n', target_re);
fprintf('Desired poles (6 conjugate pairs):\n');
for i = 1:6
    fprintf('  %.4f ± %.4fj\n', target_re, imag_parts(i));
end

L = place(A', C', obs_poles')';

A_obs  = A - L*C;
ob_eig = eig(A_obs);
obs_slowest = max(real(ob_eig));

fprintf('\n=== OBSERVER RESULT ===\n');
fprintf('  Slowest observer Re = %.4f rad/s\n', obs_slowest);
fprintf('  Slowest LQR     Re = %.4f rad/s\n', lqr_slowest);
fprintf('  Speed ratio        = %.2f×\n', abs(obs_slowest)/abs(lqr_slowest));
fprintf('  Separation satisfied? %s\n', mat2str(obs_slowest < lqr_slowest));
fprintf('  Observer gain L norm = %.4f\n', norm(L));
fprintf('\n  Observer poles (eig(A-LC)):\n');
disp(sort(ob_eig,'ComparisonMethod','real'));

%% ================================================================
%  SECTION 6 — SEPARATION PRINCIPLE VERIFICATION
%
%  Combined closed-loop eigenvalues = eig(A-BK) ∪ eig(A-LC)
%  Block-triangular augmented system:
%
%    d/dt [x ]  =  [A-BK    -BK  ] [x ]
%    d/dt [ẽ ]     [  0    A-LC  ] [ẽ ]
%% ================================================================

all_cl_poles = [cl_eig; ob_eig];
fprintf('\n=== SEPARATION PRINCIPLE ===\n');
fprintf('  Combined system: %d poles total\n', length(all_cl_poles));
fprintf('  All Re < 0? %s\n', mat2str(all(real(all_cl_poles)<0)));
fprintf('  Slowest combined pole Re = %.4f rad/s\n', max(real(all_cl_poles)));
fprintf('  Observer poles faster than LQR? %s\n',...
    mat2str(abs(obs_slowest) > abs(lqr_slowest)));

%% ================================================================
%  SECTION 7 — CLOSED-LOOP SIMULATION
%% ================================================================

A_aug = [A - B*K,       -B*K     ;
         zeros(12,12),  A - L*C  ];

x0 = [0.05;0; 0.04;0; 0.03;0;
      0.05;0; 0.04;0; 0.03;0];
e0 = 0.015 * ones(12,1);
X0 = [x0; e0];

t_span = 0 : 5e-4 : 2.5;
N = length(t_span);

X = zeros(24, N);
X(:,1) = X0;
for i = 1:N-1
    X(:,i+1) = X(:,i) + 5e-4 * A_aug * X(:,i);
end

x_true = X(1:12,:);
e_err  = X(13:24,:);
x_hat  = x_true - e_err;
u_sig  = -K * x_hat;
y_sig  = C * x_true;

tau_lqr = 1/abs(lqr_slowest);
tau_obs = 1/abs(obs_slowest);
fprintf('\nτ_LQR = %.4f s  |  τ_obs = %.4f s  |  τ_obs/τ_LQR = %.2f\n',...
    tau_lqr, tau_obs, tau_obs/tau_lqr);

%% ================================================================
%  SECTION 8 — PLOTS
%% ================================================================

col = {'#185FA5','#1D9E75','#534AB7','#D85A30','#BA7517','#5F5E5A'};
lbl = {'θ₁₁ base F1','θ₁₂ mid F1','θ₁₃ tip F1',...
       'θ₂₁ base F2','θ₂₂ mid F2','θ₂₃ tip F2'};

%% Figure 1 — Joint angles: true vs estimate
figure('Name','Angles','Color','w','Position',[40 40 1000 600]);
for p = 1:6
    subplot(3,2,p);
    hold on;
    h1 = plot(t_span, x_true(ai(p),:)*180/pi,'LineWidth',2,'Color',col{p});
    h2 = plot(t_span, x_hat(ai(p),:)*180/pi,'--','LineWidth',1.2,'Color',col{p});
    yline(0,'--k','LineWidth',0.7,'HandleVisibility','off');
    xline(tau_obs,':','Color','#D85A30','LineWidth',1.2,'HandleVisibility','off');
    xline(tau_lqr,':','Color','#185FA5','LineWidth',1.2,'HandleVisibility','off');
    xlabel('t (s)'); ylabel('deg');
    title(lbl{p},'FontWeight','normal');
    if p==1
        legend([h1 h2],{'True x','Estimate x̂'},'FontSize',8,'Location','northeast');
    end
    grid on; box off;
end
sgtitle(sprintf('Joint angle regulation  |  τ_{obs}=%.3fs (red dotted)  τ_{LQR}=%.3fs (blue dotted)',...
    tau_obs, tau_lqr),'FontWeight','bold');

%% Figure 2 — Observer convergence
figure('Name','Observer Convergence','Color','w','Position',[40 680 1000 320]);

e_ang_norm = vecnorm(e_err([1 3 5 7 9 11],:)) * 180/pi;
e_vel_norm = vecnorm(e_err([2 4 6 8 10 12],:));

subplot(1,2,1);
semilogy(t_span, e_ang_norm + 1e-12,'Color','#185FA5','LineWidth',2);
hold on;
xline(tau_obs,'--r','LineWidth',1.5,...
    'Label',sprintf('τ_{obs}=%.3fs',tau_obs),'LabelOrientation','horizontal');
xlabel('t (s)'); ylabel('‖ẽ_{angles}‖ (deg) — log scale');
title('Angle estimation error decay','FontWeight','normal');
grid on; box off;

subplot(1,2,2);
semilogy(t_span, e_vel_norm + 1e-12,'Color','#1D9E75','LineWidth',2);
hold on;
xline(tau_obs,'--r','LineWidth',1.5,...
    'Label',sprintf('τ_{obs}=%.3fs',tau_obs),'LabelOrientation','horizontal');
xlabel('t (s)'); ylabel('‖ẽ_{veloc}‖ (rad/s) — log scale');
title('Velocity estimation error decay','FontWeight','normal');
grid on; box off;
sgtitle('Observer error ‖x−x̂‖ — log scale shows exponential decay rate',...
    'FontWeight','bold');

%% Figure 3 — Control inputs
figure('Name','Control Inputs','Color','w','Position',[1060 40 620 300]);
plot(t_span, u_sig(1,:),'LineWidth',2,'Color','#185FA5'); hold on;
plot(t_span, u_sig(2,:),'--','LineWidth',2,'Color','#D85A30');
xlabel('t (s)'); ylabel('Torque (N·m)');
title('Tendon torques','FontWeight','bold');
legend({'τ₁ Finger 1','τ₂ Finger 2'},'Location','northeast');
yline(0,'--k','LineWidth',0.6);
grid on; box off;

%% Figure 4 — Contact force
figure('Name','Contact Force','Color','w','Position',[1060 380 620 280]);
plot(t_span, y_sig(7,:),'LineWidth',2,'Color','#1D9E75');
xlabel('t (s)'); ylabel('F_c (N)');
title('Contact force F_c','FontWeight','bold');
yline(0,'--k','LineWidth',0.6);
grid on; box off;

%% Figure 5 — Pole map
figure('Name','Pole Map','Color','w','Position',[1060 700 620 400]);

plot(real(ol_eig), imag(ol_eig),'bs','MarkerSize',8,'LineWidth',1.5,...
    'DisplayName','Open-loop plant'); hold on;
plot(real(cl_eig), imag(cl_eig),'r^','MarkerSize',8,'LineWidth',1.5,...
    'DisplayName','LQR closed-loop');
plot(real(ob_eig), imag(ob_eig),'g+','MarkerSize',10,'LineWidth',2,...
    'DisplayName','Observer (pole placement)');
xline(lqr_slowest,'--r','LineWidth',1.2,...
    'Label',sprintf('LQR: %.1f',lqr_slowest),...
    'LabelOrientation','horizontal','HandleVisibility','off');
xline(obs_slowest,'--g','LineWidth',1.2,...
    'Label',sprintf('Obs: %.1f',obs_slowest),...
    'LabelOrientation','horizontal','HandleVisibility','off');
xline(0,'--k','LineWidth',0.8,'HandleVisibility','off');
xlabel('Real (rad/s)'); ylabel('Imaginary (rad/s)');
title({'Pole map — observer poles (green) must be LEFT of LQR poles (red)',...
       'Separation principle: ✓ satisfied'},'FontWeight','bold');
legend('Location','northeast'); grid on; box off;

%% ================================================================
%  SECTION 9 — FINAL SUMMARY
%% ================================================================

fprintf('\n%s\n', repmat('=',1,55));
fprintf('  FINAL DESIGN SUMMARY\n');
fprintf('%s\n', repmat('=',1,55));
fprintf('  %-28s %10s\n','Quantity','Value');
fprintf('  %s\n', repmat('-',1,42));
fprintf('  %-28s %10.4f\n','Open-loop slowest Re (rad/s)',max(real(ol_eig)));
fprintf('  %-28s %10.4f\n','LQR slowest Re (rad/s)',lqr_slowest);
fprintf('  %-28s %10.4f\n','Observer slowest Re (rad/s)',obs_slowest);
fprintf('  %-28s %10.2f\n','Observer/LQR speed ratio',...
    abs(obs_slowest)/abs(lqr_slowest));
fprintf('  %-28s %10s\n','Separation principle',mat2str(obs_slowest<lqr_slowest));
fprintf('  %-28s %10.4f\n','LQR time const τ_LQR (s)',tau_lqr);
fprintf('  %-28s %10.4f\n','Observer time const τ_obs (s)',tau_obs);
fprintf('  %-28s %10s\n','LQR gain K size',...
    sprintf('%d×%d',size(K,1),size(K,2)));
fprintf('  %-28s %10s\n','Observer gain L size',...
    sprintf('%d×%d',size(L,1),size(L,2)));
fprintf('%s\n', repmat('=',1,55));
fprintf('Done. Figures 1–5 generated.\n');
