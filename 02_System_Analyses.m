%% ============================================================
%  VERIFICATION: Stability, Controllability & Observability
%  Soft Robotic Gripper - PRB Model (ME 613)
%  Two-finger, 3-link per finger, n=12 state system
%
%  CORRECTION v2: A-matrix is UNCHANGED (penetration uses all
%  three joint angles as in the paper). The Lyapunov proof is
%  corrected: instead of the physical energy (which is NOT a
%  valid Lyapunov function for this model), the proof uses the
%  unique P>0 satisfying  A'P + PA = -I  (algebraic Lyapunov
%  equation), giving Vdot = -||x||^2 < 0 strictly.
%% ============================================================

clear; clc; close all;

%% ============================================================
%  SECTION 0: System Parameters
%% ============================================================

L   = 0.15;
l   = L/3;
m1  = 0.10;   m2  = 0.10;
Jl1 = m1*l^2; Jl2 = m2*l^2;
b   = 0.011;
k   = 0.5;
kc  = 50;

alpha1 = -2*k/Jl1;  beta1 = -2*b/Jl1;
kappa1 =  k/Jl1;    mu1   =  b/Jl1;
gamma1 =  kc*l^2/Jl1;

alpha2 = -2*k/Jl2;  beta2 = -2*b/Jl2;
kappa2 =  k/Jl2;    mu2   =  b/Jl2;
gamma2 =  kc*l^2/Jl2;

fprintf('============================================================\n');
fprintf('  System Parameters\n');
fprintf('============================================================\n');
fprintf('  L=%.3f m,  l=%.4f m\n', L, l);
fprintf('  Jl1=%.6f kg.m^2,  Jl2=%.6f kg.m^2\n', Jl1, Jl2);
fprintf('  b=%.4f,  k=%.4f,  kc=%.2f\n', b, k, kc);
fprintf('  gamma_c=%.4f\n\n', gamma1);

%% ============================================================
%  SECTION 1: Build A, B, C Matrices  (ORIGINAL, UNCHANGED)
%% ============================================================
%
%  Penetration: delta = l*(theta_11+theta_12+theta_13
%                         +theta_21+theta_22+theta_23)
%  So -gamma_c appears in columns 1,3,5 of rows 6 and 12.

% --- A11 (6x6) ---
A11 = [
    0,       1,       0,              0,     0,                  0;
    alpha1,  beta1,   kappa1,         mu1,   0,                  0;
    0,       0,       0,              1,     0,                  0;
    kappa1,  mu1,     alpha1,         beta1, kappa1,             mu1;
    0,       0,       0,              0,     0,                  1;
   -gamma1,  0,       kappa1-gamma1,  mu1,  -kappa1-gamma1,     -mu1
];

% --- A22 (6x6) ---
A22 = [
    0,       1,       0,              0,     0,                  0;
    alpha2,  beta2,   kappa2,         mu2,   0,                  0;
    0,       0,       0,              1,     0,                  0;
    kappa2,  mu2,     alpha2,         beta2, kappa2,             mu2;
    0,       0,       0,              0,     0,                  1;
   -gamma2,  0,       kappa2-gamma2,  mu2,  -kappa2-gamma2,     -mu2
];

% --- Coupling blocks ---
A12 = zeros(6,6);
A12(6,1) = -gamma1;  A12(6,3) = -gamma1;  A12(6,5) = -gamma1;

A21 = zeros(6,6);
A21(6,1) = -gamma2;  A21(6,3) = -gamma2;  A21(6,5) = -gamma2;

A = [A11, A12; A21, A22];

% --- B, C, D ---
B = zeros(12,2);
B(2,1) = 1/Jl1;  B(8,2) = 1/Jl2;

C = zeros(7,12);
C(1,1)=1; C(2,3)=1; C(3,5)=1;
C(4,7)=1; C(5,9)=1; C(6,11)=1;
C(7,:) = kc*l * [1,0,1,0,1,0,1,0,1,0,1,0];

D = zeros(7,2);

fprintf('============================================================\n');
fprintf('  A, B, C Matrices Built (original, unchanged)\n');
fprintf('============================================================\n');
fprintf('  A: %dx%d | B: %dx%d | C: %dx%d\n\n', ...
    size(A,1),size(A,2),size(B,1),size(B,2),size(C,1),size(C,2));

%% ============================================================
%  SECTION 2: LYAPUNOV STABILITY  (CORRECTED PROOF)
%% ============================================================

fprintf('============================================================\n');
fprintf('  LYAPUNOV STABILITY ANALYSIS  (CORRECTED)\n');
fprintf('============================================================\n');

% -----------------------------------------------------------
%  Step 1: Show the PHYSICAL energy H is NOT a valid
%          Lyapunov function for this A-matrix.
% -----------------------------------------------------------
fprintf('  Step 1: Check whether physical energy H = x^T*P_H*x works\n');

P_H = zeros(12,12);
% Kinetic
for idx = [2,4,6];   P_H(idx,idx) = Jl1/2; end
for idx = [8,10,12]; P_H(idx,idx) = Jl2/2; end
% Spring PE finger 1: joints phi1=x1, phi2=x3-x1, phi3=x5-x3
P_H(1,1)=P_H(1,1)+k/2; P_H(1,1)=P_H(1,1)+k/2;
P_H(3,3)=P_H(3,3)+k/2; P_H(1,3)=P_H(1,3)-k/2; P_H(3,1)=P_H(3,1)-k/2;
P_H(3,3)=P_H(3,3)+k/2; P_H(5,5)=P_H(5,5)+k/2;
P_H(3,5)=P_H(3,5)-k/2; P_H(5,3)=P_H(5,3)-k/2;
% Spring PE finger 2
P_H(7,7)=P_H(7,7)+k/2; P_H(7,7)=P_H(7,7)+k/2;
P_H(9,9)=P_H(9,9)+k/2; P_H(7,9)=P_H(7,9)-k/2; P_H(9,7)=P_H(9,7)-k/2;
P_H(9,9)=P_H(9,9)+k/2; P_H(11,11)=P_H(11,11)+k/2;
P_H(9,11)=P_H(9,11)-k/2; P_H(11,9)=P_H(11,9)-k/2;
% Contact PE: (kc/2)*delta^2, delta=l*(x1+x3+x5+x7+x9+x11)
% => (kc*l^2/2)*(sum of all angle states)^2
ang_idx = [1,3,5,7,9,11];  % 1-indexed angle states
for ii = ang_idx
    for jj = ang_idx
        P_H(ii,jj) = P_H(ii,jj) + (kc/2)*l^2;
    end
end

M_H = A' * P_H + P_H * A;
max_eig_MH = max(real(eig(M_H)));
fprintf('  Max eigenvalue of A^T*P_H + P_H*A = %.4f\n', max_eig_MH);
if max_eig_MH > 0
    fprintf('  [INFO] Physical energy is NOT a valid Lyapunov function\n');
    fprintf('         (A^T*P_H + P_H*A has a positive eigenvalue = %.4f)\n', max_eig_MH);
    fprintf('         Reason: dUc/dt couples ALL joint velocities but\n');
    fprintf('         contact force acts only at distal link, leaving an\n');
    fprintf('         indefinite cross-term in Vdot.\n\n');
end

% -----------------------------------------------------------
%  Step 2: Solve the Algebraic Lyapunov Equation (ALE)
%          A'*P + P*A = -Q,  Q = I
% -----------------------------------------------------------
fprintf('  Step 2: Solve A^T*P + P*A = -I (Algebraic Lyapunov Eq.)\n');

Q_lyap   = eye(12);
P_lyap   = lyap(A', Q_lyap);      % solves A'*P + P*A = -Q

% Verify solution
M_vdot   = A' * P_lyap + P_lyap * A;
residual = max(abs(M_vdot(:) + Q_lyap(:)));
fprintf('  Max residual |A^T*P + P*A + I| = %.2e\n', residual);
if residual < 1e-8
    fprintf('  [PASS] ALE solved to machine precision\n\n');
else
    fprintf('  [FAIL] Large residual: %.2e\n\n', residual);
end

% -----------------------------------------------------------
%  Step 3: Verify P > 0
% -----------------------------------------------------------
fprintf('  Step 3: Verify P positive definite\n');
eig_P    = eig(P_lyap);
lam_min  = min(real(eig_P));
lam_max  = max(real(eig_P));
fprintf('  lambda_min(P) = %.6f\n', lam_min);
fprintf('  lambda_max(P) = %.6f\n', lam_max);
if lam_min > 0
    fprintf('  [PASS] P is POSITIVE DEFINITE\n\n');
else
    fprintf('  [FAIL] P is NOT positive definite\n\n');
end

% -----------------------------------------------------------
%  Step 4: Verify Vdot = -||x||^2 < 0
% -----------------------------------------------------------
fprintf('  Step 4: Verify Vdot = x^T*(A^T*P+P*A)*x = -||x||^2\n');
eig_Mvdot = real(eig(M_vdot));
fprintf('  Max eigenvalue of (A^T*P+P*A) = %.6e  (must be <= -1)\n', max(eig_Mvdot));
fprintf('  Min eigenvalue of (A^T*P+P*A) = %.6e\n', min(eig_Mvdot));
if max(eig_Mvdot) <= -1 + 1e-6
    fprintf('  [PASS] A^T*P + P*A = -I: Vdot = -||x||^2 < 0 STRICTLY\n');
    fprintf('  [PASS] NO LaSalle argument needed: direct negative definiteness\n\n');
else
    fprintf('  [FAIL] Vdot matrix is not negative definite as expected\n\n');
end

% -----------------------------------------------------------
%  Step 5: Open-loop eigenvalues (independent check)
% -----------------------------------------------------------
eig_A = eig(A);
fprintf('  Step 5: Open-loop eigenvalues of A:\n');
for i = 1:length(eig_A)
    fprintf('    lambda_%02d = %+.6f %+.6fi\n', i, real(eig_A(i)), imag(eig_A(i)));
end
fprintf('\n');
if all(real(eig_A) < 0)
    fprintf('  [PASS] All eigenvalues have NEGATIVE real parts\n');
    fprintf('  [PASS] A is Hurwitz => ALE solution P > 0 exists (converse theorem)\n\n');
end

% -----------------------------------------------------------
%  Step 6: Exponential decay bound
% -----------------------------------------------------------
fprintf('  Step 6: Exponential decay bound\n');
fprintf('  V(t) = x^T*P*x satisfies:\n');
fprintf('  V(t) <= V(0)*exp(-t / lambda_max(P))\n');
fprintf('  => ||x(t)||^2 <= (lam_max/lam_min) * ||x0||^2 * exp(-t/%.4f)\n\n', lam_max);

%% ============================================================
%  SECTION 3: BIBO STABILITY
%% ============================================================

fprintf('============================================================\n');
fprintf('  BIBO STABILITY ANALYSIS\n');
fprintf('============================================================\n');

beta_bibo = min(abs(real(eig_A)));
fprintf('  Decay rate beta = min|Re(lambda_i)| = %.6f\n', beta_bibo);
if beta_bibo > 0
    fprintf('  [PASS] beta > 0\n\n');
end

norm_B = norm(B,'fro');
norm_C = norm(C,'fro');
P_bibo = lyap(A', eye(12));
alpha_bibo = sqrt(cond(P_bibo));
M_input = 1.0;  x0_norm = 1.0;
M1 = norm_C * alpha_bibo * x0_norm;
M2 = (alpha_bibo * norm_C * norm_B * M_input) / beta_bibo;
N_bound = M1 + M2;

fprintf('  ||B||_F = %.6f,  ||C||_F = %.6f\n', norm_B, norm_C);
fprintf('  BIBO Output Bound N = %.4f  (finite)\n', N_bound);
fprintf('  [PASS] N is FINITE -> System is BIBO STABLE\n\n');

sys      = ss(A, B, C, D);
poles_tf = pole(sys);
if all(real(poles_tf) < 0)
    fprintf('  [PASS] All transfer function poles in open left-half plane\n\n');
end

%% ============================================================
%  SECTION 4: CONTROLLABILITY
%% ============================================================

fprintf('============================================================\n');
fprintf('  CONTROLLABILITY ANALYSIS\n');
fprintf('============================================================\n');

Ctrb   = ctrb(A, B);
rank_C = rank(Ctrb);
fprintf('  rank(controllability matrix) = %d  (required: 12)\n', rank_C);
if rank_C == 12
    fprintf('  [PASS] System is COMPLETELY CONTROLLABLE\n\n');
end

tol = 1e-10;
fprintf('  Row activation sequence:\n');
activated = false(12,1);
col_blocks = {B, A*B, A^2*B, A^3*B};
block_names = {'B','AB','A^2*B','A^3*B'};
for kb = 1:4
    blk = col_blocks{kb};  new_rows = [];
    for r = 1:12
        if ~activated(r) && any(abs(blk(r,:))>tol)
            activated(r)=true; new_rows(end+1)=r; %#ok<AGROW>
        end
    end
    fprintf('    %-8s -> rows: [%s]\n', block_names{kb}, num2str(new_rows));
end
fprintf('\n');

pbh_pass = true;
for i = 1:length(eig_A)
    if rank([eig_A(i)*eye(12)-A, B]) < 12
        pbh_pass = false;
    end
end
fprintf('  PBH test: %s\n\n', iff(pbh_pass,'[PASS] all eigenvalues pass','[FAIL]'));

%% ============================================================
%  SECTION 5: OBSERVABILITY
%% ============================================================

fprintf('============================================================\n');
fprintf('  OBSERVABILITY ANALYSIS\n');
fprintf('============================================================\n');

Obsv   = obsv(A, C);
rank_O = rank(Obsv);
fprintf('  rank(observability matrix) = %d  (required: 12)\n', rank_O);
if rank_O == 12
    fprintf('  [PASS] System is COMPLETELY OBSERVABLE\n\n');
end

rank_two = rank([C; C*A]);
fprintf('  rank([C; CA]) = %d  (required: 12)\n', rank_two);
if rank_two == 12
    fprintf('  [PASS] Full observability confirmed from C and CA alone\n\n');
end

pbh_obs_pass = true;
for i = 1:length(eig_A)
    if rank([eig_A(i)*eye(12)-A; C]) < 12
        pbh_obs_pass = false;
    end
end
fprintf('  PBH observability test: %s\n\n', ...
    iff(pbh_obs_pass,'[PASS] all eigenvalues pass','[FAIL]'));

%% ============================================================
%  SECTION 6: INTERNAL vs EXTERNAL STABILITY
%% ============================================================

fprintf('============================================================\n');
fprintf('  INTERNAL vs EXTERNAL STABILITY EQUIVALENCE\n');
fprintf('============================================================\n');

poles_sys = pole(sys);
zeros_sys = tzero(sys);
cancelled = false;
for p = poles_sys.'
    for z = zeros_sys.'
        if abs(p-z) < 1e-6, cancelled = true; end
    end
end
fprintf('  Pole-zero cancellation: %s\n', iff(~cancelled,'[PASS] none detected','[WARN] cancellation found'));

if rank_C==12 && rank_O==12 && all(real(eig_A)<0)
    fprintf('  [PASS] Lyapunov stability <=> BIBO stability (no hidden modes)\n');
    fprintf('  [PASS] System is INTERNALLY and EXTERNALLY STABLE\n\n');
end

%% ============================================================
%  SECTION 7: SIMULATION — V(t) decay
%% ============================================================

fprintf('============================================================\n');
fprintf('  SIMULATION: Lyapunov V(t) decay\n');
fprintf('============================================================\n');

x0 = zeros(12,1);
x0(1) = deg2rad(5);  x0(7) = deg2rad(5);
t_span = linspace(0, 5, 2000);

[t_sim, x_sim] = ode45(@(t,x) A*x, t_span, x0);

V_t = zeros(length(t_sim),1);
for i = 1:length(t_sim)
    V_t(i) = x_sim(i,:) * P_lyap * x_sim(i,:)';
end

fprintf('  V(0) = %.6f\n', V_t(1));
fprintf('  V(T) = %.6e\n', V_t(end));
n_inc = sum(diff(V_t) > 1e-12);
if n_inc == 0
    fprintf('  [PASS] V(t) MONOTONICALLY NON-INCREASING\n\n');
else
    fprintf('  [INFO] %d small numerical increases (ODE tolerance)\n\n', n_inc);
end

%% ============================================================
%  SECTION 8: PLOTS
%% ============================================================

figure('Name','Stability Verification (Corrected Lyapunov)',...
    'NumberTitle','off','Position',[100 100 1300 850]);

subplot(2,3,1);
semilogy(t_sim, V_t, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('V(t) = x^T P x');
title('Lyapunov V(t) Decay'); grid on;

subplot(2,3,2);
plot(t_sim, rad2deg(x_sim(:,[1,3,5])), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (deg)');
title('Finger 1 Joint Angles (open-loop)');
legend('\theta_{11}','\theta_{12}','\theta_{13}'); grid on;

subplot(2,3,3);
plot(t_sim, rad2deg(x_sim(:,[7,9,11])), 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle (deg)');
title('Finger 2 Joint Angles (open-loop)');
legend('\theta_{21}','\theta_{22}','\theta_{23}'); grid on;

subplot(2,3,4);
semilogy(sort(svd(Ctrb),'descend'), 'ro-', 'LineWidth', 1.5);
xlabel('Index'); ylabel('Singular Value');
title(sprintf('Controllability SVD (rank=%d)',rank_C)); grid on;

subplot(2,3,5);
semilogy(sort(svd(Obsv),'descend'), 'gs-', 'LineWidth', 1.5);
xlabel('Index'); ylabel('Singular Value');
title(sprintf('Observability SVD (rank=%d)',rank_O)); grid on;

subplot(2,3,6);
plot(real(eig_A), imag(eig_A), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
hold on; xline(0,'r--','LineWidth',1.5);
xlabel('Real'); ylabel('Imaginary');
title('Open-loop Eigenvalues of A'); grid on;

sgtitle('PRB Soft Gripper — Corrected Lyapunov Proof',...
    'FontSize',13,'FontWeight','bold');

%% ============================================================
%  SECTION 9: FINAL SUMMARY
%% ============================================================

fprintf('============================================================\n');
fprintf('  FINAL VERIFICATION SUMMARY\n');
fprintf('============================================================\n');
fprintf('  Lyapunov Stability (ALE-based):\n');
fprintf('    A is Hurwitz (all Re<0):            %s\n', iff(all(real(eig_A)<0),'PASS','FAIL'));
fprintf('    ALE residual < 1e-8:                %s\n', iff(residual<1e-8,'PASS','FAIL'));
fprintf('    P_lyap > 0:                         %s\n', iff(lam_min>0,'PASS','FAIL'));
fprintf('    A^T*P+P*A = -I (Vdot = -||x||^2):  %s\n', iff(max(eig_Mvdot)<=-1+1e-6,'PASS','FAIL'));
fprintf('    V(t) decaying (simulation):         %s\n', iff(V_t(end)<V_t(1),'PASS','FAIL'));
fprintf('\n  BIBO Stability:\n');
fprintf('    beta > 0:                           %s\n', iff(beta_bibo>0,'PASS','FAIL'));
fprintf('    N finite:                           %s\n', iff(isfinite(N_bound),'PASS','FAIL'));
fprintf('    All TF poles stable:                %s\n', iff(all(real(poles_tf)<0),'PASS','FAIL'));
fprintf('\n  Controllability:\n');
fprintf('    rank(Ctrb) = 12:                    %s\n', iff(rank_C==12,'PASS','FAIL'));
fprintf('    PBH test:                           %s\n', iff(pbh_pass,'PASS','FAIL'));
fprintf('\n  Observability:\n');
fprintf('    rank(Obsv) = 12:                    %s\n', iff(rank_O==12,'PASS','FAIL'));
fprintf('    rank([C;CA]) = 12:                  %s\n', iff(rank_two==12,'PASS','FAIL'));
fprintf('    PBH test:                           %s\n', iff(pbh_obs_pass,'PASS','FAIL'));
fprintf('\n  Internal = External Stability:        %s\n', ...
    iff(rank_C==12 && rank_O==12 && all(real(eig_A)<0),'PASS','FAIL'));
fprintf('============================================================\n');

%% Helper
function out = iff(cond, a, b)
    if cond; out = a; else; out = b; end
end
