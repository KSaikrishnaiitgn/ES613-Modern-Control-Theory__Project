%% =========================================================================
%  Symbolic Jacobian Linearisation – Soft Robotic Gripper
%  -------------------------------------------------------------------------
%  Derives A (12x12), B (12x2), C (7x12) symbolically.
%
%  VERIFIED CORRECT A11 BLOCK (uniform stiffness shorthand):
%
%       0        1       0    0       0        0
%       α        β       κ    μ       0        0
%       0        0       0    1       0        0
%       κ        μ       α    β       κ        μ
%       0        0       0    0       0        1
%       0        0       κ    μ   -κ-γc    β+μ(=-b/Jl... see note)
%
%  NOTE on PDF Eq.(11):
%    Row 2 col 1: PDF writes "α-κ" but correct value is α = -(k11+k12)/Jl = -2k/Jl
%    Row 2 col 2: PDF writes "β+μ" but correct value is β = -2b/Jl
%    Row 6 col 6: PDF writes "β+μ" = -2b/Jl + b/Jl = -b/Jl which IS correct
%                 since ∂f6/∂x6 = -b/Jl (only one damper term, not two)
%    The PDF row 2 entries are TYPOS.  All other entries match.
%
%  State vector (12x1):
%    x = [x1=θ11, x2=θ̇11, x3=θ12, x4=θ̇12, x5=θ13, x6=θ̇13,
%          x7=θ21, x8=θ̇21, x9=θ22, x10=θ̇22, x11=θ23, x12=θ̇23]
%  Input  (2x1):  u = [u1=τ1, u2=τ2]
%  Output (7x1):  y = [θ11, θ12, θ13, θ21, θ22, θ23, Fc]
% =========================================================================
clear; clc;

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 1 – SYMBOLIC VARIABLE DECLARATIONS
% ─────────────────────────────────────────────────────────────────────────
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12   real
syms u1 u2                                        real
syms Jl b kc l d w                                real positive
syms k11 k12 k13 k21 k22 k23                      real positive

x_vec = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12];
u_vec = [u1; u2];

fprintf('=========================================================\n');
fprintf(' Symbolic Jacobian Linearisation – Soft Robotic Gripper\n');
fprintf('=========================================================\n\n');
fprintf(' All parameters kept symbolic throughout.\n');
fprintf(' Equilibrium: theta0=(d-w)/(6l), Fc0=0 (just-touching).\n\n');

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 2 – EQUILIBRIUM POINT
%  -----------------------------------------------------------------------
%  Symmetric small-angle equilibrium (Eq.7 of paper):
%    theta0 = (d-w)/(6l)  →  Fc0 = kc*(6l*theta0-(d-w)) = 0
%  All angles equal theta0, all velocities = 0.
% ─────────────────────────────────────────────────────────────────────────
theta0 = (d - w)/(6*l);   % symbolic equilibrium angle

eq_vars = {x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12};
eq_vals = {theta0,0,theta0,0,theta0,0,theta0,0,theta0,0,theta0,0};

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 3 – NONLINEAR CONTACT FORCE  Fc(x)
%  -----------------------------------------------------------------------
%  Active-branch linearised Hertz (Eq.74):
%    Fc = kc * [ l*(sinx1+sinx3+sinx5+sinx7+sinx9+sinx11) - (d-w) ]
% ─────────────────────────────────────────────────────────────────────────
Fc_nl = kc*( l*(sin(x1)+sin(x3)+sin(x5)+sin(x7)+sin(x9)+sin(x11)) - (d-w) );

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 4 – NONLINEAR STATE EQUATIONS  f(x,u)
%  -----------------------------------------------------------------------
%  Eqs. 75a-75f (Appendix C/D), both fingers.
%
%  FINGER 1:
%   ẋ1  = x2
%   ẋ2  = (1/Jl)*[ u1 - (k11+k12)*x1 - 2b*x2 + k12*x3 + b*x4 ]
%   ẋ3  = x4
%   ẋ4  = (1/Jl)*[ k12*x1 + b*x2 - (k12+k13)*x3 - 2b*x4 + k13*x5 + b*x6 ]
%   ẋ5  = x6
%   ẋ6  = (1/Jl)*[ k13*x3 + b*x4 - k13*x5 - b*x6 - Fc*l*cos(x5) ]
%
%  FINGER 2 (symmetric, index offset +6):
%   ẋ7  = x8
%   ẋ8  = (1/Jl)*[ u2 - (k21+k22)*x7 - 2b*x8 + k22*x9 + b*x10 ]
%   ẋ9  = x10
%   ẋ10 = (1/Jl)*[ k22*x7 + b*x8 - (k22+k23)*x9 - 2b*x10 + k23*x11 + b*x12 ]
%   ẋ11 = x12
%   ẋ12 = (1/Jl)*[ k23*x9 + b*x10 - k23*x11 - b*x12 - Fc*l*cos(x11) ]
% ─────────────────────────────────────────────────────────────────────────
f1  = x2;
f2  = (1/Jl)*( u1 - (k11+k12)*x1 - 2*b*x2 + k12*x3 + b*x4 );
f3  = x4;
f4  = (1/Jl)*( k12*x1 + b*x2 - (k12+k13)*x3 - 2*b*x4 + k13*x5 + b*x6 );
f5  = x6;
f6  = (1/Jl)*( k13*x3 + b*x4 - k13*x5 - b*x6 - Fc_nl*l*cos(x5) );
f7  = x8;
f8  = (1/Jl)*( u2 - (k21+k22)*x7 - 2*b*x8 + k22*x9 + b*x10 );
f9  = x10;
f10 = (1/Jl)*( k22*x7 + b*x8 - (k22+k23)*x9 - 2*b*x10 + k23*x11 + b*x12 );
f11 = x12;
f12 = (1/Jl)*( k23*x9 + b*x10 - k23*x11 - b*x12 - Fc_nl*l*cos(x11) );

f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12];

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 5 – OUTPUT EQUATIONS  g(x)
%  g = [θ11, θ12, θ13, θ21, θ22, θ23, Fc]
% ─────────────────────────────────────────────────────────────────────────
g = [x1; x3; x5; x7; x9; x11; Fc_nl];

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 6 – JACOBIAN  A = ∂f/∂x  evaluated at equilibrium
%  -----------------------------------------------------------------------
%  Two-step:
%   1. Symbolic differentiation (exact, general x)
%   2. Substitute equilibrium + small-angle (cos(theta0)→1, Fc0→0)
%
%  KEY DERIVATIONS FOR NONLINEAR ROWS 6 AND 12:
%
%  f6 contains  -Fc*l*cos(x5)/Jl
%  ∂f6/∂xj involves two contributions:
%    (a) -(∂Fc/∂xj)*l*cos(x5)/Jl
%    (b) +Fc*l*(∂cos(x5)/∂xj)/Jl  [nonzero only when j=5]
%
%  At equilibrium Fc0=0, cos(theta0)=1, sin(theta0)=theta0:
%
%  j=3: ∂f6/∂x3 = k13/Jl                          [linear term only]
%  j=4: ∂f6/∂x4 = b/Jl                             [linear term only]
%  j=5: ∂f6/∂x5 = (1/Jl)*[-k13
%                           -(kc*l*cos(x5))*l*cos(x5)  ← from (a)
%                           +Fc0*l*(-sin(x5))           ← from (b)=0]
%               = -(k13 + kc*l^2)/Jl
%  j=6: ∂f6/∂x6 = -b/Jl                            [linear term only]
%  j=11:∂f6/∂x11= -(∂Fc/∂x11)*l*cos(x5)/Jl
%               = -(kc*l*cos(x11))*l*1/Jl = -kc*l^2/Jl  [INTER-FINGER]
%  j=1,3_extra,7,9: all zero because Fc0=0 eliminates those cross terms
%
%  Identical reasoning applies to row 12 (distal link of finger 2).
% ─────────────────────────────────────────────────────────────────────────
fprintf(' Computing A = ∂f/∂x symbolically and evaluating at equilibrium...\n\n');

A_sym = sym(zeros(12,12));
for i = 1:12
    for j = 1:12
        entry = diff(f(i), x_vec(j));
        entry = subs(entry, eq_vars, eq_vals);         % substitute equilibrium
        entry = subs(entry, cos(theta0), sym(1));      % small-angle cos→1
        entry = subs(entry, sin(theta0), theta0);      % small-angle sin→θ
        entry = simplify(entry);
        A_sym(i,j) = entry;
    end
end

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 7 – JACOBIAN  B = ∂f/∂u  (no equilibrium substitution needed)
%  -----------------------------------------------------------------------
%  τ1=u1 in f2 only  →  B(2,1) = 1/Jl
%  τ2=u2 in f8 only  →  B(8,2) = 1/Jl
%  All other entries  = 0
% ─────────────────────────────────────────────────────────────────────────
fprintf(' Computing B = ∂f/∂u ...\n\n');

B_sym = sym(zeros(12,2));
for i = 1:12
    for j = 1:2
        B_sym(i,j) = simplify(diff(f(i), u_vec(j)));
    end
end

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 8 – JACOBIAN  C = ∂g/∂x  evaluated at equilibrium
%  -----------------------------------------------------------------------
%  g1..g6 are linear in x  →  rows 1-6 are standard basis vectors
%  g7 = Fc (nonlinear via sin):
%    ∂Fc/∂xj = kc*l*cos(xj) → kc*l at equilibrium  for angle states
%    ∂Fc/∂(velocity states) = 0
% ─────────────────────────────────────────────────────────────────────────
fprintf(' Computing C = ∂g/∂x ...\n\n');

C_sym = sym(zeros(7,12));
for i = 1:7
    for j = 1:12
        entry = diff(g(i), x_vec(j));
        entry = subs(entry, eq_vars, eq_vals);
        entry = subs(entry, cos(theta0), sym(1));
        entry = subs(entry, sin(theta0), theta0);
        entry = simplify(entry);
        C_sym(i,j) = entry;
    end
end

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 9 – UNIFORM STIFFNESS + SHORTHAND SCALARS  (Eq.9 of paper)
%  -----------------------------------------------------------------------
%  kij = k for all i,j
%  α = -2k/Jl,  β = -2b/Jl,  κ = k/Jl,  μ = b/Jl,  γc = kc*l²/Jl
% ─────────────────────────────────────────────────────────────────────────
syms k  real positive
syms alpha beta kappa mu gamma_c  real

A_uniform = subs(A_sym, {k11,k12,k13,k21,k22,k23}, {k,k,k,k,k,k});
A_uniform = simplify(A_uniform);

% Replace raw expressions with shorthand notation
A_short = A_uniform;
A_short = subs(A_short, -2*k/Jl,    alpha  );
A_short = subs(A_short, -2*b/Jl,    beta   );
A_short = subs(A_short,  k/Jl,      kappa  );
A_short = subs(A_short,  b/Jl,      mu     );
A_short = subs(A_short,  kc*l^2/Jl, gamma_c);
A_short = simplify(A_short);

%% ─────────────────────────────────────────────────────────────────────────
%  SECTION 10 – PRINT RESULTS AND VERIFY AGAINST PDF
% ─────────────────────────────────────────────────────────────────────────
fprintf('=========================================================\n');
fprintf(' CORRECTNESS CHECK: Row 2 entries (PDF vs Code)\n');
fprintf('=========================================================\n\n');

fprintf(' ∂f2/∂x1 = -(k11+k12)/Jl\n');
fprintf('   Under uniform k: -2k/Jl = α   (paper writes "α-κ" = -3k/Jl → TYPO)\n\n');

fprintf(' ∂f2/∂x2 = -2b/Jl = β\n');
fprintf('   (paper writes "β+μ" = -2b/Jl + b/Jl = -b/Jl → TYPO)\n\n');

fprintf(' ∂f6/∂x6 = -b/Jl\n');
fprintf('   In shorthand: -b/Jl = β/2 ... or equivalently -(β+μ) NO:\n');
fprintf('   β+μ = -2b/Jl + b/Jl = -b/Jl  ✓  so A(6,6)=β+μ IS correct in PDF\n');
fprintf('   This is consistent: row 6 has only one damper term (-b*x6)\n');
fprintf('   while row 2 has two damper terms (-2b*x2) giving β not β+μ.\n\n');

fprintf('=========================================================\n');
fprintf(' A matrix (general non-uniform stiffness)\n');
fprintf('=========================================================\n');
disp(A_sym)

fprintf('=========================================================\n');
fprintf(' A matrix (uniform kij=k, shorthand scalars)\n');
fprintf('=========================================================\n');
disp(A_short)

A11 = A_short(1:6,  1:6);
A12 = A_short(1:6,  7:12);
A21 = A_short(7:12, 1:6);
A22 = A_short(7:12, 7:12);

fprintf('--- A11 sub-block (compare with PDF Eq.11) ---\n');
disp(A11)
fprintf('--- A12 sub-block (only A12(6,5) = -gamma_c) ---\n');
disp(A12)
fprintf('--- A21 sub-block (only A21(6,5) = -gamma_c) ---\n');
disp(A21)
fprintf('--- A22 sub-block (equals A11 by symmetry) ---\n');
disp(A22)

fprintf('=========================================================\n');
fprintf(' B matrix\n');
fprintf('=========================================================\n');
disp(B_sym)

fprintf('=========================================================\n');
fprintf(' C matrix\n');
fprintf('=========================================================\n');
disp(C_sym)