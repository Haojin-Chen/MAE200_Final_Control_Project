%% Compute State matrixes
s.h=0.01; s.N=T/s.h; s.mc=10; t=[0:s.N]*s.h;               % STEP 0: initialize simulation
s.m1=1; s.L1=1;    s.ell1=s.L1; s.I1=s.m1*s.ell1^2/3; % system, & derived parameters
s.m2=0.5; s.L2=0.5;  s.ell2=s.L2; s.I2=s.m2*s.ell2^2/3; alpha=0.1;
s.B=[0; 0; 0; 1; 0; 0]; s.Q=diag([0 0 0 0 0 0]); s.R=0; s.QT=diag([5 40 10 .1 60 10]);

g = 9.81;

E44 = s.m1 + s.m2 + s.mc; E45 = -s.m1*s.L1; E46 = -s.m2*s.L2;
E54 = -s.m1*s.L1 ; E55 = s.I1 + s.m1*s.L1^2;
E64 = -s.m2*s.L2 ; E66 = s.I2 + s.m2*s.L2^2;


E_bar = [[1 0 0  0   0   0];...
     [0 1 0  0   0   0];
     [0 0 1  0   0   0];
     [0 0 0 E44 E45 E46];
     [0 0 0 E54 E55  0 ];
     [0 0 0 E64  0  E66]]

A52 = s.m1*g*s.L1
A63 = s.m2*g*s.L2
A_bar = [[0  0  0  1   0   0];...
         [0  0  0  0   1   0];
         [0  0  0  0   0   1];
         [0  0  0  0   0   0];
         [0 A52 0 0   0   0];
         [0  0  A63  0   0   0]]

B_bar = [0 0 0  1  0  0]'
%% Compute normal A, B, and C

A = E_bar^-1 * A_bar;
B = E_bar^-1 * B_bar;
C = [[1 0 0  0   0   0];...
     [0 1 0  0   0   0];...
     [0 0 1  0   0   0];...
     [0 0 0  0   0   0];...
     [0 0 0  0   0   0];...
     [0 0 0  0   0   0]];...
% 
% C = [[1 0 0  0   0   0];...
%      [0 1 0  0   0   0];...
%      [0 0 1  0   0   0]];

D = [0 0 0 0 0 0]'

E = eye(6);
%% Solve for Controller 
 
% 
% Solve raccatti K

Q = eye(6)
R = 1
[X,K_inf,L] = icare(A,B,Q,R)

% Solve for raccatti L using iCARE (use a-bk)

Q1 = eye(6)
Q2 = eye(6)
[P,L_inf,L] = icare(A', C', Q1, Q2)



%%

(A-L_inf*C)



eig()