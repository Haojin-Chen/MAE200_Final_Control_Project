clear all; clc
%{
The following code is used to design an optimal controller and estimator for 2 pendulem of
different length and mass on a single cart. The only acutator on this car
the motor that drives the linear position of the cart. 

We are attempting to swing up both the pendulem from the inital position of
all pendulem down to both pendulem up and stabilize

%}

%% Determine u_k and state
clear all; clc
load_prev = 1;
%{
Use code provided in NR Ch21 to generate x_k state, u_k, using adjoint
method for the swing up phase. This phase last about T = 3 seconds because 

0) Within this code, tune QT to drive error at T close to zero
1) Warm start: use an inital guess of u_k being all zeros
2) iterate the provided function with a new u_k until the state reach zero

%}

% T = 3; %approx sqrt(g/l) 
% 
% %initial guess of u_k being all zero 
% [u_k(:,1),x_k] = OptimizeTraj(T);
% 
% state_at_T(:,1) = x_k(1:6, end);
% 
% number_of_trial = 1;
% trial_number = 1:number_of_trial
% for i = trial_number
%     [u_k(:,i+1), x_k] = OptimizeTraj(T, u_k(:,i))
%     state_at_T(:,i+1) = x_k(1:6, end)
% end

%load previous trajectory if it is better
if load_prev == 1 
    clear all; clc;
    load("SavePoint3.mat")
    figure(1); clf; subplot(3,1,1); plot(t,x_k(1,:),'r-',t,x_k(2,:),'b-',t,x_k(3,:),'g-'); hold on;
              subplot(3,1,2); plot(t,x_k(4,:),'r-',t,x_k(5,:),'b-',t,x_k(6,:),'g-'); hold on;
              subplot(3,1,3); plot(t,u_k,'r--'); hold on;
end



%% Find A(t), E(t), B(t) for Each x_k

%compute for the state matrix E_t, A_t, B_t, C_t, D_t wrt time per timestep
%h

s.h=0.01; s.N=T/s.h; s.mc=10; t=[0:s.N]*s.h;               % STEP 0: initialize simulation
s.m1=1; s.L1=1;    s.ell1=s.L1; s.I1=s.m1*s.ell1^2/3; % system, & derived parameters
s.m2=0.5; s.L2=0.5;  s.ell2=s.L2; s.I2=s.m2*s.ell2^2/3; alpha=0.1;
s.B=[0; 0; 0; 1; 0; 0]; s.Q=diag([0 0 0 0 0 0]); s.R=0; s.QT=diag([5 40 10 .1 60 10]);

for i = 1:length(x_k)
    A_temp = Compute_A(x_k(:,i) ,s);
    A_t(i) = {A_temp};
end

for i = 1:length(x_k)
    E_temp = Compute_E(x_k(1:6,i),s);
    E_t(i) = {E_temp};
end

for i = 1:length(x_k)
    B_temp = [0; 0; 0; 1; 0; 0];
    B_t(i) = {B_temp};
end

for i = 1:length(x_k)
    C_temp = [[1 0 0  0   0   0];...
              [0 1 0  0   0   0];...
              [0 0 1  0   0   0];...
              [0 0 0  0   0   0];...
              [0 0 0  0   0   0];...
              [0 0 0  0   0   0]];
    C_t(i) = {C_temp};
end

%% Compute K_t using A, E, B, C

%{
 compute controller gain K(t) using differential ricatti equation, which
 can be solved using RK4 marching

%}
h = s.h; %timestep

R = .1; %controller cost
%Q = eye(6) %state deviation cost
 Q = [[1 0 0  0   0   0];...
      [0 5 0  0   0   0];...
      [0 0 5  0   0   0];...
      [0 0 0  1   0   0];...
      [0 0 0  0   1   0];...
      [0 0 0  0   0   1]];

X_t = {}; %allocate cell frame to store X at every time frame
X_t(length(u_k)) = {eye(6)} %X at T

%backward march to find X_t using RK4
for i = length(u_k):-1:2

f1 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i},R,Q);
f2 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i} - f1*h/2,R,Q);
f3 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i} - f2*h/2,R,Q);
f4 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i} - f3*h,R,Q);

X_t{i-1} = X_t{i} - h*(f1/6 + f2/3 + f3/3 + f4/6);

end

%find K_t using X_t
for i = 1:length(u_k)
    K_t(i) = {-inv(R)*B_t{i}'*X_t{i}}
end

%% Compute L_t using A, E, B, C

%{
 compute observer gain K(t) using differential ricatti equation, which
 can be solved using RK4 marching

%}

Q1 = eye(6); %covariance of state disturbance
Q2 = eye(6); %convariance of sensor noise 

P_t = {}; %allocate space for storage of P

P_t(1) = {eye(6)}; %intial estimation covariance

%forward march to find L(t) using RK4
for i = 1:1:length(u_k)

f1 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i},R,Q);
f2 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i} + f1*h/2,R,Q);
f3 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i} + f2*h/2,R,Q);
f4 = KRiccati(E_t{i},A_t{i},B_t{i},X_t{i} + f3*h,R,Q);

P_t{i+1} = P_t{i} + h*(f1/6 + f2/3 + f3/3 + f4/6);

end

%find P_t
for i = 1:length(u_k)
    L_t(i) = {-inv(P_t{i})*C_t{i}'*inv(Q2)}
end

%% Validation of K(t)

eigControl = {};
ControlSys = {};
for i = 1:length(u_k)
    ControlSys(i) = {inv(E_t{i})*A_t{i} + inv(E_t{i})* B_t{i}*K_t{i}};
   [evecK, eigControl{i}] = eig(ControlSys{i});
end

eigEstimate = {};
EstimateSys = {};
for i = 1:length(u_k)
    EstimateSys(i) = {A_t{i} + L_t{i}*C_t{i}};
   [evecL, eigEstimate{i}] = eig(EstimateSys{i});
end

%forward march of state x using control system
x_state(:, 1) = [0 pi pi 0 0 0]'; 

for i = 1:1:length(u_k)-1

f1 = ControlSys{i} * x_state(:, i);
f2 = ControlSys{i} * (x_state(:, i) + f1*h/2);
f3 = ControlSys{i} * (x_state(:, i) + f2*h/2);
f4 = ControlSys{i} * (x_state(:, i) + f3*h);

x_state(:,i+1) = x_state(i) + h*(f1/6 + f2/3 + f3/3 + f4/6);
%x_state(:,i+1) = x_state(i) + h*f1
end

for i = 1:length(u_k)
    u_input(i) = K_t{i} * x_state(:,i);
end

figure(1); subplot(3,1,1); plot(t,x_state(1,:),'r-',t,x_state(2,:),'b-',t,x_state(3,:),'g-');
                    subplot(3,1,2); plot(t,x_state(4,:),'r-',t,x_state(5,:),'b-',t,x_state(6,:),'g-');
                 subplot(3,1,3); plot(t, u_input);


