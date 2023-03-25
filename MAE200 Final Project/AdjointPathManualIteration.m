clear all; clc
load("SavePoint3.mat")
figure(1); clf; subplot(3,1,1); plot(t,x_k(1,:),'r-',t,x_k(2,:),'b-',t,x_k(3,:),'g-');
              subplot(3,1,2); plot(t,x_k(4,:),'r-',t,x_k(5,:),'b-',t,x_k(6,:),'g-');
              subplot(3,1,3); plot(t,u_k,'r--');
%% Manual Tuning u_k and x_k
T = 3;
[u_k, x_k] = OptimizeTraj(T, u_k);


t = 0:0.01:T
figure(1); clf; subplot(3,1,1); plot(t,x_k(1,:),'r-',t,x_k(2,:),'b-',t,x_k(3,:),'g-');
              subplot(3,1,2); plot(t,x_k(4,:),'r-',t,x_k(5,:),'b-',t,x_k(6,:),'g-');
              subplot(3,1,3); plot(t,u_k,'r--'); 


%%
%save("SavePoint.mat")
%save("SavePoint2.mat")
save("SavePoint4.mat")