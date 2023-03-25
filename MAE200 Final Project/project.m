clear all; clc
T = 3; %approx sqrt(g/l) 

theta1_tolerance = 5*(pi/180);
theta2_tolerance = 5*(pi/180);

%initial guess of u_k being all zero 
[u_k(:,1),x_k] = Optimize_u_k(T)

state_at_T(:,1) = x_k(1:6, end)

number_of_trial = 16
trial_number = 1:number_of_trial -1 
for i = trial_number
    [u_k(:,i+1), x_k] = Optimize_u_k(T, u_k(:,i))
    state_at_T(:,i+1) = x_k(1:6, end)
end

figure(2) 
subplot(2,1,1); plot(1:number_of_trial,state_at_T(1,:),'r-',...
                     1:number_of_trial,state_at_T(2,:),'b-',...
                     1:number_of_trial,state_at_T(3,:),'g-',...
                     1:number_of_trial,theta1_tolerance*ones(1,number_of_trial), 'b--',...
                     1:number_of_trial,0*ones(1,number_of_trial), 'r--');
subplot(2,1,2); plot(1:number_of_trial,state_at_T(4,:),'r-',...
                     1:number_of_trial,state_at_T(5,:),'b-',...
                     1:number_of_trial,state_at_T(6,:),'g-');
