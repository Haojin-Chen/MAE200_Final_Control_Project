h = 0.01;  % set the step size
t = 0:h:3;  % set the interval of t

% X = X(:);
X = zeros(length(t),36);   % X is transformed from "n"-by-"n" to "n^2"-by-1 for ease of calc
X(1,:) = ones(1,36);       % set the intial value for X
n = length(t)-1;
X_dot =@(t,X)(A'X + XA - XB(R^-1)B'X + Q); %insert function to be solved
for i = 1:n

    k1 = X_dot(t(i),X(i,:));
    k2 = X_dot(t(i)+.5h,X(i,:)+.5k1h);
    k3 = X_dot(t(i)+.5h,X(i,:)+.5k2h);
    k4 = X_dot(t(i)+h,X(i,:)+k3h);
    X(i+1,:) = X(i,:)+((k1+2k2+2k3+k4)/6)h;
    X = reshape(X, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"

end

for j=1:length(t)
    X{j} = reshape(X(j,:), size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
end