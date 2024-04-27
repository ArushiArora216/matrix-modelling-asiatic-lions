%Defining all parameters using arrays.
%p is the array of survivorship and d is the array for resepctive stage duration
p=[0.57 0.87 0.90 0.90];
d=[1 1 1 13];

%Preallocate P and G; also define F
P=zeros(1,4);
G=zeros(1,4);
F = [0 0 0 1.7445];

%Calculating each of the P_i's and G_i's
for i=1:4
    P(i)=((1-p(i)^(d(i)-1))/(1-p(i)^d(i)))*p(i);
end

for i=1:3
    G(i)=((p(i)^d(i))*(1-p(i)))/(1-p(i)^d(i));
end

%Create population projection matrix A
A=[P(1) F(2) F(3) F(4); G(1) P(2) 0 0; 0 G(2) P(3) 0; 0 0 G(3) P(4)];
disp(A)
lambda = max(eig(A));
fprintf('The initial estimate of the intrinsic growth rate of the population is %f \n', lambda)

% Calculating the intrinsic growth rate of the population
[V,D,W] = eig(A);
[lambda, ind] = max(diag(D));
fprintf('The intrinsic growth rate of the population is %f \n', lambda)