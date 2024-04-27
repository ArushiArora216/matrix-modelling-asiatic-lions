format long g
%Defining all parameters using arrays.
%p is the array of survivorship and d is the array for resepctive stage duration
surv1 = 0.57;
surv2 = 0.87;
surv3 = 0.90;
surv4 = 0.90;
p=[surv1 surv2 surv3 surv4];
d=[1 1 1 13];
F = [0 0 0 1.7445];

l=1.2575; % This is the lambda as obtained from Crouse's model.
gamma=zeros(1,4);
for i=1:4
    gamma(i)=(((p(i)/l)^d(i))-((p(i)/l)^(d(i)-1)))/(((p(i)/l)^d(i))-1);
end

P=zeros(1,4);
G=zeros(1,4);

for i=1:4
    P(i)=p(i)*(1-gamma(i));
end

for i=1:3
    G(i)=p(i)*gamma(i);
end
G(4)=0;

%Create population projection matrix A
A=[P(1) F(2) F(3) F(4); G(1) P(2) 0 0; 0 G(2) P(3) 0; 0 0 G(3) P(4)];
fprintf('The population projection matrix is: \n')
disp(A)

% Calculating the intrinsic growth rate of the population
[V,D,W] = eig(A);
[lambda, ind] = max(diag(D));
fprintf('The intrinsic growth rate of the population is %f \n', lambda)

% Finding w (stable stage distribution)
right_eigenvector=V(:, ind);
right_eigenvector_sum = sum(right_eigenvector);
w=(100/right_eigenvector_sum)*right_eigenvector;
fprintf('The stable stage distribution of the population is: \n')
disp(w) 

% Finding v (reproductive output)
left_eigenvector=W(:, ind);
v=(1/left_eigenvector(1)).*left_eigenvector;
fprintf('The stage-specific reproductive output of the population is: \n')
disp(v)

% Sensitivity analyses
% Sensitivity w.r.t. the stage survival probability P_i
Sensitivity_P=zeros(1,4);
for i=1:4
    Sensitivity_P(i)=(A(i,i)/lambda)*((v(i)*w(i))/dot(v,w));
end

%Sensitivity w.r.t. G_i
Sensitivity_G=zeros(1,3);
for i=1:3
    Sensitivity_G(i)=(A(i+1,i)/lambda)*((v(i+1)*w(i))/dot(v,w));
end
Sensitivity_G(4)=0;

%Sensitivity w.r.t. fecundity
Sensitivity_F=zeros(1,4);
Sensitivity_F(4)=(F(4)/lambda)*((v(1)*w(3))/dot(v,w));

categories = {'Cubs','Juveniles','Subadults', 'Adults'};

plot(1:4, Sensitivity_P, 'o-', 'DisplayName', 'P_{i}'); hold on;
plot(1:4, Sensitivity_G, 's-', 'DisplayName', 'G_{i}'); 
plot(1:4, Sensitivity_F, 'd-', 'DisplayName', 'F_{i}'); 

hold off; 

xticks(1:4);
xticklabels(categories);

xlabel('Life Stage');
ylabel('Sensitivity');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
grid on
legend;

%Sensitivity w.r.t. stage-specific survivorship rate p(i)
syms p_i real
syms d_i real

gamma_i = ((p_i/l)^d_i - (p_i/l)^(d_i-1)) / ((p_i/l)^d_i - 1);
G_i = p_i * gamma_i;
P_i = p_i * (1-gamma_i);
Diff_G_p = diff(G_i, p_i);
PD_G_p_Values = zeros(1,4);
p_i_values=[surv1 surv2 surv3 surv4];
d_i_values=[1 1 1 13];
for i=1:3
    PD_G_p = subs(Diff_G_p, {p_i,d_i}, {p_i_values(i), d_i_values(i)});
    PD_G_p_Values(i) = double(PD_G_p);
end
PD_G_p_Values(4)=0;

Diff_P_p = diff(P_i, p_i);
PD_P_p_Values = zeros(1,4);
for i=1:4
    PD_P_p = subs(Diff_P_p, {p_i,d_i}, {p_i_values(i), d_i_values(i)});
    PD_P_p_Values(i) = double(PD_P_p);
end

Sensitivity_p=zeros(1,4);
for i=1:3
    Sensitivity_p(i)=((p(i)/lambda)*((v(i)*w(i))/dot(v,w))*PD_P_p_Values(i))+((p(i)/lambda)*((v(i+1)*w(i))/dot(v,w))*PD_G_p_Values(i));
end
Sensitivity_p(4)=((p(4)/lambda)*((v(4)*w(4))/dot(v,w))*PD_P_p_Values(4));

%Sensitivity w.r.t. stage duration i.e. d_i
Diff_G_d = diff(G_i, d_i);
PD_G_d_Values = zeros(1,4);
for i=1:3
    PD_G_d = subs(Diff_G_d, {p_i,d_i}, {p_i_values(i), d_i_values(i)});
    PD_G_d_Values(i) = double(PD_G_d);
end
PD_G_d_Values(4)=0;

Diff_P_d = diff(P_i, d_i);
PD_P_d_Values = zeros(1,4);
for i=1:4
    PD_P_d = subs(Diff_P_d, {p_i,d_i}, {p_i_values(i), d_i_values(i)});
    PD_P_d_Values(i) = double(PD_P_d);
end

Sensitivity_d=zeros(1,4);
for i=1:3
    Sensitivity_d(i)=((d(i)/lambda)*((v(i)*w(i))/dot(v,w))*PD_P_d_Values(i))+((d(i)/lambda)*((v(i+1)*w(i))/dot(v,w))*PD_G_d_Values(i));
end
Sensitivity_d(4)=((d(4)/lambda)*((v(4)*w(4))/dot(v,w))*PD_P_d_Values(4));

% Plotting the sensitivities on the same graph
figure;
plot(1:4, Sensitivity_p, 's-', 'DisplayName','p_{i}'); hold on
plot(1:4, Sensitivity_d, 'o--', 'DisplayName','d_{i}');

hold off
xticks(1:4);
xticklabels(categories)
xlabel('Life Stage');
ylabel('Sensitivity');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
grid on
legend;

% How lamda changes with cub survival
p=[surv1 surv2 surv3 surv4];
p1_range = 0:0.05:1.0;
lambda_values_p1 = zeros(size(p1_range));
p1_values = p1_range;

for idx = 1:length(p1_range)
    p(1) = p1_range(idx); % Update p(1) for this iteration

    % Recalculate gamma, P and G based on the new p(1)
    gamma = zeros(1,4);
    for i = 1:4
        gamma(i) = (((p(i)/l)^d(i)) - ((p(i)/l)^(d(i)-1))) / (((p(i)/l)^d(i)) - 1);
    end
    P = [p(1)*(1-gamma(1)), p(2)*(1-gamma(2)), p(3)*(1-gamma(3)), p(4)*(1-gamma(4))];
    G = [p(1)*gamma(1), p(2)*gamma(2), p(3)*gamma(3)];

    % Update the population projection matrix A
    A=[P(1) F(2) F(3) F(4); G(1) P(2) 0 0; 0 G(2) P(3) 0; 0 0 G(3) P(4)];

    % Recalculate the intrinsic growth rate lambda
    [V,D,W] = eig(A);
    lambda_values_p1(idx) = max(diag(D));
end

figure; 
plot(p1_values, lambda_values_p1, '-o');
xlabel('Cub survival');
ylabel('Intrinsic growth rate (\lambda)');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
grid on;

%How lambda changes with juvenile survival
p=[surv1 surv2 surv3 surv4];
p2_range = 0:0.05:1.0;
lambda_values_p2 = zeros(size(p2_range));
p2_values = p2_range;

for idx = 1:length(p2_range)
    p(2) = p2_range(idx); 

    % Recalculate gamma, P and G based on the new p(2)
    gamma = zeros(1,4);
    for i = 1:4
        gamma(i) = (((p(i)/l)^d(i)) - ((p(i)/l)^(d(i)-1))) / (((p(i)/l)^d(i)) - 1);
    end
    P = [p(1)*(1-gamma(1)), p(2)*(1-gamma(2)), p(3)*(1-gamma(3)), p(4)*(1-gamma(4))];
    G = [p(1)*gamma(1), p(2)*gamma(2), p(3)*gamma(3)];

    % Update the population projection matrix A
    A=[P(1) F(2) F(3) F(4); G(1) P(2) 0 0; 0 G(2) P(3) 0; 0 0 G(3) P(4)];

    % Recalculate the intrinsic growth rate lambda
    [V,D,W] = eig(A);
    lambda_values_p2(idx) = max(diag(D));
end

figure; 
plot(p2_values, lambda_values_p2, '-o');
xlabel('Juvenile survival');
ylabel('Intrinsic growth rate (\lambda)')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
grid on;

%How lambda changes with subadult survival
p=[surv1 surv2 surv3 surv4];
p3_range = 0:0.05:1.0;
lambda_values_p3 = zeros(size(p3_range));
p3_values = p3_range;

for idx = 1:length(p3_range)
    p(3) = p3_range(idx); 

    % Recalculate gamma, P and G based on the new p(3)
    gamma = zeros(1,4);
    for i = 1:4
        gamma(i) = (((p(i)/l)^d(i)) - ((p(i)/l)^(d(i)-1))) / (((p(i)/l)^d(i)) - 1);
    end
    P = [p(1)*(1-gamma(1)), p(2)*(1-gamma(2)), p(3)*(1-gamma(3)), p(4)*(1-gamma(4))];
    G = [p(1)*gamma(1), p(2)*gamma(2), p(3)*gamma(3)];

    % Update the population projection matrix A
    A=[P(1) F(2) F(3) F(4); G(1) P(2) 0 0; 0 G(2) P(3) 0; 0 0 G(3) P(4)];

    % Recalculate the intrinsic growth rate lambda
    [V,D,W] = eig(A);
    lambda_values_p3(idx) = max(diag(D));
end

figure; 
plot(p3_values, lambda_values_p3, '-o');
xlabel('Subadult survival');
ylabel('Intrinsic growth rate (\lambda)')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
grid on;

% How lambda changes with adult survival
p=[surv1 surv2 surv3 surv4];
p4_range = 0:0.05:1.0;
lambda_values_p4 = zeros(size(p4_range));
p4_values = p4_range;

for idx = 1:length(p4_range)
    p(4) = p4_range(idx); 

    % Recalculate gamma, P and G based on the new p(4)
    gamma = zeros(1,4);
    for i = 1:4
        gamma(i) = (((p(i)/l)^d(i)) - ((p(i)/l)^(d(i)-1))) / (((p(i)/l)^d(i)) - 1);
    end
    P = [p(1)*(1-gamma(1)), p(2)*(1-gamma(2)), p(3)*(1-gamma(3)), p(4)*(1-gamma(4))];
    G = [p(1)*gamma(1), p(2)*gamma(2), p(3)*gamma(3)];

    % Update the population projection matrix A
    A=[P(1) F(2) F(3) F(4); G(1) P(2) 0 0; 0 G(2) P(3) 0; 0 0 G(3) P(4)];

    % Recalculate the intrinsic growth rate lambda
    [V,D,W] = eig(A);
    lambda_values_p4(idx) = max(diag(D));
end

figure; 
plot(p4_values, lambda_values_p4, '-o');
xlabel('Adult survival');
ylabel('Intrinsic growth rate (\lambda)')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
grid on;

% ERROR ANLYSIS
% Del(lambda)/del(a_ij)
Diff_lambda_P1 = ((v(1)*w(1))/dot(v,w));
Diff_lambda_P2 = ((v(2)*w(2))/dot(v,w));
Diff_lambda_P3 = ((v(3)*w(3))/dot(v,w));
Diff_lambda_P4 = ((v(4)*w(4))/dot(v,w));
Diff_lambda_G1 = ((v(2)*w(1))/dot(v,w));
Diff_lambda_G2 = ((v(3)*w(2))/dot(v,w));
Diff_lambda_G3 = ((v(4)*w(3))/dot(v,w));
Diff_lambda_F4 = ((v(1)*w(4))/dot(v,w));

LS= 2.39;
BI= 1.37;

% d(a_ij)
dP1=abs(PD_P_p_Values(1)*0.03) + abs(PD_P_d_Values(1)*0);
dP2=abs(PD_P_p_Values(2)*0.04) + abs(PD_P_d_Values(2)*0);
dP3=abs(PD_P_p_Values(3)*0.04) + abs(PD_P_d_Values(3)*0);
dP4=abs(PD_P_p_Values(4)*0.12) + abs(PD_P_d_Values(4)*0);
dG1=abs(PD_G_p_Values(1)*0.03) + abs(PD_G_d_Values(1)*0);
dG2=abs(PD_G_p_Values(2)*0.04) + abs(PD_G_d_Values(2)*0);
dG3=abs(PD_G_p_Values(3)*0.04) + abs(PD_G_d_Values(3)*0);
dF4=abs((1/BI)*0.12)+abs(((-1*LS)/BI^2)*0.25);

% Therefore, total error in lambda, given by del(lambda):
d_lambda = abs(Diff_lambda_P1*dP1) + abs(Diff_lambda_P2*dP2) + abs(Diff_lambda_P3*dP3) + abs(Diff_lambda_P4*dP4) + abs(Diff_lambda_G1*dG1) + abs(Diff_lambda_G2*dG2) + abs(Diff_lambda_G3*dG3) + abs(Diff_lambda_F4*dF4);
fprintf('Error estimate in lambda: \n')
disp(d_lambda)