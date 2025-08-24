%% Home assignment 3: 
% FMSN50 Monte Carlo and Empirical Methods
% Kajsa Hansson Willis & Victoria Lagerstedt

clear; 
close all;
load coal_mine_disasters.mat
load atlantic.txt

%% Part 1: Coal mine disasters—constructing a complex MCMC algorithm

figure;
stairs(tau, 1:length(tau));


%% Hybrid MCMC algorithm
t1 = 1658;
tend = 1980;
yeardiff = tend-t1;
d = 5; 
epsilon = 1e-10;

N = 10000;
burnin = 1000;
interval = yeardiff/d;
psi = 0.05; % Experiment with this
rho = 5; % Experiment with this
t = zeros(d+1,N+1);
t(1,1) = t1;
t(d+1,1) = tend;
for i=2:d
    t(i,1) = t1+(i-1)*interval;
end
ninit = zeros(d,1);
lambda = zeros(d,1);
for i=1:d % filling up the n-vector
    ninit(i) = nTau(tau,t(i),t(i+1));
end
lambdaVector = zeros(d,N);
thetaVector = zeros(1,N); 
nVector = zeros(d,N+1);
num = 0;
nVector(:,1) = ninit;
for j=1:N

    if j >= N
    disp("Reached max iterations")
    end
    theta = gamrnd(2, 1 / psi);
    for i=1:d % filling up the initial lambda-vector
        lambda(i) = gamrnd(2, 1 ./ theta); % lambda_i
    end
    X_theta = gamrnd(2*d+2, 1/(psi+sum(lambda)));
    lambdaVector(1,j) = gamrnd(nVector(1,j) + 2, 1 ./ (X_theta + (t(2,j) - t(1,j))));
    for D=2:d 
        lambda(D) = gamrnd(nVector(D,j) + 2, 1 ./ (X_theta + (t(D+1,j) - t(D,j))));   
        lambdaVector(D,j) = lambda(D);
        % MH-sampling:
        e = unifrnd(-rho,rho,d,1);
        e(1) = 0;
        tstar=t;
        tstar(1:d,j)=tstar(1:d,j)+e;
        if all(diff(tstar(:,j)) > epsilon) 
            nnew = arrayfun(@(i) nTau(tau, tstar(i,j), tstar(i+1,j)), 1:d)';
            n = arrayfun(@(i) nTau(tau, t(i,j), t(i+1,j)), 1:d)';
            prob = alpha(d, lambda, e, n, nnew,t(:,j),tstar);
            test2(D,j) = t(D+1,j) - t(D,j);
            
            if rand <= prob
                t(:,j+1) = tstar(:,j);
                nVector(:,j+1) = nnew; % Updating disaster counts
                num = num +1;
            else
                nVector(:,j+1) = nVector(:,j); % If rejection
                t(:,j+1) = t(:,j);
            end
        else
            nVector(:,j+1) = nVector(:,j); % If rejection
            t(:,j+1) = t(:,j);
        end
    end
    thetaVector(j) = X_theta;
end


% Removing the burnin part
thetaVectorTrue = thetaVector(:,burnin+1:end);
lambdaVectorTrue = lambdaVector(:,burnin+1:end);
tVector = t(:,burnin+1:end);
nVectorTrue = nVector(:,burnin+1:end);
accept_rate = num/(N*(d-1)); % number of accepted ones divided by possible ones
fprintf('Metropolis-Hastings Acceptance Rate: %.2f%%\n', accept_rate * 100);
figure;
plot(thetaVectorTrue(1,:));
title('Trace Plot of θ');
xlabel('Iteration');
ylabel('θ');
figure;
hold on
for i=1:d
    plot(lambdaVectorTrue(i,:));
end
hold off
title('Trace plot of lambda');
xlabel('Iteration');
ylabel('lambda');
legend('first breakpoint', 'second breakpoint','third breakpoint','fourth breakpoint') 
figure;
hold on
for i=1:d+1
    plot(tVector(i,:));
end
hold off
title('Trace plot of the breakpoints');
xlabel('Iteration');
ylabel('Breakpoints');



%% Part 2: Parametric bootstrap for the 100-year Atlantic wave

load atlantic.txt
figure;
plot(atlantic);

%% 2b two-sided confidence intervals 

[beta,mu] = est_gumbel(atlantic);
inv = @(u,muV, betaV) -betaV*log(-log(u))+muV;
n = length(atlantic);
B = 2000;
beta_hat = beta;
mu_hat = mu;
beta_boot = zeros(1,B);
mu_boot = zeros(1,B);

for b = 1:B % bootstrap
    %ystar =  inv(rand(n, 1), mu_hat, beta_hat);
    I = randsample(n,n,true);
    ystar = atlantic(I);
    [beta_boot(b),  mu_boot(b)] = est_gumbel(ystar);
end
mu_delta = sort(mu_boot - mu_hat); % sorting to obtain quantiles
beta_delta = sort(beta_boot - beta_hat); % sorting to obtain quantiles
alpha = 0.01; % CB level
mu_L = mu_hat - mu_delta(ceil((1 - alpha/2)*B));
mu_U = mu_hat - mu_delta(ceil(alpha*B/2));
beta_L = beta_hat - beta_delta(ceil((1 - alpha/2)*B));
beta_U = beta_hat - beta_delta(ceil(alpha*B/2));
fprintf('99%% Confidence Interval for mu: [%.3f, %.3f]\n', mu_L, mu_U);
fprintf('99%% Confidence Interval for beta: [%.3f, %.3f]\n', beta_L, beta_U);


figure;
subplot(2,1,1);
histogram(mu_boot, 50, 'Normalization', 'pdf', 'FaceColor', 'b');
hold on;
y = ylim;
plot([mu_L mu_L], y, 'g--', 'LineWidth', 2);
plot([mu_U mu_U], y, 'm--', 'LineWidth', 2);
plot([mu_hat mu_hat], y, 'r--', 'LineWidth', 2);
hold off;
xlabel('\mu Estimates');
ylabel('Density');
title('Bootstrap Distribution of \mu');
legend('Bootstrap Estimates', 'Lower Bound', 'Upper Bound', 'Original Estimate');
subplot(2,1,2);
histogram(beta_boot, 50, 'Normalization', 'pdf', 'FaceColor', 'b');
hold on;
y = ylim;
plot([beta_L beta_L], y, 'g--', 'LineWidth', 2);
plot([beta_U beta_U], y, 'm--', 'LineWidth', 2);
plot([beta_hat beta_hat], y, 'r--', 'LineWidth', 2);
hold off;
xlabel('\beta Estimates');
ylabel('Density');
title('Bootstrap Distribution of \beta');
legend('Bootstrap Estimates', 'Lower Bound', 'Upper Bound', 'Original Estimate');


%% 2c one-sided confidence interval 

T = 3*14*100; 
val = inv(1-1/T,mu, beta); 
disp(val);
est = zeros(1, B);
for b = 1:B
    est(b) = inv(1 - 1 / T, mu_boot(b), beta_boot(b));
end
delta = sort(est-val);
U = val - delta(ceil((alpha) * B));
disp(U)
