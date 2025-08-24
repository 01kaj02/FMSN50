%% Home Assignment 1
% Victoria Lagerstedt & Kajsa Hansson Willis

clear;
clc;

load('powercurve_V164.mat');


%% Power production of a wind turbine
% 2a)
months = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"];

k = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];
lambda = [10.6 9.7 9.2 8.0 7.8 8.1 7.8 8.1 9.1 9.9 10.6 10.6];
l = length(k);

N = 10000;
Nroot = sqrt(N);
WD = zeros(N,l);
pValue = 1.96; % for 95% two-sided intervals 

mikro = zeros(1,l);
conf = zeros(2,l);

for i = 1:l
    WD(:,i) = wblrnd(lambda(i), k(i), N, 1);  
    mikro(:,i) = mean(P(WD(:,i)));
    stdev = std(P(WD(:,i)))/Nroot;
    conf(:,i) = [mikro(i)-pValue*stdev, mikro(i) + pValue*stdev];
    disp(months(i));
    fprintf('The confidence interval using draws from a Weibull distribution: %.4f [%.4f, %.4f]\n', mikro(i), conf(1,i), conf(2,i));
end


%%  The truncated version considered in problem 1a

%Instead of using the built in wblrnd here we use random variables from a
%uniform distribution and then transform them to a Weibull distribution. 

WDtrun = zeros(N,l);
mikro_trunc = zeros(1,l);
conf_trunc = zeros(2,l);

for i = 1:l
    x = unifrnd(0, 1, [N, 1]); 

    Fx = 1-exp(-(x/lambda(i)).^k(i)); 
    Fa = 1-exp(-(3.5/lambda(i)).^k(i));
    Fb = 1-exp(-(25/lambda(i)).^k(i)); 
    Fdif = Fb - Fa;


    WDTrun = lambda(i)*(-log(1-Fa-((Fb-Fa)*x))).^(1/k(i));
    WDtrun(:,i) = WDTrun;

    mikro_trunc(i) = mean(Fdif*P(WDtrun(:,i)));
    std_trunc = std(Fdif*P(WDtrun(:,i)))/Nroot; 
    conf_trunc(:,i) = [mikro_trunc(:,i)-pValue*std_trunc, mikro_trunc(:,i) + pValue*std_trunc];
    disp(months(i));
    fprintf('The confidence interval using the truncated version: %.4f [%.4f, %.4f]\n', mikro_trunc(i), conf_trunc(1, i), conf_trunc(2, i));
end


%% 2b) Wind as control variate 

mikroPCV = zeros(1,l);
confPCV = zeros(2,l);

for i=1:l
    m = gamma(1+(1/k(i)))*lambda(i); % The theoretical values 
    PCV1 = P(WD(:,i));

    COV = cov(WD(:,i), PCV1);
    c = - COV(1,2)/var(WD(:,1)); 

    PCV2 = c*(WD(:,i) - m);
    PCV = PCV1 + PCV2;
    mikroPCV(i) = mean(PCV);

    varPCV = var(PCV);
    stdPCV = sqrt(varPCV)/Nroot;

    confPCV(:,i) = [mikroPCV(i)-pValue*stdPCV, mikroPCV(i)+pValue*stdPCV];
    disp(months(i))
    fprintf('The confidence interval using wind as a control variate: %.4f\n [%.4f, %.4f]\n', mikroPCV(i), confPCV(1,i), confPCV(2,i));
end 

%% 2c) Importance sampling

figure
hold on
v_test = linspace(3.5,25,N);
f = (k(1)/lambda(1)).* (v_test/lambda(1)).^(k(1)-1).*exp(-(v_test/lambda(1)).^k(1));
testing = P(v_test)'.*f;
plot(testing)
g = 5000000*gampdf(v_test,10,1.3);
plot(g)
legend('Power times f','Instrumental distribution');
hold off

figure;
plot(testing./g);
legend("Should be constant") % Close enough! 

g = @(x) gampdf(x, 10,1.3); % The scaling factor is no longer necessary  

tau = zeros(1,l);
conf_imp = zeros(2,l);

for i = 1:l
    f = @(x) wblpdf(x, lambda(i), k(i));
    phi = @(x) P(x);
    X = gamrnd(10, 1.3, [1, N]); 
    target_func = @(x) (f(x) .* phi(x)') ./ g(x); 
    Y = target_func(X); 
    tau(i) = mean(Y);
    sd = std(Y);
    conf_imp(:,i) = [tau(i) - pValue * sd / Nroot, tau(i) + pValue * sd / Nroot];
    disp(months(i));
    fprintf('The confidence interval using importance sampling: %.4f [%.4f, %.4f]\n', tau(i), conf_imp(1, i), conf_imp(2, i));
end

%% 2d) Antithetic sampling
mean_d = zeros(1,l);
confAT = zeros(2,l);
WDtrun = zeros(N,l);
WDtrun2 = zeros(N,l);
Prob = zeros(1,l);
v_2d = linspace(3.5,25,N);
maxEffect = max(P(v_2d));
for i=1:l
    x = unifrnd(0,1,N/4,1); 
    Fx = 1-exp(-(x/lambda(i)).^k(i)); 
    Fa = 1-exp(-(3.5/lambda(i)).^k(i));
    Fb = 1-exp(-(14/lambda(i)).^k(i)); 
    Fdif = Fb - Fa;
    WDTrun = lambda(i)*(-log(1-Fa-((Fb-Fa)*x))).^(1/k(i));
    x1 = WDTrun;
    x2 = 1-x; 
    Fx2 = 1-exp(-(x2/lambda(i)).^k(i)); 
    Fa2 = 1-exp(-(3.5/lambda(i)).^k(i));
    Fb2 = 1-exp(-(14/lambda(i)).^k(i)); 
    Fdif2 = Fb2 - Fa2;
    WDTrun2 = lambda(i)*(-log(1-Fa2-((Fb2-Fa2)*x2))).^(1/k(i));
    x12 = WDTrun2;
    mean_d1 = mean(P(x1)+P(x12))/2;
    var_d1 = (1/4)*(var(P(x1))+var(P(x12)))+(1/2)*cov(P(x1),P(x12));
    std_d1 = sqrt(var_d1(1:1)/(N/2));
    FaP = 1-exp(-(14/lambda(i)).^k(i));
    FbP = 1-exp(-(25/lambda(i)).^k(i)); 
    Prob2(i) = FbP - FaP; % F(25)-F(14)
    Prob1(i) = FaP - Fa2; % F(14)-F(3.5)
    mean_d2 = Prob2(i)*maxEffect;
    var_d2 = 0;
    std_d2 = 0;
    mean_d(i) = mean_d2 + (mean_d1*Prob1(i)); 
    %std_d = (std_d1 + std_d2)/2;
    std_d = (std_d1 + std_d2);
    confAT(:,i) = [mean_d(i)-pValue*std_d/Nroot, mean_d(i) + pValue*std_d/Nroot];
    disp(months(i));
    fprintf('The confidence interval using antithetic sampling: %.4f [%.4f, %.4f]\n', mean_d(i), confAT(1,i), confAT(2,i));
end 
figure;
histogram(x1);
figure;
histogram(x12);
fprintf('The correlation for december: %.4f\n', corr(x1,x12));
%We only need half of the samples as we need in normal monte carlo. 

%% 2e) Probability

% The probability P(P(V)>0)is with the constraints given in the assignment = P(3.5≤V≤25)
PZeroPlus = zeros(1,l);

for i = 1:l
    Fa = 1-exp(-(3.5/lambda(i)).^k(i));
    Fb = 1-exp(-(25/lambda(i)).^k(i)); 

    PZeroPlus(i) = Fb - Fa; %P(P(V)>0) = P(3.5≤V≤25)
    disp(months(i));
    fprintf('Probability (percentage) that the turbine delivers power: %.4f\n', PZeroPlus(i)*100);
end

%% 2f) Confidence interval for the average ratio of actual wind turbine output to total windpower

constant = 0.5*1.225*pi*(164^2)/4;

for i = 1:l
    mean_V = gamma(1+(3/k(i)))*lambda(i)^3;
    Ptot_E = constant*mean_V; 
        
    P_E = mean(P(WD(:,i)));
    std_P = std(P(WD(:,i)));

    ratio = P_E/Ptot_E;

    % Since there is no variability in Ptot the standard deviation primarily
    % comes from the sampling in the P_E function. With large sample sizes the
    % central limit theorem holds so we can assume normal distribution. 
    stdR = (std_P / sqrt(N)) / Ptot_E; 

    confI = [ratio - pValue*stdR, ratio + pValue*stdR];
    disp(months(i));
    fprintf('Estimation of the average ratio of actual wind turbine output to total wind power (average power coefficient) with confidence interval: %.4f\n[%.4f, %.4f]\n', ratio, confI(1),confI(2));

end 

%% 2g) Capacity factor and availability factor


t = 12; % Time period is monthly 

% Capacity factor
cf_a = sum(mikro)/((9.5*10^6)*t);
cf_a2 = sum(mikro_trunc)/(9.5*10^6*t);
cf_b = sum(mikroPCV)/(9.5*10^6*t);
cf_c = sum(tau)/(9.5*10^6*t);
cf_d = sum(mean_d)/(9.5*10^6*t);

fprintf('Capacity factor for part a with the standard Monte Carlo: %.4f\n', cf_a);
fprintf('Capacity factor for part a with the truncated Monte Carlo: %.4f\n', cf_a2);
fprintf('Capacity factor for part a with wind as a control variate: %.4f\n', cf_b);
fprintf('Capacity factor for part a with importance sampling: %.4f\n', cf_c);
fprintf('Capacity factor for part a with antithetic sampling: %.4f\n', cf_d);

% Availability factor
af = sum(PZeroPlus)/t;
fprintf('Availability factor: %.4f\n', af);

%% 3. Combined power production of two wind turbines

% 3a) 

k3 = 1.96;
lambda3 = 9.13;

% As we need to estimate the expected amount of combined power generated by
% both turbines E(P(V1)+P(V2)) we can exploit that they have the same
% distribution, which means the expected value is the same for both and it
% is equivalent to 2*E(P(V)). This means we can use the importance sampling
% from 2c as a foundation and scale it with a reasonable factor. 

scalingFactor = 2;
v_test = linspace(3.5,25,N);
g_plot = scalingFactor*5000000*gampdf(v_test, 10,1.3); % Same as in 2c multiplied by scaling factor 

figure
hold on
f = (k3/lambda3).* (v_test/lambda3).^(k3-1).*exp(-(v_test/lambda3).^k3);
testing = 2*P(v_test)'.*f; 
plot(testing)
plot(g_plot)
legend('Power times f','Instrumental distribution');
hold off

figure;
plot(testing./g_plot);
legend("Should be constant") % Not perfect but satisfactory

g = @(x) gampdf(x, 10,1.3);
phi = @(x) 2*P(x); % Since it is 2*P(V)
X = gamrnd(10, 1.3, [1, N]);
weights = f / g(X); 
target_func = @(x) weights .* phi(x)';
Y = target_func(X);
tau_3a = mean(Y);
sd_3a = std(Y);
conf_3a = [tau_3a - pValue*sd_3a/Nroot, tau_3a + pValue*sd_3a/Nroot];
fprintf('The confidence interval using importance sampling: %.4f [%.4f, %.4f]\n', tau_3a, conf_3a(1), conf_3a(2));



%% 3b) Covariance
% The Covariance C(P(V1),P(V2))
% It holds that C(P(V1),P(V2))=E(P(V1)P(V2))-E(P(V1))E(P(V2))
% With importance sampling: E(P(V1)P(V2))= Eg(P(V1)*P(V2))f(V1,V2)/g(V1,V2)
% Determine E(P(V)) and V(P(V)) for the univariate case with korrekt k and
% lambda
N3b = 100;
g = @(x) gampdf(x, 10,1.3);
f = @(x) wblpdf(x, lambda3, k3);
phi = @(x) P(x);
X = gamrnd(10, 1.3, [1, N3b]);
target_func = @(x) (f(x) .* phi(x)') ./ g(x); 
Y = target_func(X); 
tau_3b = mean(Y);
var_3b = var(Y);
v1_test = linspace(3.5,25,N3b);
v2_test = linspace(3.5,25,N3b);
f = bivWei(v1_test,v2_test,k3,lambda3);
[X, Y] = meshgrid(v1_test, v2_test);
g = (10^17)*gampdf(v1_test,10,1.3)'*gampdf(v2_test,25,0.9);
figure;
surf(X, Y, f.*P(v1_test).*P(v2_test), 'FaceAlpha', 0.6);
hold on;
surf(X, Y, g, 'FaceAlpha', 0.6);
colormap jet;
colorbar;
xlabel('v1')
ylabel('v2')
hold off;
v1_new = wblrnd(lambda3, k3, N3b, 1);
v2_new = wblrnd(lambda3, k3, N3b, 1);
g = gamrnd(10,1.3,N3b,1)*gamrnd(25,0.9,N3b,1)';
f_3b = bivWei(v1_new',v2_new',k3,lambda3);
EXY = mean(mean((f_3b*(P(v1_new)*P(v2_new)'))./g));
EX = mean(P(v1_new));
covar = EXY-(EX^2);
fprintf('The covariance of power production in two wind turbines: %.4f\n', covar);

%% 3b without importance sampling

N3b2 = 10000; % Number of samples


w1 = wblrnd(lambda3, k3, [1 N3b2]);
w2 = wblrnd(lambda3, k3, [1 N3b2]);

% Compute power output for each turbine
P1 = P(w1); % Power output for turbine 1
P2 = P(w2); % Power output for turbine 2

% Compute expectations using Monte Carlo
E_P1 = mean(P1); % E[P(V1)]
E_P2 = mean(P2); % E[P(V2)]
E_P1P2 = mean(P1 .* P2); % E[P(V1) P(V2)]

% Compute covariance
cov_P1_P2 = E_P1P2 - (E_P1 * E_P2);

% Display result
fprintf('The covariance between power production in two wind turbines: %.4f\n', cov_P1_P2);


%% 3c) Variance

%V(P(V1)+P(V2)) = V(P(V1))+V(P(V2))+2*C(P(V1),P(V2))
%D(P(V1)+P(V2))=sqrt(V(P(V1)+P(V2))
N3c = 100000;
wind_3c = wblrnd(lambda3, k3, N3c, 1);
sing_var = var(P(wind_3c));
 
var_3c = 2*(sing_var) + 2*(covar); % Leveraging information from 3b
std_3c = sqrt(var_3c);
fprintf('The variability of power production in two wind turbines: %.4f\n', var_3c);
fprintf('The standard deviation of power production in two wind turbines: %.4f\n', std_3c);



%% 3d) Confidence interval
%Analysing which instrumental distribution to use
N3b = 500;
v1_test = linspace(0,50,N3b);
v2_test = linspace(0,50,N3b);
[X, Y] = meshgrid(v1_test, v2_test);
C = 9.5*10^6;
Z1 = P(v1_test);
Z2 = P(v2_test)';
P_XY = bivWeiNew(X,Y,k3,lambda3);
valid = (Z1 + Z2 > C);
prob = P_XY .* valid;
prob(prob == 0) = NaN;
g = 0.8*gampdf(v1_test,10,1.3)'*gampdf(v2_test,10,1.3);
figure;
surf(X, Y, prob,'EdgeColor', 'none', 'FaceAlpha', 0.5); 
hold on;
surf(X, Y, g,'EdgeColor', 'none', 'FaceAlpha', 0.5); 
shading interp; 
colormap(parula);
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Probability');
title('P(V1)+P(V2) > 9.5MW');
hold off;
N3b = 100000;
results = zeros(100,1);
for i=1:100
    g1 = gamrnd(lambda3, k3, N3b,1);
    g2 = gamrnd(lambda3, k3, N3b,1);
    g = gampdf(g1,10,1.3).*gampdf(g2,10,1.3);
    Z1 = P(g1);
    Z2 = P(g2);
    P_XY = bivWeiNew(g1,g2,k3,lambda3);
    weights = P_XY ./ g;
    valid = (Z1 + Z2 > C);
    P_rob = sum(weights .* valid) ./ sum(weights);
    results(i,1)=P_rob;
    
end
tau = mean(results);
sd = std(results);
   
pValue = 1.96;  
conf = [tau - pValue * sd, tau + pValue * sd];
fprintf('The confidence interval for P_rob_test: %.4f [%.4f, %.4f]\n', tau, conf(1), conf(2));


%% 3d) Opposite probability
N3b = 100;
v1_test = linspace(0,50,N3b);
v2_test = linspace(0,50,N3b);
[X, Y] = meshgrid(v1_test, v2_test);
C = 9.5*10^6;
Z1 = P(v1_test);
Z2 = P(v2_test)';
P_XY = bivWeiNew(X,Y,k3,lambda3);
valid = (Z1 + Z2 < C);
prob = P_XY .* valid;
prob(prob == 0) = NaN;
g2 = 1.4*wblpdf(v1_test,lambda3,k3)'*wblpdf(v2_test,lambda3,k3);
figure;
surf(X, Y, prob, 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
hold on;
mesh(X, Y, g2, 'EdgeColor', 'none', 'FaceAlpha', 0.5)
shading interp; 
colormap(parula);
colorbar;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Probability');
title('3D Probability Heat Map');
hold off;
N3b = 100000;
results2 = zeros(100,1);
for i=1:100
    v1 = wblrnd(lambda3, k3, N3b,1);
    v2 = wblrnd(lambda3, k3, N3b,1);
    g2 = wblpdf(v1,lambda3,k3)'*wblpdf(v2,lambda3,k3);
    Z1 = P(v1);
    Z2 = P(v2);
    P_XY2 = bivWeiNew(v1,v2,k3,lambda3);
    weights = P_XY2 ./ g2;
    valid = (Z1 + Z2 > C);
    P_rob2 = sum(weights .* valid)/sum(weights);
    results2(i,1)=P_rob2;
end
tau = mean(results2);
sd = std(results2);
    
pValue = 1.96;  
conf = [tau - pValue * sd, tau + pValue * sd];
fprintf('The confidence interval for P_rob_test: %.4f [%.4f, %.4f]\n', tau, conf(1), conf(2));

