%% Home assignment 2: Self-avoiding walks in Z^d & Filter estimation of noisy population measurements
% FMSN50 Monte Carlo and Empirical Methods
% Kajsa Hansson Willis & Victoria Lagerstedt

clear; 
close all;
load population_2024.mat

%% Startup

N = 1*10^5; 

n = 11; 
d = 2; 
powerLawCorrection = @(n) n^(43/32-1); % for d=2 

%% Part 1: Self-avoiding walks in Z^d

% Problem 3 

X_vector3 = zeros(n+1,N); 
Y_vector3 = zeros(n+1,N);

count_vector = zeros(1, n); 
est = zeros(1, n);

for i=1:N
    SAW = true;
    for p=1:n
        rnd = randi([1 4],1); 
        if rnd == 1
            X_vector3(p+1, i) = X_vector3(p,i) + 1; 
            Y_vector3(p+1, i) = Y_vector3(p,i);
        elseif rnd == 2
            X_vector3(p+1,i) = X_vector3(p,i) - 1;
            Y_vector3(p+1,i) = Y_vector3(p,i); 
        elseif rnd == 3
            X_vector3(p+1,i) = X_vector3(p,i);
            Y_vector3(p+1,i) = Y_vector3(p,i) + 1;
        elseif rnd == 4
            X_vector3(p+1,i) = X_vector3(p,i);
            Y_vector3(p+1,i) = Y_vector3(p,i) - 1;
        end 
         if any(X_vector3(1:p, i) == X_vector3(p+1, i) & Y_vector3(1:p, i) == Y_vector3(p+1, i))
            SAW = false;
            count_vector(p:end) = count_vector(p:end) + 1;
            break; 
         end
    end
end 

for p = 1:n
    Nsa = N - count_vector(p);
    cnEst = Nsa / N; 
    walksTot = 4^p; 
    est(p) = walksTot * cnEst; 
end

disp("Estimated c_n(2) values for each n:"); % Excluding 1 
disp(est);

fprintf("The estimated c%d(2) from problem 3 is: %.5f\n", n, est(n));



%% Problem 4

X_vector4 = zeros(n+1,N); 
Y_vector4 = zeros(n+1,N);
SAW_vector4 = zeros(1,4);
z_vector = zeros(n,N);
for i=1:N
    SAW = true;
    for p=1:n
        SAW_vector4 = zeros(1,4);
        if any(X_vector4(1:p, i) == (X_vector4(p, i)+1) & Y_vector4(1:p, i) == Y_vector4(p, i))
            SAW_vector4(1) = 1;
        end 
        if any(X_vector4(1:p, i) == (X_vector4(p, i)-1) & Y_vector4(1:p, i) == Y_vector4(p, i))
            SAW_vector4(2) = 1;
        end 
        if any(X_vector4(1:p, i) == (X_vector4(p, i)) & Y_vector4(1:p, i) == (Y_vector4(p, i)+1))
            SAW_vector4(3) = 1;
        end 
         if any(X_vector4(1:p, i) == (X_vector4(p, i)) & Y_vector4(1:p, i) == (Y_vector4(p, i)-1))
            SAW_vector4(4) = 1;
         end 
        noFreeNeighbors = 4-sum(SAW_vector4);
        z_vector(p,i) = noFreeNeighbors;
        if noFreeNeighbors == 0
            X_vector4(p,i) = X_vector4(p-1,i); 
            Y_vector4(p,i) = Y_vector4(p-1,i);
            break   
        end
        rnd = randi([1 4],1);  
        while SAW_vector4(rnd) == 1
            rnd = randi([1 4],1); 
        end
        if rnd == 1 
            X_vector4(p+1, i) = X_vector4(p,i) + 1; 
            Y_vector4(p+1, i) = Y_vector4(p,i);
        end 
        if rnd == 2
            X_vector4(p+1,i) = X_vector4(p,i) - 1;
            Y_vector4(p+1,i) = Y_vector4(p,i);
        end 
        if rnd == 3
            X_vector4(p+1,i) = X_vector4(p,i);
            Y_vector4(p+1,i) = Y_vector4(p,i) + 1;
        end
        if rnd == 4
            X_vector4(p+1,i) = X_vector4(p,i);
            Y_vector4(p+1,i) = Y_vector4(p,i) - 1;
        end
    end
end 

%SIS algorithm:
w_vectorSIS = zeros(n,N);
g = @(n,N) 1/z_vector(n-1,N);
z = 1; 
for i=1:N
    w_vectorSIS(1,i) = 1;
    for k = 1:n-1
        if z_vector(k+1,i) == 0 % indicator function part 
            w_vectorSIS(k+1,i) = 0;
        else
            w_vectorSIS(k+1,i) = (z/g(k+1,i))*w_vectorSIS(k,i); 
        end
    end
end
c_vector = mean(w_vectorSIS,2);
fprintf("The estimated c%d(2) from problem 4 is: %.5f\n", n, c_vector(n)); 
%% Problem 5 preparations (for problem 6) 

replicates = 10;
cn_tot = zeros(replicates,n+1); % This is filled with values by redoing Problem 5 ten times and changing the index 

%% Problem 5

X_vector5 = zeros(n+1,N); 
Y_vector5 = zeros(n+1,N);
SAW_vector5 = zeros(1,4);
z_vector5 = zeros(n+1,N);
X_vectorSIS5 = zeros(n+1,N);
w_vectorSIS5 = zeros(n+1,N);
g5 = @(n,N) 1/z_vector5(n-1,N);
z5 = 1;
c_vector5 = zeros(1,n); 
c_vector5(1,1) = 1;
estSAW5 = zeros(1,n);
for p=1:n
    for i=1:N
        w_vectorSIS5(1,i) = 1;
        SAW_vector5 = zeros(N,4);
        if any(X_vector5(1:p, i) == (X_vector5(p, i)+1) & Y_vector5(1:p, i) == Y_vector5(p, i))
            SAW_vector5(i,1) = 1;
        end 
        if any(X_vector5(1:p, i) == (X_vector5(p, i)-1) & Y_vector5(1:p, i) == Y_vector5(p, i))
            SAW_vector5(i,2) = 1;
        end 
        if any(X_vector5(1:p, i) == (X_vector5(p, i)) & Y_vector5(1:p, i) == (Y_vector5(p, i)+1))
            SAW_vector5(i,3) = 1;
        end 
        if any(X_vector5(1:p, i) == (X_vector5(p, i)) & Y_vector5(1:p, i) == (Y_vector5(p, i)-1))
            SAW_vector5(i,4) = 1;
        end 
        noFreeNeighbors5 = 4-sum(SAW_vector5(i,:));
        z_vector5(p,i) = noFreeNeighbors5;
        if noFreeNeighbors5 == 0 % indicator function part 
            X_vector5(p,i) = X_vector5(p-1,i); 
            Y_vector5(p,i) = Y_vector5(p-1,i);
        else
            rnd = randi([1 4],1);  
            while SAW_vector5(i,rnd) == 1
                rnd = randi([1 4],1); 
            end
            if rnd == 1 
                X_vector5(p+1, i) = X_vector5(p,i) + 1; 
                Y_vector5(p+1, i) = Y_vector5(p,i);
            end 
            if rnd == 2
                X_vector5(p+1,i) = X_vector5(p,i) - 1;
                Y_vector5(p+1,i) = Y_vector5(p,i);
            end 
            if rnd == 3
                X_vector5(p+1,i) = X_vector5(p,i);
                Y_vector5(p+1,i) = Y_vector5(p,i) + 1;
            end
            if rnd == 4
                X_vector5(p+1,i) = X_vector5(p,i);
                Y_vector5(p+1,i) = Y_vector5(p,i) - 1;
            end
        end
        if z_vector5(p,i) == 0 
           w_vectorSIS5(p+1,i) = 0;
        else
           X_vectorSIS5(p+1,i) = 1/z_vector5(p,i);
           w_vectorSIS5(p+1,i) = z5/(1/z_vector5(p,i));
        end
    end
    indices = randsample(1:N, N, true, w_vectorSIS5(p+1,:)); 
    X_vector5 = X_vector5(:, indices); 
    Y_vector5 = Y_vector5(:, indices);
    c_vector5(1,p+1) = c_vector5(1,p)*mean(w_vectorSIS5(p,:));
    disp(c_vector5(1,p));
    w_vectorSIS5(p+1,i) = 1;
end
cn_tot(9,:) = c_vector5; % this is adjusted so that cn_tot is filled with 10 values 

fprintf("The estimated c%d(2) from problem 5 is: %.5f\n", n-1, c_vector5(n+1));


%% Problem 6

params = zeros(replicates, 3);
cn_tot2 = cn_tot;
cn_tot2(:,1) = []; % removing the first column of ones  

for i=1:replicates
    y =  log(cn_tot2(i,:));
    n6 = length(cn_tot2(i,:));
   

    %y = log(c_vector5); %for just one values 
    %n6 = length(c_vector5);
    x1 = (1:n6)';
    x2 = log(x1); 

    B = y(:);
    A = [ones(n6,1), x1, x2]; 
    [X, int] = regress(B,A);
    params(i, :) = X';    
end
disp(rank(A)); % should be complete 

beta1 = params(:, 1);
beta2 = params(:, 2);
beta3 = params(:, 3);

Ad_tot = exp(beta1);
Ad = mean(Ad_tot);
mud_tot = exp(beta2);
mud = mean(mud_tot);
gammad_tot = 1+beta3;
gammad = mean(gammad_tot);

fprintf('A_2 = %.4f\n', Ad);
fprintf('μ_2 = %.4f\n', mud);
fprintf('γ_2 = %.4f\n', gammad); % should be 43/32(1,34375)

%% Problem 9
d = 3; % Dimension
X_vector9 = zeros(n+1,N); 
Y_vector9 = zeros(n+1,N);
Third_vector9 = zeros(n+1,N); % Accounting for the added dimension
SAW_vector9 = zeros(1,6);
z_vector9 = zeros(n+1,N);
w_vectorSIS9 = zeros(n+1,N);
g9 = @(n,N) 1/z_vector9(n-1,N);
z9 = 1;
c_vector9 = zeros(1,n); 
c_vector9(1,1) = 1;
estSAW9 = zeros(1,n); 
for p=1:n
    for i=1:N
        w_vectorSIS9(1,i) = 1;
        SAW_vector9 = zeros(N,6);
        if any(X_vector9(1:p, i) == (X_vector9(p, i)+1) & Y_vector9(1:p, i) == Y_vector9(p, i) & Third_vector9(1:p, i) == Third_vector9(p, i))
            SAW_vector9(i,1) = 1;
        end 
        if any(X_vector9(1:p, i) == (X_vector9(p, i)-1) & Y_vector9(1:p, i) == Y_vector9(p, i) & Third_vector9(1:p, i) == Third_vector9(p, i))
            SAW_vector9(i,2) = 1;
        end 
        if any(X_vector9(1:p, i) == (X_vector9(p, i)) & Y_vector9(1:p, i) == (Y_vector9(p, i)+1) & Third_vector9(1:p, i) == Third_vector9(p, i))
            SAW_vector9(i,3) = 1;
        end 
        if any(X_vector9(1:p, i) == (X_vector9(p, i)) & Y_vector9(1:p, i) == (Y_vector9(p, i)-1) & Third_vector9(1:p, i) == Third_vector9(p, i))
            SAW_vector9(i,4) = 1;
        end 
         if any(X_vector9(1:p, i) == (X_vector9(p, i)) & Y_vector9(1:p, i) == (Y_vector9(p, i)) & Third_vector9(1:p, i) == (Third_vector9(p, i)+1))
            SAW_vector9(i,5) = 1;
        end 
        if any(X_vector9(1:p, i) == (X_vector9(p, i)) & (Y_vector9(1:p, i) == (Y_vector9(p, i))) & (Third_vector9(1:p, i) == (Third_vector9(p, i)-1)))
            SAW_vector9(i,6) = 1;
        end 
        noFreeNeighbors9 = 6-sum(SAW_vector9(i,:));
        z_vector9(p,i) = noFreeNeighbors9;
        if noFreeNeighbors9 == 0 
            X_vector9(p,i) = X_vector9(p-1,i); 
            Y_vector9(p,i) = Y_vector9(p-1,i);
        else
            rnd = randi([1 6],1);  
            while SAW_vector9(i,rnd) == 1
                rnd = randi([1 6],1); 
            end
            if rnd == 1 
                X_vector9(p+1, i) = X_vector9(p,i) + 1; 
                Y_vector9(p+1, i) = Y_vector9(p,i);
                Third_vector9(p+1, i) = Third_vector9(p,i);
            end 
            if rnd == 2
                X_vector9(p+1,i) = X_vector9(p,i) - 1;
                Y_vector9(p+1,i) = Y_vector9(p,i);
                Third_vector9(p+1, i) = Third_vector9(p,i);
            end 
            if rnd == 3
                X_vector9(p+1,i) = X_vector9(p,i);
                Y_vector9(p+1,i) = Y_vector9(p,i) + 1;
                Third_vector9(p+1, i) = Third_vector9(p,i);
            end
            if rnd == 4
                X_vector9(p+1,i) = X_vector9(p,i);
                Y_vector9(p+1,i) = Y_vector9(p,i) - 1;
                Third_vector9(p+1, i) = Third_vector9(p,i);
            end
            if rnd == 5
                X_vector9(p+1,i) = X_vector9(p,i);
                Y_vector9(p+1,i) = Y_vector9(p,i);
                Third_vector9(p+1, i) = Third_vector9(p,i)+1;
            end
            if rnd == 6
                X_vector9(p+1,i) = X_vector9(p,i);
                Y_vector9(p+1,i) = Y_vector9(p,i);
                Third_vector9(p+1, i) = Third_vector9(p,i)-1;
            end
        end
        if z_vector9(p,i) == 0 
           w_vectorSIS9(p+1,i) = 0;
        else
           w_vectorSIS9(p+1,i) = z9/(1/z_vector9(p,i));
        end
    end
    indices = randsample(1:N, N, true, w_vectorSIS9(p+1,:)); 
    X_vector9 = X_vector9(:, indices); 
    Y_vector9 = Y_vector9(:, indices);
    Third_vector9 = Third_vector9(:, indices);
    c_vector9(1,p+1) = c_vector9(1,p)*mean(w_vectorSIS9(p,:));
    w_vectorSIS9(p+1,i) = 1; 
end
fprintf("The estimated c%d(3) is: %.5f\n", n-1, c_vector9(n+1));
y9 = log(c_vector9(2:end)); 
n9 = length(y9);
x19 = (1:n9)';  
x29 = log(x19);
B9 = y9(:);
A9 = [ones(n9,1), x19, x29]; 
[X9, int] = regress(B9,A9);
params9 = X9';    
disp(rank(A9));
beta19 = params9(:, 1);
beta29 = params9(:, 2);
beta39 = params9(:, 3);
Ad9 = exp(beta19);
mud9 = exp(beta29);
gammad9 = 1+beta39;
fprintf('A_3 = %.4f\n', Ad9); 
fprintf('μ_3 = %.4f\n', mud9); % Should be between 3 and 5 for d=3 
fprintf('γ_3 = %.4f\n', gammad9);


%% Problem 10
clear; 
load population_2024.mat

N = 1000;  
n = 100;
A = 0.8;
B = 3.8;
C = 0.6;
D = 0.99; 
G = 0.8;
H = 1.25;
tau = zeros(1,n+1); % vector of filter means
w = zeros(N,1);
p = @(x,y) unifpdf(y,G*x,H*x); % observation density 

for k = 1:n+1 
    if k == 1
        part = C + rand(N, 1)*(D-C); % Initialization U(0.6,0.99)
    else
        part = (A + rand(N, 1)*(B-A)).*part.*(1-part); %U(0.8,3.8)
    end
    
    w = p(part, Y(k)); 
    tau(k) = sum(part.*w) / sum(w); 
    [xx,I]=sort(part); 
    cw=cumsum(w(I))/sum(w);  

    Ilower=find(cw>=0.025,1); 
    Iupper=find(cw>=0.975,1); 
    taulower(k)=xx(Ilower); 
    tauupper(k)=xx(Iupper); 
    ind = randsample(N, N, true, w); 
    part = part(ind);
end
x = 0:n;

% Plot tau
figure()
hold on;
fill([x, fliplr(x)], [tauupper, fliplr(taulower)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(0:n, tau, 'b', 'DisplayName', 'Estimated size(X)');
plot(0:n, X, 'r', 'DisplayName', 'True size(X)'); 
xlabel('Generation (k)');
ylabel('Population size');
legend();
grid on;
hold off;

