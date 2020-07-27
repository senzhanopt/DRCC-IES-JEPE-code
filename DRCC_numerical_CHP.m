%% parameters: parameter units in $USD, kg, MW, hour, celsius degree, Volt
% single cell in kW for numerical concern
clear;
tic
load cluster;
load wind;

%system
discount = 0.05;
lifetime = 20;
CRF = discount*(1+discount)^lifetime/((1+discount)^lifetime-1);
epsilon = 0.05; %confidence level for DRCC
T = 24;  
Z = 2; %no. of WTs
windCap = 18; %MW, single wind farm
k_r = transpose(weights);  %weights of representative days
n_r = length(k_r); %representative days
m = zeros(Z, n_r, T);
for r = 1:n_r
    for t = 1:T % Z * n_r * T, predicted wind generation vector
        m(:,r,t) = windCap * transpose(windOutput_mean(:,:,r,t));
    end
end
Sigma = windCap^2 .* Sigma; % covariance (Z * Z) * (n_r * T), wind prediction error covariance
Sigma_coupling = windCap^2 .* Sigma_coupling; % (2Z * 2Z) * (n_r * T) with first being 0
N = 5; %P2HH stack <= (2^N - 1)*100*1.99E-3 =6.19;
M = 10; %big M

%electrolyser 
c_el = CRF * 750 * 1.13 * 2.03 + 93.8*2.03;  %$/cell:single cell 2.03 kW , 750eur/kw                                    
c_H2 = 1.13 * 5; %$/kg 
P_H2 =  [0.745,1.491, 1.497, 0.749 ]; %P2HH operation corner points, kW
P_Heat = [0.184,0.539, 0.472, 0.164 ];%kW
Temp_P2HH = [60, 60, 80, 80];%Celsius degree
n_P2HH_corner = length(P_H2);
U_tn = 1.494;%V
F = 96485;
Temp_init = 80;%Celsius degree
C = (2.5*0.66*4200+2500*0.66*0.3/56*65.87)/(3.6E6); %specific heat, kW/Celsius degree
T_min = 60;
T_max = 80;
R_eqv = 0.167E3; %C / kW
Ta = 25;
c_converter = 155E3 * discount*(1+discount)^10/((1+discount)^10-1); % $/MW
eta_converter = 0.95;

%CHP
SU = 14; 
RU = 14; %MW
RD = 14;
SD = 14;  %MW
c_su = 50; %start-up cost     $                           
c_sd = 50; %shut-down       $
P = [5, 3.97, 13.9, 15.5]; %MW
Q = [0, 5.14, 10.7, 0]; %MW
Cost = 20 * 1.13 .* (0.5*Q+P); % CHP corner point cost, natural
n_CHP_corner = length(Cost);
min_on = 2; %min ON time, hour
min_off = 2; %min OFF time, hour

%power grid
c_trans = 5;    %$/(MW*MW)                                    
% H_CHP; % R (L * U)
% L = length(H_CHP);
% H_wind;
% H_P2HH;
% H_trans;
% H_D;
% fmax; %MW

%electric boiler
c_EB = CRF * 178.84E3;
eta_EB = 0.9; %power to heat ratio

%demand
d_p = 42  .* clustered_e_dh_demands(:,1:T);  %MW
d_q = 11.5 .* clustered_e_dh_demands(:,T+1:2*T); %MW

%STORAGE
eta_comp = 1.1E-3; % MWh/kg
c_comp = 12E3 * CRF; % per kg/h
c_tank = 822 * CRF; %per kg
m_tank_max = 1500;
cputime_start = cputime;
%% define variables
u_su = binvar(n_r, T, 'full'); %start-up
u_sd = binvar(n_r, T, 'full'); %shut-down
u = binvar(n_r, T, 'full'); %ON/OFF
x = sdpvar(n_CHP_corner, n_r, T, 'full'); %CHP corner point weights
alpha = sdpvar(n_CHP_corner, n_r, T, 'full');
p_trans = sdpvar(n_r, T, 'full'); 
beta = sdpvar(n_r, T, 'full');%trans 
C_times_x = sdpvar(n_r, T, 'full'); %CHP cost 
%% objective
Sigma_sum = ones(n_r, T); % 1^{\top} * Sigma_{r,t} * 1
for i = 1:n_r
    for j = 1:T
        Sigma_sum(i,j) = sum(sum(Sigma(:,:,i,j)));
    end
end
obj =  c_su * (k_r * sum(u_su, 2)) +  c_sd * (k_r * sum(u_sd, 2)) ...
     + k_r * sum(C_times_x, 2)...
    + c_trans * (k_r * sum(p_trans.*p_trans + Sigma_sum .* (beta.*beta), 2));
%% constraints
cons = [ ];
for i = 1:n_r
    for j = 1:T
        cons = [cons, C_times_x(i,j) == Cost * x(:,i,j)];
    end
end
% Unit commitment 
cons = [cons, u(:,1) - u_su(:,1) <= 0];    
for t = 2:T
    cons = [cons, -u(:,t-1) + u(:,t) - u_su(:,t) <= 0];
    cons = [cons, u(:,t-1) - u(:,t) - u_sd(:,t) <= 0];
end
for tao = 2: min_on
    cons = [cons, u(:,1) - u(:, tao) <= 0];
end
for t = 2:T-1
    for tao = t+1: min(T, min_on + t - 1)
        cons = [cons, -u(:,t-1) + u(:,t) - u(:, tao) <= 0];
    end
    for tao = t+1: min(T, min_off + t - 1)
        cons = [cons, u(:,t-1) - u(:,t) + u(:, tao) <= 1];
    end
end
% x
for r = 1:n_r
    for t = 1:T
        cons = [cons, ones(1,n_CHP_corner) * x(:,r,t) == u(r,t)];
        cons = [cons, ones(1,n_CHP_corner) * alpha(:,r,t) == 0];
        for k = 1:n_CHP_corner
                cons = [cons, cone(alpha(k, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                        sqrt(epsilon/(1-epsilon))*x(k, r, t))];
                cons = [cons, cone(alpha(k, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                        sqrt(epsilon/(1-epsilon))*(1-x(k, r, t)))];
        end
    end
end
% ramp up at t = 1
for r = 1:n_r
    t = 1;
    cons = [cons, cone(P * alpha(:, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                sqrt(epsilon/(1-epsilon))*(SU - P * x(:,r,t)))]; 
end
%ramp up at t >= 2
for r = 1:n_r
    for t = 2:T
        cons = [cons, cone( sqrtm(Sigma_coupling(:,:,r,t))*( P*alpha(:,r,t)* [zeros(Z,1); ones(Z,1)]...
           - P*alpha(:,r,t-1)*[ones(Z,1); zeros(Z,1)] ),sqrt(epsilon/(1-epsilon))*( P*x(:,r,t-1)...
           -P*x(:,r,t)+SU*(1-u(r,t-1)) + RU*u(r,t-1) ) )];
        cons = [cons, cone( sqrtm(Sigma_coupling(:,:,r,t))*( P*alpha(:,r,t)* [zeros(Z,1); ones(Z,1)]...
           - P*alpha(:,r,t-1)*[ones(Z,1); zeros(Z,1)] ),sqrt(epsilon/(1-epsilon))*( - P*x(:,r,t-1)...
           + P*x(:,r,t)+SD*(1-u(r,t)) + RD*u(r,t) ) )];
    end
end
% power, heat balance
for r = 1:n_r
    for t = 1:T
        cons = [cons, p_trans(r,t) + P*x(:,r,t) + ones(1,Z)*m(:,r,t) ...
            == d_p(r,t) ];
        cons = [cons, beta(r,t) + P*alpha(:,r,t)+1 == 0];
        cons = [cons, Q*x(:,r,t) == d_q(r,t)];
        cons = [cons, Q*alpha(:,r,t) == 0];
    end
end

%% solving model and result presentation
options = sdpsettings( 'solver', 'gurobi', 'gurobi.TimeLimit', 300);
result = optimize ( cons, obj, options);
cputime_used = cputime - cputime_start;
obj = value(obj);
u_su = value(u_su); %start-up
u_sd = value(u_sd); %shut-down
u = value(u); %ON/OFF
x = value(x); %CHP corner point weights
alpha = value(alpha);
p_trans = value(p_trans);
beta = value(beta);%trans 
C_times_x = value(C_times_x);
toc
save('DRCC_CHP.mat',"obj", "u_su", "u_sd", ...
    "u", "x", "alpha","p_trans","beta",...
    "C_times_x");


T=24;
n_r = 10;
P = [5, 3.97, 13.9, 15.5]; %MW
Q = [0, 5.14, 10.7, 0]; %MW
plot(Q,P);
xlabel("Q(MW)");
ylabel("P(MW)")
hold on;
for r = 1:n_r
    for t = 1:T
        P_DRCC(r,t) = P * x(:,r,t);
        Q_DRCC(r,t) = Q * x(:,r,t);
    end
end
scatter(Q_DRCC(:),P_DRCC(:));
hold off;
saveas(gcf, 'DRCC_CHP','epsc');
disp("inverse flow");
weights' * sum(abs(p_trans),2)
disp("obj");
obj






        


