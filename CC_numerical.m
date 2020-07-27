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
epsilon = 0.01; %confidence level for DRCC
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
n_el = intvar(1);
P_EB = sdpvar(1);
m_comp = sdpvar(1);
P_converter = intvar(1);
m_tank = sdpvar(1);
z = binvar(N, 1, 'full');
u_su = binvar(n_r, T, 'full'); %start-up
u_sd = binvar(n_r, T, 'full'); %shut-down
u = binvar(n_r, T, 'full'); %ON/OFF
x = sdpvar(n_CHP_corner, n_r, T, 'full'); %CHP corner point weights
alpha = sdpvar(n_CHP_corner, n_r, T, 'full');
p_trans = sdpvar(n_r, T, 'full'); 
beta = sdpvar(n_r, T, 'full');%trans 
p_EB = sdpvar(n_r, T, 'full'); 
rho = sdpvar(n_r, T, 'full');%EB p 
q_EB = sdpvar(n_r, T, 'full'); 
nu = sdpvar(n_r, T, 'full');%EB q
n_H2 = sdpvar(n_r, T, 'full'); % kg
gamma = sdpvar(n_r, T, 'full');
p_P2HH = sdpvar(n_r, T, 'full');%MW
Delta = sdpvar(n_r, T, 'full');
p_p2hh = sdpvar(n_r, T, 'full');
delta = sdpvar(n_r, T, 'full');
q_EXC = sdpvar(n_r, T, 'full');
Lambda = sdpvar(n_r, T, 'full');
q_exc = sdpvar(n_r, T, 'full');
lambda = sdpvar(n_r, T, 'full');
h_P2HH = sdpvar(n_r, T, 'full');
Pi = sdpvar(n_r, T, 'full');
h_p2hh = sdpvar(n_r, T, 'full');
pi = sdpvar(n_r, T, 'full');
q_p2hh = sdpvar(n_r, T, 'full');
kappa = sdpvar(n_r, T, 'full');
y = sdpvar(n_P2HH_corner, n_r, T, 'full');
ksi = sdpvar(n_P2HH_corner, n_r, T, 'full');
Temp = sdpvar(n_r, T, 'full'); 
mu = sdpvar(n_r, T, 'full');
e1 = sdpvar(N, n_r, T, 'full');
e2 = sdpvar(N, n_r, T, 'full');
e3 = sdpvar(N, n_r, T, 'full');
e4 = sdpvar(N, n_r, T, 'full');
e5 = sdpvar(N, n_r, T, 'full');
e6 = sdpvar(N, n_r, T, 'full');
C_times_x = sdpvar(n_r, T, 'full'); %CHP cost 
m_H2 = sdpvar(n_r, T, 'full');
upsilon = sdpvar(n_r, T, 'full');
%% objective
Sigma_sum = ones(n_r, T); % 1^{\top} * Sigma_{r,t} * 1
for i = 1:n_r
    for j = 1:T
        Sigma_sum(i,j) = sum(sum(Sigma(:,:,i,j)));
    end
end
obj = c_el*n_el + m_tank*c_tank + c_EB * P_EB + c_converter * P_converter...
    + c_comp * m_comp - c_H2 * (k_r * sum(n_H2, 2)) ...
    + c_su * (k_r * sum(u_su, 2)) +c_sd * (k_r * sum(u_sd, 2)) + k_r * sum(C_times_x, 2)...
    + c_trans * (k_r * sum(p_trans.*p_trans + Sigma_sum .* (beta.*beta), 2))...
    ;
%% constraints
cons = [ ];
cons = [cons, m_tank <= m_tank_max];
cons = [cons, m_tank >= 0, P_EB >= 0, m_comp >= 0, P_converter >= 0];
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
                        1/norminv(1-epsilon)*x(k, r, t))];
                cons = [cons, cone(alpha(k, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                        1/norminv(1-epsilon)*(1-x(k, r, t)))];
        end
    end
end
% ramp up at t = 1
for r = 1:n_r
    t = 1;
    cons = [cons, cone(P * alpha(:, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                1/norminv(1-epsilon)*(SU - P * x(:,r,t)))]; 
end
%ramp up at t >= 2
for r = 1:n_r
    for t = 2:T
        cons = [cons, cone( sqrtm(Sigma_coupling(:,:,r,t))*( P*alpha(:,r,t)* [zeros(Z,1); ones(Z,1)]...
           - P*alpha(:,r,t-1)*[ones(Z,1); zeros(Z,1)] ),1/norminv(1-epsilon)*( P*x(:,r,t-1)...
           -P*x(:,r,t)+SU*(1-u(r,t-1)) + RU*u(r,t-1) ) )];
        cons = [cons, cone( sqrtm(Sigma_coupling(:,:,r,t))*( P*alpha(:,r,t)* [zeros(Z,1); ones(Z,1)]...
           - P*alpha(:,r,t-1)*[ones(Z,1); zeros(Z,1)] ),1/norminv(1-epsilon)*( - P*x(:,r,t-1)...
           + P*x(:,r,t)+SD*(1-u(r,t)) + RD*u(r,t) ) )];
    end
end
% power, heat balance
for r = 1:n_r
    for t = 1:T
        cons = [cons, 1E2 * ( p_trans(r,t) + P*x(:,r,t) + ones(1,Z)*m(:,r,t) )...
            == 1E2 * ( p_P2HH(r,t)/eta_converter+p_EB(r,t) + d_p(r,t) + n_H2(r,t)*eta_comp )];
        cons = [cons, 1E2 * ( beta(r,t) + P*alpha(:,r,t)+1) == 1E2 * (Delta(r,t)/eta_converter + rho(r,t)+gamma(r,t)*eta_comp)];
        cons = [cons, Q*x(:,r,t) + q_EXC(r,t)+q_EB(r,t) == d_q(r,t)];
        cons = [cons, Q*alpha(:,r,t) + Lambda(r,t) + nu(r,t)== 0];
    end
end
% electric boiler
for r = 1:n_r
    for t = 1:T
        cons = [cons, cone( rho(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1),...
            1/norminv(1-epsilon)*p_EB(r,t) )];
        cons = [cons, cone( rho(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1),...
            1/norminv(1-epsilon)*(P_EB-p_EB(r,t)) )];
    end
end
cons = [cons, q_EB == eta_EB * p_EB, nu == eta_EB * rho];
% bilinear terms
binary_interger_set = ones(1, N);
for i = 1:N
    binary_interger_set(i) = 100 * 2^(i-1);
end
cons = [cons, n_el == z' * binary_interger_set'];
for r = 1:n_r
    for t = 1:T
        cons = [cons, 1000*p_P2HH(r,t) == binary_interger_set*e1(:,r,t)];
        cons = [cons, e1(:,r,t) <= M * z, e1(:,r,t) >= -M * z];
        cons = [cons, e1(:,r,t) <= M * (ones(N,1)-z) + p_p2hh(r,t)*ones(N,1),...
           e1(:,r,t) >= -M * (ones(N,1)-z) + p_p2hh(r,t)*ones(N,1)];
        cons = [cons,1000* Delta(r,t) == binary_interger_set*e2(:,r,t)];
        cons = [cons, e2(:,r,t) <= M * z, e2(:,r,t) >= -M * z];
        cons = [cons, e2(:,r,t) <= M * (ones(N,1)-z) + delta(r,t)*ones(N,1),...
           e2(:,r,t) >= -M * (ones(N,1)-z) + delta(r,t)*ones(N,1)];
        cons = [cons, 1000*h_P2HH(r,t) == binary_interger_set*e3(:,r,t)];
        cons = [cons, e3(:,r,t) <= M * z, e3(:,r,t) >= -M * z];
        cons = [cons, e3(:,r,t) <= M * (ones(N,1)-z) + h_p2hh(r,t)*ones(N,1),...
           e3(:,r,t) >= -M * (ones(N,1)-z) + h_p2hh(r,t)*ones(N,1)];
        cons = [cons, 1000*Pi(r,t) == binary_interger_set*e4(:,r,t)];
        cons = [cons, e4(:,r,t) <= M * z, e4(:,r,t) >= -M * z];
        cons = [cons, e4(:,r,t) <= M * (ones(N,1)-z) + pi(r,t)*ones(N,1),...
           e4(:,r,t) >= -M * (ones(N,1)-z) + pi(r,t)*ones(N,1)];
        cons = [cons, 1000*q_EXC(r,t) == binary_interger_set*e5(:,r,t)];
        cons = [cons, e5(:,r,t) <= M * z, e5(:,r,t) >= -M * z];
        cons = [cons, e5(:,r,t) <= M * (ones(N,1)-z) + q_exc(r,t)*ones(N,1),...
           e5(:,r,t) >= -M * (ones(N,1)-z) + q_exc(r,t)*ones(N,1)];
        cons = [cons, 1000*Lambda(r,t) == binary_interger_set*e6(:,r,t)];
        cons = [cons, e6(:,r,t) <= M * z, e6(:,r,t) >= -M * z];
        cons = [cons, e6(:,r,t) <= M * (ones(N,1)-z) + lambda(r,t)*ones(N,1),...
           e6(:,r,t) >= -M * (ones(N,1)-z) + lambda(r,t)*ones(N,1)];
    end
end
% exchange heat
for r = 1:n_r
    for t = 1:T
        cons = [cons, cone( lambda(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z, 1),...
            1/norminv(1-epsilon)*q_exc(r,t))];
    end
end
% p2hh
cons = [cons, p_p2hh == h_p2hh + q_p2hh];
cons = [cons, delta == pi + kappa];
cons = [cons, n_H2 == 3.6E6/(U_tn*F) * h_P2HH];
cons = [cons, gamma == 3.6E6/(U_tn*F) * Pi];
% y
for r = 1:n_r
    for t = 1:T
        cons = [cons, ones(1, n_P2HH_corner)*y(:,r,t) == 1];
        cons = [cons, ones(1, n_P2HH_corner)*ksi(:,r,t) == 0];
        for k = 1:n_P2HH_corner
                cons = [cons, cone(ksi(k, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                        1/norminv(1-epsilon)*y(k, r, t))];
                cons = [cons, cone(ksi(k, r, t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
                        1/norminv(1-epsilon)*(1-y(k, r, t)))];
        end
    end
end
%p2hh
for r = 1:n_r
    for t = 1:T
        cons = [cons, P_H2 * y(:,r,t) == h_p2hh(r,t)];
        cons = [cons, P_H2 * ksi(:,r,t) == pi(r,t)];
        cons = [cons, P_Heat * y(:,r,t) == q_p2hh(r,t)];
        cons = [cons, P_Heat * ksi(:,r,t) == kappa(r,t)];
        cons = [cons, Temp_P2HH * y(:,r,t) == Temp(r,t)];
        cons = [cons, Temp_P2HH * ksi(:,r,t) == mu(r,t)];
    end
end
cons = [cons, Temp(:,1) == Temp_init];
cons = [cons, mu(:,1) == 0];
%temperature evolution
for r = 1:n_r
    for t = 2:T
        cons = [cons, cone( mu(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
           1/norminv(1-epsilon)*(Temp(r,t)-T_min) )];
        cons = [cons, cone( mu(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1), ...
           1/norminv(1-epsilon)*(T_max-Temp(r,t)) )];   
    end
end
for  r = 1:n_r
    for t = 1:T-1
        cons = [cons, Temp(r,t+1) == Temp(r,t) + 1/C * (q_p2hh(r,t)...
            -q_exc(r,t)-1/R_eqv * (Temp(r,t) - Ta) )];
        cons = [cons, mu(r,t+1) == mu(r,t) + 1/C *(kappa(r,t)-lambda(r,t) - 1/R_eqv * mu(r,t))];
    end
end
%temperature at end of the day
for r = 1:n_r
    t = T;
    cons = [cons, cone( ( (1-1/(R_eqv*C))*mu(r,t)+1/C*kappa(r,t)-1/C*lambda(r,t) ) * sqrtm(Sigma(:,:,r,t))*ones(Z,1),...
        1/norminv(1-epsilon)*((1-1/(R_eqv*C))*Temp(r,t)+1/C*(q_p2hh(r,t)-q_exc(r,t)+1/R_eqv*Ta)-T_min  ) )];
    cons = [cons, cone( ( (1-1/(R_eqv*C))*mu(r,t)+1/C*kappa(r,t)-1/C*lambda(r,t) ) * sqrtm(Sigma(:,:,r,t))*ones(Z,1),...
        1/norminv(1-epsilon)*(-(1-1/(R_eqv*C))*Temp(r,t)-1/C*(q_p2hh(r,t)-q_exc(r,t)+1/R_eqv*Ta)+T_max  ) )];   
end
% % tank capacity evolution
cons = [cons, m_H2(:,1)==n_H2(:,1), gamma(:,1) == upsilon(:,1)];
for r = 1:n_r
    for t = 2:T
        cons = [cons, m_H2(r,t)== m_H2(r,t-1)+n_H2(r,t),...
            upsilon(r,t) == upsilon(r,t-1) + gamma(r,t)];
    end
end
% end of the day
for r = 1:n_r
    T = 24;
    cons = [cons, cone( sqrtm(Sigma(:,:,r,T))*(upsilon(r,T)*ones(Z,1)),...
        1/norminv(1-epsilon)* (m_tank-m_H2(r,T)))];
end
%converter
for r = 1:n_r
    for t = 1:T
        cons = [cons, cone( Delta(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1),...
            1/norminv(1-epsilon) * (P_converter * eta_converter - p_P2HH(r,t)) )];
    end
end
%compressor
for r = 1:n_r
    for t = 1:T
        cons = [cons, cone( gamma(r,t)*sqrtm(Sigma(:,:,r,t))*ones(Z,1),...
            1/norminv(1-epsilon) * (m_comp - n_H2(r,t)) )];
    end
end
% cons = [cons, P_EB == 0];
% cons = [cons, n_el == 0];
%% solving model and result presentation
options = sdpsettings( 'solver', 'gurobi','gurobi.TimeLimit', 10800);
result = optimize ( cons, obj, options);
cputime_used = cputime - cputime_start;
obj = value(obj);
n_el = value(n_el);
P_EB = value(P_EB);
z = value(z);
u_su = value(u_su); %start-up
u_sd = value(u_sd); %shut-down
u = value(u); %ON/OFF
x = value(x); %CHP corner point weights
alpha = value(alpha);
p_trans = value(p_trans);
beta = value(beta);%trans 
p_EB = value(p_EB);
rho = value(rho);%EB p 
q_EB = value(q_EB);
nu = value(nu);%EB q
n_H2 = value(n_H2); % kg
gamma = value(gamma);
p_P2HH = value(p_P2HH);%MW
Delta = value(Delta);
p_p2hh = value(p_p2hh);
delta = value(delta);
q_EXC = value(q_EXC);
Lambda = value(Lambda);
q_exc = value(q_exc);
lambda = value(lambda);
h_P2HH = value(h_P2HH);
Pi = value(Pi);
h_p2hh = value(h_p2hh);
pi = value(pi);
q_p2hh = value(q_p2hh);
kappa = value(kappa);
y = value(y);
ksi = value(ksi);
Temp = value(Temp);
mu = value(mu);
e1 = value(e1);
e2 = value(e2);
e3 = value(e3);
e4 = value(e4);
e5 = value(e5);
e6 = value(e6);
C_times_x = value(C_times_x);
m_H2 = value(m_H2);
upsilon = value(upsilon);
m_tank = value(m_tank);
P_converter = value(P_converter);
m_comp = value(m_comp);
toc

save('CC_P2HH_EB_99.mat',"obj", "n_el", "P_EB", "z","u_su","u_sd", ...
    "u", "x", "alpha","p_trans","beta","p_EB", "rho","q_EB", "nu", "n_H2",...
    "gamma","p_P2HH", "Delta", "p_p2hh","delta","q_EXC", "Lambda","q_exc",...
    "lambda","h_P2HH","Pi", "h_p2hh","pi","q_p2hh","kappa","y","ksi",...
    "Temp","mu","e1","e2","e3","e4","e5","e6","C_times_x","m_H2","upsilon", "m_tank",...
    "P_converter", "m_comp");

T=24;
n_r = 10;
P = [5, 3.97, 13.9, 15.5]; %MW
Q = [0, 5.14, 10.7, 0]; %MW
plot(Q,P);
xlim([0 12]);
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
saveas(gcf, 'CC_P2HH_EB_99', 'epsc');

disp("inverse flow");
weights' * sum(abs(p_trans),2)
disp("obj");
obj





% for r = 1:n_r
%     gamma_vec = [];
%     for t = 1:T
%         gamma_vec = [gamma_vec; gamma(r,t)*ones(Z,1)];
%     end
%     cons = [cons, cone( 2E4*sqrtm(Sigma_wholeDay(:,:,r)) * gamma_vec ,...
%         1/norminv(1-epsilon)* 2E4*(m_tank-n_H2(r,:)*ones(T,1))) ];
% end

% %power grid
% for r = 1:n_r
%     for t = 1:T
%         for l = 1:L
%             cons = [cons, cone( sqrtm(Sigma(:,:,r,t)) * ( (H_CHP(l)*P*alpha(:,r,t)...
%                 +H_trans(l)*beta(r,t)-H_P2HH(l)*Delta(r,t)) * ones(Z,1) + (H_wind(l,:))' ) ,...
%                 1/norminv(1-epsilon) * ( -H_CHP(l)*P*x(:,r,t)-H_trans(l)*p_trans(r,t)...
%                 -H_wind(l,:)*m(:,r,t)+H_P2HH(l)*p_P2HH(r,t)+H_D(l,:)*d_p(:,r,t)+fmax(l) ) )];
%             cons = [cons, cone( sqrtm(Sigma(:,:,r,t)) * ( (H_CHP(l)*P*alpha(:,r,t)...
%                 +H_trans(l)*beta(r,t)-H_P2HH(l)*Delta(r,t)) * ones(Z,1) + (H_wind(l,:))' ) ,...
%                 1/norminv(1-epsilon) * (H_CHP(l)*P*x(:,r,t)+H_trans(l)*p_trans(r,t)...
%                 +H_wind(l,:)*m(:,r,t)-H_P2HH(l)*p_P2HH(r,t)-H_D(l,:)*d_p(:,r,t)+fmax(l)) )];
%         end
%     end

% end


        





