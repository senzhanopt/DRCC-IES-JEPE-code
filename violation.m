clear;clc;
load DRCC_P2HH_EB_60
load wind
T = 24;
C = (2.5*0.66*4200+21000*0.66*0.3/56*65.87)/(3.6E6); %specific heat, kW/Celsius degree
R_eqv = 0.167E3; %C / kW
SU = 14; 
RU = 14; %MW
RD = 14;
SD = 14;  %MW
P = [5, 3.97, 13.9, 15.5]; %MW
Q = [0, 5.14, 10.7, 0]; %MW
n_r = 10;
windCap = 18; %,MW for each wind farm
windOutput_error_mw = windCap.* windOutput_error;
for r = 1:n_r
    for t  = 1:T
        for k = 1:4
            viol_x_big(k,r,t)=sum((x(k,r,t)+alpha(k,r,t)* sum(windOutput_error_mw(:,:,r,t), 2))>=-1E-10)/1000;
            viol_x_less(k,r,t)=sum((x(k,r,t)+alpha(k,r,t)* sum(windOutput_error_mw(:,:,r,t), 2))<=1+1E-10)/1000;
        end
    end
end
viol_x = mean(viol_x_big+viol_x_less, "all")/2;
for r = 1:n_r
    t = 1;
    viol_start(r) = sum( P * x(:,r,t) + P * alpha(:,r,t) * sum(windOutput_error_mw(:,:,r,t), 2) <= SU+1E-10)/1000;
end  
for r = 1:n_r
    for t = 2:T
        viol_rampup(r,t) = sum(P * x(:,r,t) + P * alpha(:,r,t) * sum(windOutput_error_mw(:,:,r,t), 2) - ...
            P * x(:,r,t-1) + P * alpha(:,r,t-1) * sum(windOutput_error_mw(:,:,r,t-1), 2) <= SU*(1-u(r,t-1))...
            + RU * u(r,t-1) )/1000;
        viol_rampdown(r,t) = sum(P * x(:,r,t) + P * alpha(:,r,t) * sum(windOutput_error_mw(:,:,r,t), 2) - ...
            P * x(:,r,t-1) + P * alpha(:,r,t-1) * sum(windOutput_error_mw(:,:,r,t-1), 2) >= -SD*(1-u(r,t))...
            - RD * u(r,t) )/1000;
    end
end
for r = 1:n_r
    for t  = 1:T
        viol_EB_big(r,t) = sum(p_EB(r,t) + rho(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) >= -1E-10)/1000;
        viol_EB_less(r,t) = sum(p_EB(r,t) + rho(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) <= 1E-10+P_EB)/1000;
        viol_exc(r,t) =  sum(q_exc(r,t) + lambda(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) >= -1E-10)/1000;
    end
end
for r = 1:n_r
    for t  = 1:T
        for k = 1:4
            viol_y_big(k,r,t)=sum((y(k,r,t)+ksi(k,r,t)* sum(windOutput_error_mw(:,:,r,t), 2))>=-1E-10)/1000;
            viol_y_less(k,r,t)=sum((y(k,r,t)+ksi(k,r,t)* sum(windOutput_error_mw(:,:,r,t), 2))<=1+1E-10)/1000;
        end
    end
end
viol_y = mean(viol_y_big+viol_y_less, "all")/2;
for r = 1:n_r
    for t = 2:T
        violTemp_big(r,t) = sum(Temp(r,t) + mu(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) >= 60-1E-10)/1000;
        violTemp_less(r,t) = sum(Temp(r,t) + mu(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) <= 80+1E-10)/1000;
    end
end
for r = 1:n_r
    t = T;
    viol_endTemp_big(r) = sum(Temp(r,t) + mu(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) + 1/C * ...
        (q_p2hh(r,t)+kappa(r,t)*sum(windOutput_error_mw(:,:,r,t), 2)-q_exc(r,t)-lambda(r,t)*sum(windOutput_error_mw(:,:,r,t), 2)...
        -1/R_eqv * (Temp(r,t) + mu(r,t)*sum(windOutput_error_mw(:,:,r,t), 2)-25) )...
        >= 60-1E-10)/1000;
    viol_endTemp_less(r) = sum(Temp(r,t) + mu(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) + 1/C * ...
        (q_p2hh(r,t)+kappa(r,t)*sum(windOutput_error_mw(:,:,r,t), 2)-q_exc(r,t)-lambda(r,t)*sum(windOutput_error_mw(:,:,r,t), 2)...
        -1/R_eqv * (Temp(r,t) + mu(r,t)*sum(windOutput_error_mw(:,:,r,t), 2)-25))...
        <= 80+1E-10)/1000;     
end
for r = 1:n_r
    t = T;
    viol_tank(r) =  sum(m_H2(r,t) + upsilon(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) <= m_tank+1E-10)/1000;
end
for r = 1:n_r
    for t  = 1:T
        viol_conv(r) =  sum(p_P2HH(r,t) + Delta(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) <= P_converter*0.95+1E-10)/1000;
    end
end
for r = 1:n_r
    for t  = 1:T
        viol_comp(r) =  sum(n_H2(r,t) + gamma(r,t)*sum(windOutput_error_mw(:,:,r,t), 2) <= m_comp+1E-10)/1000;
    end
end