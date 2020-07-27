clear; clc;
load wind_profile;
z = 2;
n_sce = 1000;
n_r = 10;
T = 24;
windOutput = zeros(n_sce, z, n_r, T);
for r = 1:n_r
    for t = 1:T
        for n = 1:n_sce
            windOutput(n,1,r,t) = Windfarm1(t, 10*(n-1)+r); 
            windOutput(n,2,r,t) = Windfarm2(t, 10*(n-1)+r);
        end
    end
end
windOutput_mean = zeros(1,z,n_r,T);
for r = 1:n_r
    for t = 1:T
        windOutput_mean(:,:,r,t) = mean(windOutput(:,:,r,t));
    end
end
windOutput_error = zeros(n_sce, z, n_r, T);
for r = 1:n_r
    for t = 1:T
        for n = 1:n_sce
            windOutput_error(n,:,r,t) =...
                windOutput(n,:,r,t)-windOutput_mean(:,:,r,t);
        end
    end
end
Sigma = zeros(z,z,n_r,T);
for r = 1:n_r
    for t = 1:T
        Sigma(:,:,r,t) = (transpose(windOutput_error(:,:,r,t))...
            *windOutput_error(:,:,r,t))./n_sce;
    end
end
Sigma_coupling = zeros(2*z,2*z,n_r,T);
for r = 1:n_r
    for t = 2:T
        Sigma_coupling(:,:,r,t) = (transpose([windOutput_error(:,:,r,t-1),windOutput_error(:,:,r,t)])...
            *[windOutput_error(:,:,r,t-1),windOutput_error(:,:,r,t)])./n_sce;
    end
end
Sigma_wholeDay = zeros(T*z, T*z, n_r);
for r = 1:n_r
    wind_error_day = []; % 1000 * 48
    for t = 1:24
        wind_error_day = [wind_error_day, windOutput_error(:,:,r,t)];
    end
    Sigma_wholeDay(:,:,r) = transpose(wind_error_day) * wind_error_day/n_sce;
end
save('wind.mat',"windOutput_mean", "windOutput_error", "Sigma", "Sigma_coupling");