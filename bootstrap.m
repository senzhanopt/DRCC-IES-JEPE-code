clear; clc;
load wind_profile;
z = 2;
n_sce = 1000;
n_r = 10;
T = 24;
% bootstrap for out of sample analysis
Windfarm1_sim = zeros(T, n_r * n_sce);
Windfarm2_sim = zeros(T, n_r * n_sce);
for r = 1:n_r
    for sce = 1:n_sce
        random_number = randi(n_sce);
        Windfarm1_sim(:,n_r*(sce-1)+r) = Windfarm1(:, n_r*(random_number-1)+r);
        Windfarm2_sim(:,n_r*(sce-1)+r) = Windfarm2(:, n_r*(random_number-1)+r);
    end
end
windOutput = zeros(n_sce, z, n_r, T);
for r = 1:n_r
    for t = 1:T
        for n = 1:n_sce
            windOutput(n,1,r,t) = Windfarm1(t, 10*(n-1)+r); 
            windOutput(n,2,r,t) = Windfarm2(t, 10*(n-1)+r);
        end
    end
end
windOutput_sim = zeros(n_sce, z, n_r, T);
for r = 1:n_r
    for t = 1:T
        for n = 1:n_sce
            windOutput_sim(n,1,r,t) = Windfarm1_sim(t, 10*(n-1)+r); 
            windOutput_sim(n,2,r,t) = Windfarm2_sim(t, 10*(n-1)+r);
        end
    end
end
windOutput_mean = zeros(1,z,n_r,T);
for r = 1:n_r
    for t = 1:T
        windOutput_mean(:,:,r,t) = mean(windOutput(:,:,r,t));
    end
end
windOutput_error_sim = zeros(n_sce, z, n_r, T);
for r = 1:n_r
    for t = 1:T
        for n = 1:n_sce
            windOutput_error_sim(n,:,r,t) =...
                windOutput_sim(n,:,r,t)-windOutput_mean(:,:,r,t);
        end
    end
end
save('wind_sim.mat',"windOutput_error_sim");

