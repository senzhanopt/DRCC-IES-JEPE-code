clear; clc; 
%% K means
D = 365;
T = 24;
k = 10; % no. of centroids
eload = xlsread("eload.xlsx","E1:E8760"); %DK Electricity demand 2019
DHload = xlsread("reference_data_set-v1.0.0.xlsx", "2.2 Demand Neighbors",...
    "E8:E8767"); %DK DH demand 
eload_normalized = eload ./ max(eload);
eload_normalized_reshaped = transpose(reshape(eload_normalized, [T, D]));
DHload_normalized = DHload ./ max(DHload);
DHload_normalized_reshaped = transpose(reshape(DHload_normalized, [T, D]));
% subplot(2,1,1);
% plot(1:D*T, eload_normalized);
% title("Normalized electric load");
% xlabel("T(hour)")
% subplot(2,1,2);
% plot(1:D*T, DHload_normalized);
% title("Normalized district heating load");
% xlabel("T(hour)")
% saveas(gcf, 'hourlyLoads', 'epsc');

e_dh_loads = zeros(D, 2*T);
for d = 1:D
    e_dh_loads(d,:) = [eload_normalized_reshaped(d,:),...
        DHload_normalized_reshaped(d,:)];
end
sumd_best = 1000;
idx_best = [];
Centroids_best = [];
for i = 1:1000
    [idx, Centroids, sumd] = kmeans(e_dh_loads, k);
    if sum(sumd) <= sumd_best
        sumd_best = sum(sumd);
        idx_best = idx;
        Centroids_best = Centroids;
    end
end
clustered_e_dh_demands = Centroids_best;
weights = zeros(k, 1);
for i = 1:k
    weights(i) = sum(idx_best == i);
end
save('cluster.mat', "clustered_e_dh_demands", 'weights');

%% PCA representation
mu = mean(e_dh_loads);
e_dh_loads_zeromean = e_dh_loads - mu;
Centroids_best_zeromean = Centroids_best - mu;
coeff = pca(e_dh_loads_zeromean);
components = 2;
coeff = transpose(coeff(:, 1:components));
z = zeros(D, components);
for d = 1:D
    z(d,:) = transpose(coeff * transpose(e_dh_loads_zeromean(d,:)));
end
z_c = zeros(k, components);
for i = 1:k
    z_c(i,:) = transpose(coeff * transpose(Centroids_best_zeromean(i,:)));
end
col = ['r','g','b','y','m','c','w','k','g','r'];
for i = 1:k
    z1 = z(idx == i, :);
    scatter(z1(:,1), z1(:,2),'MarkerFaceColor',col(i));
%     plot(z_c(i,1), z_c(i,2),'+', 'MarkerFaceColor',col(i),'LineWidth', 10);
    hold on;
end
hold off;
saveas(gcf, 'kmeans','epsc');