%% PPE HS 20 - Kriging
% Yue Pan & Fandré Josianne

function [semivar_func] = fit_semivariogram(data,function_type,compare_or_not)
%fit_semivariogram fit the empirical semivariogram
% Input:   
% - data: (3 x n_data)
% - function_type: string, selected from 'spherical', 'exponential' and 'squared_exponential'  
% - compare_or_not: plot each model's variogram and do the comparison or % not (0 or 1)
% Output:
% - semivar_function: function handle

% reference: https://vsp.pnnl.gov/help/Vsample/Kriging_Variogram_Model.htm

% The figure shows an experimental variogram with a variogram model fitted to it. Each red square is a lag of the experimental variogram. The x-axis represents the distance between pairs of points, and the y-axis represents the calculated value of the variogram, where a greater value indicates less correlation between pairs of points. This particular variogram shows a spatial relationship well suited for geostatistical analysis since pairs of points are more correlated the closer they are together and become less correlated the greater the distance between points.
% This graphic also illustrates three important parameters that control the fit of the variogram model. The nugget is the y-intercept of the variogram. In practical terms, the nugget represents the small-scale variability of the data. A portion of that short range variability can be the result of measurement error. The range is the distance after which the variogram levels off. The physical meaning of the range is that pairs of points that are this distance or greater apart are not spatially correlated. The sill is the total variance contribution, or the maximum variability between pairs of points.

%% Default setting
compare_or_not=0;
if nargin<3
   compare_or_not=1;
end

%% Prepare
% calculate the distance and value difference between each pixel pair
diff_mat = pdist2(data(3,:)',data(3,:)');     % delta z
dist_mat = pdist2(data(1:2,:)',data(1:2,:)'); % dist(x,y)
% Mat2Vec
diff_vec = reshape(diff_mat,[],1);
dist_vec = reshape(dist_mat,[],1);

total_sample_num = size(data,2); %n_data
disp(['Total variogram sample number is:', num2str(total_sample_num)]);

%% sampled semivariance

% distance bin seperation
bin_number = max(10, round(total_sample_num/50)); % number of variogram lags
gamma = zeros(1,bin_number);
disp(['Variogram lag number is:', num2str(bin_number)]);

% set range
range = 0.8 * max(dist_vec);
dist_step_bin = range/bin_number;

% calculate semivariances
% semivariance:
% r_ij = 0.5 * E[(z_i - z_j)^2]= sigma^2 - C_ij
for i=1:bin_number
    cur_bin_dist_min = (i-1) * dist_step_bin;
    cur_bin_dist_max = i * dist_step_bin;
    diff_vec_cur_bin = diff_vec(dist_vec > cur_bin_dist_min & dist_vec<cur_bin_dist_max);
    diff_squared_sum = sum(diff_vec_cur_bin.^2);
    cur_bin_sample_count = length(diff_vec_cur_bin);
    
    gamma(i) = 0.5* diff_squared_sum/(cur_bin_sample_count+1e-10); % E(dZ^2) = sum(dZ^2)/count % 1e-10 for make sure the denominator would not be zero
end

bin_cent_vec = dist_step_bin/2 : dist_step_bin : dist_step_bin*bin_number; % vec(begin, step, end)
sampled_semivariogram = [bin_cent_vec; gamma]; % variogram matrix fed to the Fit_model_to_data function


%% Fit model to data (using the provided function)
if compare_or_not
    semivar_func_sphe = Fit_model_to_data(sampled_semivariogram(:,:),'spherical');
    semivar_func_exp = Fit_model_to_data(sampled_semivariogram(:,:),'exponential');
    semivar_func_squexp = Fit_model_to_data(sampled_semivariogram(:,:),'squared_exponential');
else
    switch function_type
        case 'spherical'
            semivar_func_sphe = Fit_model_to_data(sampled_semivariogram(:,:),'spherical');
        case 'exponential'
            semivar_func_exp = Fit_model_to_data(sampled_semivariogram(:,:),'exponential');
        case 'squared_exponential'
            semivar_func_squexp = Fit_model_to_data(sampled_semivariogram(:,:),'squared_exponential');
    end
end

%% Plot the fitting result of different models
if compare_or_not
    figure(2);
    clf;
    set(gca, 'Fontname', 'Times New Roman','FontSize',14);
    scatter(bin_cent_vec, gamma, 20, [1,0,0], 'filled');
    hold on
    grid on
    plot(bin_cent_vec, semivar_func_sphe(bin_cent_vec),'LineWidth',2);
    plot(bin_cent_vec, semivar_func_exp(bin_cent_vec),'LineWidth',2);
    plot(bin_cent_vec, semivar_func_squexp(bin_cent_vec),'LineWidth',2);
    legend('empirical','spherical','exponential','squared exponential','Location','southeast','Fontname','Times New Roman','FontSize',12);
    ylabel('Gamma','Fontname', 'Times New Roman','FontSize',14);
    xlabel('Distance','Fontname', 'Times New Roman','FontSize',14);
    title('Semivariogram and fitted model','Fontname', 'Times New Roman','FontSize',16);
end

%% Select and output the best model
switch function_type
    case 'spherical'
        semivar_func = semivar_func_sphe;
    case 'exponential'
        semivar_func = semivar_func_exp; 
    case 'squared_exponential'
        semivar_func = semivar_func_squexp;
end

end

