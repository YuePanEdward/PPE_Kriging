%% PPE HS 20 - Kriging
% Yue Pan & Fandré Josianne
% Kriging interpolation for filling the missing value in a DEM datatset

clear; close all; clc;

%% Load data
load('dtm_basedata_coarse_lossy.mat')
load('Field_values_init.mat')

%% Preprocessing
% DTM to list
dtm = From_grid_to_list(dtm_basedata_coarse_lossy); % convert to [(m*n) x 3] matrix, each row is [x,y,z]
dtm_data = [dtm(2,:);dtm(1,:);dtm(3,:)]; % re-arrange x,y value

min_z = min(dtm(3,:)); max_z = max(dtm(3,:)); std_z = std(dtm(3,:));
min_z_thre = min_z-std_z; max_z_thre=max_z+std_z; z_thre= [min_z_thre, max_z_thre];

%sampled_data = zeros(size(dtm_basedata_coarse_lossy));  % empty matrix for storing the interpolation result

%% Kriging (with spherical, exponential and squ. exponential semivariogram model)
% Calculate missing values via ordinary Kriging
%used_model = 'spherical'; % select from 'exponential', 'squared exponential' and  'spherical'
% spherical model
[Field_values_sphe, Field_variances_sphe, scale] = kriging_ppe(dtm_data, Field_values_init, 'spherical');
% exponential model
[Field_values_exp, Field_variances_exp, scale] = kriging_ppe(dtm_data, Field_values_init, 'exponential');
% squared exponential model
[Field_values_expsq, Field_variances_expsq, scale] = kriging_ppe(dtm_data, Field_values_init, 'squared_exponential');

%% Baseline methods (linear, cubic, nature, nearest interpolation)
[Nearest_Neighbor, Linear_Interp, Natural_Neighbor, Cubic] = Make_comparative_interpolation(dtm_data, Field_values_init);

%% Assign interpolation results
%[missing_rows, missing_cols]=find(isnan(dtm_basedata_coarse_lossy));    % missing pixel coordinates for interpolating
%missing_ind = sub2ind(size(dtm_basedata_coarse_lossy),missing_rows,missing_cols); % index of the missing position
%sampled_data(missing_ind) = Field_values(missing_ind); %image with only the missing part

%% Plots
figure(1)
% Input image with missing parts
subplot(3,4,1)
imagesc(dtm_basedata_coarse_lossy, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Input Image','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation result
subplot(3,4,2)
imagesc(Field_values_sphe, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Kriging estimation (spherical)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,3)
imagesc(Field_values_exp, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Kriging estimation (exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,4)
imagesc(Field_values_expsq, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Kriging estimation (squ. exponential)','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation variance
subplot(3,4,6)
imagesc(Field_variances_sphe)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Kriging variance (spherical)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,7)
imagesc(Field_variances_exp)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Kriging variance (exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,8)
imagesc(Field_variances_expsq)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Kriging variance (squ. exponential)','Fontname','Times New Roman','FontSize',14);

% Baseline intepolation methods

% Natural neighbor
subplot(3,4,9)
imagesc(Natural_Neighbor, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Natural Neighbor','Fontname','Times New Roman','FontSize',14);

% Nearest neighbor
subplot(3,4,10)
imagesc(Nearest_Neighbor, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Nearest Neighbor','Fontname','Times New Roman','FontSize',14);

% Linear interpolation
subplot(3,4,11)
imagesc(Linear_Interp, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Linear Interpolation','Fontname','Times New Roman','FontSize',14);

% Cubic spline
subplot(3,4,12)
imagesc(Cubic, z_thre)
colorbar
axis equal
xlim([1,40])
ylim([1,40])
title('Cubic Spline','Fontname','Times New Roman','FontSize',14);


