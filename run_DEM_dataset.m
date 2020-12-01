%% PPE HS 20 - Kriging
% Yue Pan & Fandré Josianne
% Kriging interpolation for filling the missing value in a DEM datatset

clear; close all; clc;

%% Load data
data_path = ['.' filesep 'test_data' filesep];

load([data_path 'dtm_basedata_coarse_lossy.mat']);
load([data_path 'Field_values_init.mat']);

%% Preprocessing
% DTM to list
dtm = From_grid_to_list(dtm_basedata_coarse_lossy); % convert to [(m*n) x 3] matrix, each row is [x,y,z]
dtm_data = [dtm(2,:);dtm(1,:);dtm(3,:)]; % re-arrange x,y value

min_z = min(dtm(3,:)); max_z = max(dtm(3,:)); std_z = std(dtm(3,:));
min_z_thre = min_z-0.5*std_z; max_z_thre=max_z+0.5*std_z; z_thre= [min_z_thre, max_z_thre];
size_x = size(dtm_basedata_coarse_lossy,2); size_y = size(dtm_basedata_coarse_lossy,1); 

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

%% Plots
figure(1);
clf;

% Input image with missing parts
subplot(3,4,1)
imagesc(dtm_basedata_coarse_lossy, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Input Image','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation result
subplot(3,4,2)
imagesc(Field_values_sphe, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging estimation (spherical)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,3)
imagesc(Field_values_exp, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging estimation (exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,4)
imagesc(Field_values_expsq, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging estimation (squ. exponential)','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation variance
subplot(3,4,6)
imagesc(Field_variances_sphe)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging variance (spherical)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,7)
imagesc(Field_variances_exp)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging variance (exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,8)
imagesc(Field_variances_expsq)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging variance (squ. exponential)','Fontname','Times New Roman','FontSize',14);

% Baseline intepolation methods

% Natural neighbor
subplot(3,4,9)
imagesc(Natural_Neighbor, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Natural Neighbor','Fontname','Times New Roman','FontSize',14);

% Nearest neighbor
subplot(3,4,10)
imagesc(Nearest_Neighbor, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Nearest Neighbor','Fontname','Times New Roman','FontSize',14);

% Linear interpolation
subplot(3,4,11)
imagesc(Linear_Interp, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Linear Interpolation','Fontname','Times New Roman','FontSize',14);

% Cubic spline
subplot(3,4,12)
imagesc(Cubic, z_thre)
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Cubic Spline','Fontname','Times New Roman','FontSize',14);


