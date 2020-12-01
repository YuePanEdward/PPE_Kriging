%% PPE HS 20 - Kriging
% Yue Pan & Fandré Josianne
% Kriging interpolation for the sparse dataset.

clear; close all; clc;

%% Load data
data_path = ['.' filesep 'test_data' filesep];

load([data_path 'Data_sparse_HW3.mat']);
load([data_path 'Field_values_init.mat']);
load([data_path 'Original_image.mat']); % ground truth

%% Preprocessing
% list to grid field
field_data = Data_sparse_HW3;
field = From_list_to_grid(field_data,Field_values_init);

min_z = min(field_data(3,:)); max_z = max(field_data(3,:)); std_z = std(field_data(3,:));
min_z_thre = min_z-0.5*std_z; max_z_thre=max_z+0.5*std_z; z_thre= [min_z_thre, max_z_thre];
diff_z_thre = [0, 0.5*(max_z - min_z)];
size_x = size(field,2); size_y = size(field,1); 
interpolate_count=size_x*size_y-size(field_data,2);

%% Kriging (with spherical, exponential and squ. exponential semivariogram model)
% Calculate missing values via ordinary Kriging
%used_model = 'spherical'; % select from 'exponential', 'squared exponential' and  'spherical'
% spherical model
[Field_values_sphe, Field_variances_sphe, scale] = kriging_ppe(field_data, Field_values_init, 'spherical');
% exponential model
[Field_values_exp, Field_variances_exp, scale] = kriging_ppe(field_data, Field_values_init, 'exponential');
% squared exponential model
[Field_values_expsq, Field_variances_expsq, scale] = kriging_ppe(field_data, Field_values_init, 'squared_exponential');

%% Baseline methods (linear, cubic, nature, nearest interpolation)
[Nearest_Neighbor, Linear_Interp, Natural_Neighbor, Cubic] = Make_comparative_interpolation(field_data, Field_values_init);

%% Plots
figure(1);
clf;

% Input image with missing parts
subplot(3,4,1)
imagesc(field, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Input Image','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation result
subplot(3,4,2)
imagesc(Field_values_sphe, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging estimation (spherical)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,3)
imagesc(Field_values_exp, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging estimation (exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,4)
imagesc(Field_values_expsq, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging estimation (squ. exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,5)
imagesc(Original_image, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Underlying function (ground truth)','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation variance
subplot(3,4,6)
imagesc(Field_variances_sphe);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging variance (spherical)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,7)
imagesc(Field_variances_exp);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging variance (exponential)','Fontname','Times New Roman','FontSize',14);

subplot(3,4,8)
imagesc(Field_variances_expsq);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Kriging variance (squ. exponential)','Fontname','Times New Roman','FontSize',14);

% Baseline intepolation methods

% Natural neighbor
subplot(3,4,9)
imagesc(Natural_Neighbor, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Natural Neighbor','Fontname','Times New Roman','FontSize',14);

% Nearest neighbor
subplot(3,4,10)
imagesc(Nearest_Neighbor, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Nearest Neighbor','Fontname','Times New Roman','FontSize',14);

% Linear interpolation
subplot(3,4,11)
imagesc(Linear_Interp, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Linear Interpolation','Fontname','Times New Roman','FontSize',14);

% Cubic spline
subplot(3,4,12)
imagesc(Cubic, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Cubic Spline','Fontname','Times New Roman','FontSize',14);


%% Plot the deviation from the ground truth

figure(3);
clf;

% Input image with missing parts
subplot(2,4,1)
imagesc(Original_image, z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title('Underlying function (ground truth)','Fontname','Times New Roman','FontSize',14);

% Kriging interpolation result
rmse_sphe = sqrt(nanmean((Field_values_sphe-Original_image).^2, 'all'));
subplot(2,4,2)
imagesc(abs(Field_values_sphe-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Kriging (spherical) RMSE:', num2str(rmse_sphe)],'Fontname','Times New Roman','FontSize',14);

rmse_exp = sqrt(nanmean((Field_values_exp-Original_image).^2, 'all'));
subplot(2,4,3)
imagesc(abs(Field_values_exp-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Kriging (exp.) RMSE:', num2str(rmse_exp)], 'Fontname','Times New Roman','FontSize',14);

rmse_expsq = sqrt(nanmean((Field_values_expsq-Original_image).^2, 'all'));
subplot(2,4,4)
imagesc(abs(Field_values_expsq-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Kriging (squ. exp.) RMSE:', num2str(rmse_expsq)],'Fontname','Times New Roman','FontSize',14);

% Baseline intepolation methods

% Natural neighbor
rmse_naturenei = sqrt(nanmean((Natural_Neighbor-Original_image).^2, 'all'));
subplot(2,4,5)
imagesc(abs(Natural_Neighbor-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Natural Neighbor RMSE:', num2str(rmse_naturenei)], 'Fontname','Times New Roman','FontSize',14);

% Nearest neighbor
rmse_nearnei = sqrt(nanmean((Nearest_Neighbor-Original_image).^2, 'all'));
subplot(2,4,6)
imagesc(abs(Nearest_Neighbor-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Nearest Neighbor RMSE:', num2str(rmse_nearnei)], 'Fontname','Times New Roman','FontSize',14);

% Linear interpolation
rmse_li = sqrt(nanmean((Linear_Interp-Original_image).^2, 'all'));
subplot(2,4,7)
imagesc(abs(Linear_Interp-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Linear Interpolation RMSE:', num2str(rmse_li)], 'Fontname','Times New Roman','FontSize',14);

% Cubic spline
rmse_cs = sqrt(nanmean((Cubic-Original_image).^2, 'all'));
subplot(2,4,8)
imagesc(abs(Cubic-Original_image),diff_z_thre);
colorbar
axis equal
xlim([1,size_x])
ylim([1,size_y])
title(['Cubic Spline RMSE:', num2str(rmse_cs)], 'Fontname','Times New Roman','FontSize',14);


