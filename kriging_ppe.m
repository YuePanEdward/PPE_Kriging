%% PPE HS 20 - Kriging
% Yue Pan & Fandré Josianne

function [Field_values, Field_variances, scale] = kriging_ppe(Data, Field_values_init,function_type)
%kriging_ppe: Do interpolation using Kriging (project parameter estimation course)
% Input:   
% -Data: (3 x n_data) matrix. This matrix contains in the first two rows the coordinates of 
% the data points and in the third row it contains their elevation values. n_data is the 
% number of data points  
% -Field_values_init: (n x m) matrix. This Matrix contains nothing interesting in the 
% beginning and then will later be filled with the interpolated values. It is only of use 
% to define the dimension of the output.
% -function_type: string. selected from 'exponential', 'squared exponential' and
% 'spherical'. It defines how the semivariogram is fitted
%                       
% Output:
% -Field_values: (n x m) matrix. Each entry is a BLUP interpolation result or a data value.
% -Field_variances: (n x m) matrix. Each entry is the estimated error variances.
% -scale: It represents the appropriate conversion factors  between Pixel/Matrix entries and 
% real world length scales

%% find missing values, fill the existing values into the grid with normalization 
[Field_values, Grid_parameters, coord_to_index, Data_to_estimate] = From_list_to_grid(Data, Field_values_init);

% output the scale (output 3)
scale = Grid_parameters(:,3);

%% Fit semivariogram function with regards to distance
semivar_function = fit_semivariogram(Data,function_type);


%% Calculate variogram 
% calculate standard deviation from the whole dataset, taking nan as missing value
sigma_0 = nanstd(Data(3,:));

known_count = size(Data,2);               % n
unknown_count = size(Data_to_estimate,2); % m

% Get distance pair matrix for known positions 
distMatC = pdist2(Data(1:2,:)',Data(1:2,:)'); % [n x n]

% Get distance pair matrix between known and missing positions
distMatc = pdist2(Data_to_estimate',Data(1:2,:)'); % [m x n]

% Covariance and semi-variogram matrix of existing points
Gamma_mat = semivar_function(distMatC); % semivar matrix [n x n]
C_mat = sigma_0^2 - Gamma_mat;          % covariance matrix [n x n]

% Covariance and semi-variogram vectors (matrix) for new points
gamma_vecs = semivar_function(distMatc)'; % semivar column vectors [n x m]
c_vecs = sigma_0^2 - gamma_vecs;          % covariance column vectors [n x m] 

one_vec =  ones(known_count, 1);   % [n x 1]
Z_vec = Data(3,:)';                % [n x 1]
inv_Gamma_mat = pinv(Gamma_mat);   % [n x n], pinv use the Moore-Penrose pseudo inverse, which can ease the situation of ill-condition matrix

%% Interpolation
% Refer to the paper 'The origins of Kriging'
% BLUP
Z_s_vec = gamma_vecs'*inv_Gamma_mat*Z_vec + (1 - gamma_vecs'*inv_Gamma_mat*one_vec) * pinv(one_vec'*inv_Gamma_mat*one_vec) * (one_vec'*inv_Gamma_mat*Z_vec); % BLUP [m x 1]

% MSE of the prediction
Z_s_var_mat = gamma_vecs'*inv_Gamma_mat*gamma_vecs - (1-gamma_vecs'*inv_Gamma_mat*one_vec).^2 * pinv(one_vec'*inv_Gamma_mat*one_vec); % Variance [m x m]
Z_s_var_vec = diag(Z_s_var_mat)'; % diagonal elements [1 x m]

% add the x,y coordinates
XYZ_s_vec = [Data_to_estimate; Z_s_vec']; % [3 x m]
XYZ_s_var_vec = [Data_to_estimate; Z_s_var_vec]; % [3 x m]

Field_values= Write_list_in_grid(Field_values, coord_to_index, XYZ_s_vec,3);   % output 1
Field_variances= Write_list_in_grid(Field_values_init, coord_to_index, XYZ_s_var_vec,3); % output 2

end

