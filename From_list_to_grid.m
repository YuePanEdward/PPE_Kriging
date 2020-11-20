function [ Field_values, Grid_parameters, coord_to_index, Data_to_estimate]=From_list_to_grid(Data,Field_values_init)

% [ Field_values, Grid_parameters, Data_to_estimate]=From_list_to_grid(Data,Field_values_init)
%
% This function takes as Inputs two Matrices containing the data to be
% placed in a Grid (Data, (3,n_data)) and a Matrix "Field_values_init"
% (n,m) filled with NaN's serving as an initialization for the embedding
% of the list into a grid.
%
% Output 1 consists of a grid "Field_values". It is scaled, s.t. all the
% locations associated with the data lie inside and its entries are NaN,
% if there is no data for this location and Z(X,Y) otherwise. Scaling and
% finding the corresponding values is done with the help of a mapping
% function, which clears up the relation between grid points and
% coordinates.
% Output 2 "Grid_parameters" gives some hints on how this mapping function
% looks like by archiving which coordinate difference corresponds to which
% step size in the matrix.
% Output 3 "coord_to_index" is a function handle representing the function
% which maps coordinates onto indices.
% Output 4 "Data_to_estimate" is a matrix containing all the coordinates of
% those grid points, for which there is no data available. The coordinates
% of those grid_points will be use for interpolation later on.
%
% The formats are:
%       Data_to_estimate     = [X_s1 .......        X_sn_est]    (2,n_est)
%                              [Y_s1 .......        Y_sn_est]
%       Data                 = [X_p1 .......        X_pn_data]   (3,n_data)
%                              [Y_p1 .......        Y_pn_data]
%                              [Z_p1 .......        Z_pn_data]
%       Field_values_init    = [NaN ....            ... NaN ]         (n,m)
%                              [NaN ....            ... NaN ]
%       Field_values         = [NaN ....    Z(X_13) ... NaN ]         (n,m)
%                              [NaN ....       ... Z(X_n,m) ]
%       Grid_parameters      = [min_X max_x delta_X]                  (2,3)
%                              [min_Y max_y delta_Y]
%       coord_to_index       = function handle

% Get size of the Matrix, in which the list should be embedded
[n,m]=size(Field_values_init);

% Extract maximum and minimum coordinates from data and plug them into
% matrix minmax_mat = [x_min  x_max ]
%                     [y_min  y_max ]
%                     [z_min  z_max ]

minmax_mat=zeros(3,2);
minmax_mat(:,1)=min(Data,[],2);
minmax_mat(:,2)=max(Data,[],2);

% Make the mapping from coordinates to grid points coordinate_to_index and
% its inverse index_to_coordinate

% 1. Find the step sizes and archive everything in "Grid_parameters"
delta_y=(minmax_mat(2,2)-minmax_mat(2,1))/(n-1);
delta_x=(minmax_mat(1,2)-minmax_mat(1,1))/(m-1);
Grid_parameters=zeros(2,3);
Grid_parameters(:)=[minmax_mat(1,1),minmax_mat(2,1), minmax_mat(1,2),...
    minmax_mat(2,2), delta_x, delta_y];


% 2. Make index_to_coord. This is a function handle working with the syntax
% [x,y]=index_to_coord(i,j)
index_to_coord=@(i,j) deal(minmax_mat(1,1)+delta_x*(j-1),minmax_mat(2,1)+delta_y*(i-1));

% 3. Make coord_to_index. This is a function handle working with the syntax
% [i,j]=coord_to_index(x,y)
coord_to_index=@(x,y) deal(round((y-minmax_mat(2,1))/(delta_y)+1),...
    round((x-minmax_mat(1,1))/(delta_x)+1));

% Extract coordinates from data and use mapping to find corresponding
% indices. Write the values from the data points into the grid at the
% corresponding locations
Field_values=Write_list_in_grid(Field_values_init,coord_to_index,Data,3);

% Extract the indices of those Grid points, where there is no data.
Nan_index=find(isnan(Field_values));
[row_ind, col_ind]=ind2sub(size(Field_values),Nan_index);

% Calculate those coordinates and write them into the Matrix designating
% them for later estimation
Data_to_estimate=zeros(2,size(Nan_index,1));
for k=1:size(Nan_index,1)
    [x,y]=index_to_coord(row_ind(k),col_ind(k));
    Data_to_estimate([1,2],k)=[x,y];
end






end















