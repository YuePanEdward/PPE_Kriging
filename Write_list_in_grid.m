function [Filled_grid]=Write_list_in_grid(Grid,coord_to_index,Data_list,Row_data)

% [Filled_grid]=Write_list_in_grid(Grid,coordinate_to_index,Data_list)
% 
% This function takes as input a list consisting of coordinates and values,
% a grid into which these values should be written as well as the necessary
% function (handle) that performs the mapping from coordinates to indices.
% This function is normally made by "From_list_to_grid".
% Row data is just an index specifying in which row of the Data_list the
% actual values lie that should be mapped on the grid.
%
% The output of this function is just a grid "Filled_grid" of exactly the
% same dimension as the Grid but now with additionally the values injected,
% that are specified in Data_list.
%
% The formats are:
%       Grid                 = [NaN .......          Z_p95   ]        (n,m)    
%                              [NaN ... Z_p43             NaN]
%       Data_list            = [X_s1 .......         X_sn_est]   (3,n_data)    
%                              [Y_s1 .......         Y_sn_est]
%                              [Z_s1 .......         Z_sn_est]
%                              [sigma_s1 .......  sigma_n_est]
%       Filled_Grid          = [Z_s1 ....   Z_s12   ... Z_p95]        (n,m)
%                              [Z_s24 ....   Z_p43  Z_sn_est ]
%       coordinate_to_index  = function handle


% Extract coordinates from Data_list and use mapping to find corresponding
% indices.
index_list=zeros(size(Data_list,2),2);
for k=1:size(Data_list,2)
    [i,j]=coord_to_index(Data_list(1,k),Data_list(2,k));
    index_list(k,[1,2])=[i,j];
end

% Fill the grid by accessing the elements corresponding to data points and
% filling them with the corresponding data
Filled_grid=Grid;
Filled_grid(sub2ind(size(Grid),index_list(:,1),index_list(:,2)))=Data_list(Row_data,:);


end