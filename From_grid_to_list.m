function [Data]=From_grid_to_list(Grid_data)

% [Data]=From_grid_to_list(Grid_data)
%
% This simple auxiliary function takes as input a Matrix Grid_data
% containing mostly NaNs and at some points measurement values.
% It will extract the row and column indices as well as the values of the
% measurements and put them into a List "Data", which can serve as input to
% Kriging_PPE.m
% 
% The formats are:
%       Grid_data            = [NaN .......         15]          (n,m)    
%                              [NaN .......        NaN]
%                              [12  .......        NaN]
%                              [NaN ......11       NaN]
%       Data                 = [X_p1 .......        X_pn_data]   (3,n_data)    
%                              [Y_p1 .......        Y_pn_data]
%                              [Z_p1 .......        Z_pn_data]

[rows, cols]=find(~isnan(Grid_data));
vals=zeros(size(rows));
for k=1:size(rows,1)
    vals(k)=Grid_data(rows(k),cols(k));
end
Data= [rows'; cols'; vals'];



end