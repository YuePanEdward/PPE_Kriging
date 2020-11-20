function [Nearest_Neighbor, Linear_Interp, Natural_Neighbor, Cubic]=Make_comparative_interpolation(Data, Field_values_init)


% 1: Convert Data list to grid and back to get rows and cols instead of
% coordinates
Data_grid=From_list_to_grid(Data,Field_values_init);
Data_list=From_grid_to_list(Data_grid);

y=Data_list(1,:);
x=Data_list(2,:);
z=Data_list(3,:);

% 2: Apply interpolation functions
[X, Y] = meshgrid(1:size(Data_grid,2),1:size(Data_grid,1));
Cubic = griddata(x,y,z,X,Y,'cubic');
Linear_Interp=griddata(x,y,z,X,Y,'linear');
Natural_Neighbor=griddata(x,y,z,X,Y,'natural');
Nearest_Neighbor=griddata(x,y,z,X,Y,'nearest');
end