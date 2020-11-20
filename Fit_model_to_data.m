function [semivar_function]=Fit_model_to_data(experimental_semivariogram,var_type)

% [semivar_function]=Fit_model_to_data(experimental_semivariogram)
%
% This function takes as input the experimental semivariogram consisting of
% the lags and associated semivariance estimates. It outputs a model for
% the semivariogram function. This model is found by assuming that the
% underlying true semivariogram is of type "exponential" and then
% by finding those model parameters that minimize an objective function
% similar to the squared residuals.
%
% The formats are:
%                               
% experimental_semivariogram  = [Lag_bin1 .......  lag_binn_bin]   (2,n_bin)
%                               [semi_bin1  ....  semi_binn_bin]
%        variogram_type       =  Type of variogram function: either 
%                                'exponential' , 'squared_exponential' or 
%                                 'spherical'
%        semivar_function     =  function handle

% Take internal short name and delete empty bins
exp_sem=experimental_semivariogram;
exp_sem(:,isnan(exp_sem(2,:)))=[];
lag_exp=exp_sem(1,:);

% First define the model to fit into our experimental semivariogram. This 
% Function takes as input a row vector params =[range,sill] and the lag at
% which to evaluate the function.

if strcmp(var_type,'exponential')==1
semivar_fun_init = @(params,lag) params(2)*(1-exp(-abs(lag)./params(1)));
elseif strcmp(var_type,'squared_exponential')==1
semivar_fun_init = @(params,lag) params(2)*(1-exp(-(lag./params(1)).^2));
elseif strcmp(var_type,'spherical')==1
semivar_fun_init = @(params,lag) params(2)*((3*lag)./(2*params(1))-(lag.^3)./(2*params(1).^3)).*(lag<params(1))+params(2).*(lag>=params(1));
else
    disp('Unknown Convariance function. Choose either "exponential" or "squared exponential"')
end

% Prepare everything for optimization step. This includes finding some 
% initial starting parameters somewhat close to the solution and defining
% the objective function to be minimized.

% Initial value sill and range
params(2)=max(exp_sem(2,:));                       % just take the maximum
params(1)=exp_sem(1,find(exp_sem(2,:)>=0.5*params(2),1,'first')); % take
                                  % lag, where first hit 0.5 the maximum

                                  
% Please note that there are better and stochastically justified of fitting
% for semivariance functions than simple least squares. More professional
% are the versions in the textbooks from Cressie or Chiles & Delfiner


objectfun = @(params)...
    sum(((semivar_fun_init(params,lag_exp)-exp_sem(2,:)).^2));

% Use fminsearch to find those parameters [range,sill]=params_fin that
% minimize the objective function. Display error messages.
[params_fin,~,~,output] = fminsearch(objectfun,params);
disp(output)

% Define the model by using the optimizing parameters
if strcmp(var_type,'exponential')==1
semivar_function= @(lag) params_fin(2)*(1-exp(-abs(lag)./params_fin(1)));
elseif strcmp(var_type,'squared_exponential')==1
semivar_function= @(lag) params_fin(2)*(1-exp(-(lag./params_fin(1)).^2));
elseif strcmp(var_type,'spherical')==1   % Not correct; needs to have fixed limit
semivar_function = @(lag) params_fin(2)*((3*lag)./(2*params_fin(1))-(lag.^3)./(2*params_fin(1).^3)).*(lag<params_fin(1))+params_fin(2).*(lag>=params_fin(1));
else
    disp('Unknown Convariance function. Choose either "exponential" or "squared exponential"')
end
end