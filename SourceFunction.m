function [S_weather_1,S_weather_2]=SourceFunction(model,t,Delta,droughtState)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% MODEL         - this specifies the weather model to use: {High Frequency, Low Frequency}. 
% T             - this is the value of time to evaluate the weather model at. 
% DELTA         - this takes the cell array of differences between the nodes, the cells and time.
% DROUGHTSTATE  - this input tells the source to be in a drought state or not. 

% OUTPUTS:
% S_WEATHER_1   - this is the source function for the rain evaluated at the current time step. 
% S_WEATHER_2   - this is the source function for the rain evaulated at the next time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sum(strcmp(model, 'High Frequency') || strcmp(model, 'Low Frequency') || strcmp(model,'Low Frequency Stochastic') || strcmp(model,'Simple')) ~= 1
    error('Please select either High Frequency, Low Frequency or Drought for the model.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delt    = Delta{3}; % the time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAIN MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch model
   
    % This is for the coarse function that describes the rainfall.
    case 'High Frequency'
        data=importdata('High_frequency_data.csv'); % import all of the coefficients calculated in R
        parms=data.data(1:19);
        % the functional form low frequency
        f_weather = @(t) parms(1) + parms(2)*t + parms(4)*sin(1/2*pi*t) + parms(5)*cos(1/2*pi*t) + parms(6)*sin(1/4*pi*t) + ...
           parms(7)*cos(1/4*pi*t) + parms(8)*sin(1/6*pi*t) + parms(9)*cos(1/6*pi*t) + parms(10)*sin(1/8*pi*t) + ...
           parms(11)*sin(1/10*pi*t) + parms(12)*sin(1/2000*pi*t) + parms(13)*cos(1/20000*pi*t) + ...
           parms(14)*sin(1/40000*pi*t) + parms(15)*cos(1/80000*pi*t) + parms(16)*sin(1/120000*pi*t) + parms(17)*sin(1/150000*pi*t) + parms(18)*cos(1/180000*pi*t) + ...
           parms(19)*sin(1/190000*pi*t);
        S_weather_1 = (-1)*f_weather(t)*10^-3;
        S_weather_2 = (-1)*f_weather((t+delt))*10^-3;

    
    % This is for the smooth function that describes the rainfall.
    case 'Low Frequency'
        data=importdata('Low_frequency_data.csv');% import all of the coefficients calculated in R
        parms = data.data(1:10);
        data.textdata(9:14);
        % the functional form of low frequency
        f_weather = @(t) parms(1) +parms(2)*t + parms(3)*sin(1/2000*pi*t) + parms(4)*cos(1/20000*pi*t) +...
            parms(5)*sin(1/40000*pi*t) + parms(6)*cos(1/80000*pi*t) +parms(7)*sin(1/120000*pi*t)+ ...
            parms(8)*sin(1/150000*pi*t)+parms(9)*cos(1/180000*pi*t)+parms(10)*sin(1/190000*pi*t);
        S_weather_1 = (-1)*f_weather(t)*10^-3;
        S_weather_2 = (-1)*f_weather(t+delt)*10^-3;
        
   % This calls the stochastic rainfall function. 
   case 'Low Frequency Stochastic'
        data=importdata('Low_frequency_data.csv'); % import all of the coefficients calculated in R
        parms = data.data(1:10);
        data.textdata(9:14);
        f_weather = @(t) parms(1) +parms(2)*t + parms(3)*sin(1/2000*pi*t) + parms(4)*cos(1/20000*pi*t) +...
            parms(5)*sin(1/40000*pi*t) + parms(6)*cos(1/80000*pi*t) +parms(7)*sin(1/120000*pi*t)+ ...
            parms(8)*sin(1/150000*pi*t)+parms(9)*cos(1/180000*pi*t)+parms(10)*sin(1/190000*pi*t);
        
        % add a normally distributed amplitude of the rainfall by sampling
        % at the mean of the low frequency model.
        S_weather_1 = normrnd(f_weather(t),20000)*10^-3;
        % check if the result is less than zero, if it is then make it 0.
        if S_weather_1 < 0 
            S_weather_1 = 0;
        else
            S_weather_1 = (-1)*S_weather_1;
        end
        
        % do the same check for the time + delt step
        % check if the result is less than zero, if it is then make it 0.
        S_weather_2 = normrnd(f_weather(t+delt),20000)*10^-3;
        if S_weather_2 < 0 
            S_weather_2 = 0;
        else
            S_weather_2 = (-1)*S_weather_2;
        end
        
        
    case 'Simple'    
    f_weather      = @(x) -((1.65*10^-3)+(1.65*10^-3)*cos(2*pi*x/365));
    S_weather_1    = (-1)*f_weather(t);
    S_weather_2    = (-1)*f_weather(t+delt);

end

% check the drought state
switch droughtState
    case 'FullyIn'
        S_weather_1 = 0;
        S_weather_2 = 0;
    case 'GoingOut'
        S_weather_2 = 0;
    case 'FullyOut' 
        
    case 'GoingIn'
        S_weather_1 = 0;
    otherwise
        eString = strcat('String "',droughtState,'" is not accepted');
        error(eString)
end


end % end function





