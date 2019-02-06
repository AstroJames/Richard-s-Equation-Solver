function StatisticsandVisualisation(hStorage,runTime,timeSteps,FEvals,tInterval,z,espis,type,choice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes in the output from the model and performs either
% visualisations or statistics on the output depending on the user inupt
%
% INPUTS:
%   hStorage    = the matrix of time iterations across the grid
%   runTime     = the runtime for the model
%   timeSteps   = the number of time steps the model took
%   FEvals      = the number of evals of F(x)
%   tInterval   = the time domain whihch the model was solved over
%   z           = the z grid from the 1D problem
%   espis       = the water budget residuals for every time step
%   type        = type specifies whether you want to visualise the solution
%                   or look at the statistics of the solution
%   choice      = chooses the type of visualisation for the
%                   visualisation type
%
% OUTPUTS:
%   the function outputs either statistics and visualisations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first switch case chooses between visualisation and statistics.
switch type
    case 'Visualisation'
        
        switch choice 
            % The second switch case chooses between what profile to
            % visualise.
            case 'PressureHead'
                disp('You have selected to visualise the pressure head.')
                figure
                for i = 1:size(hStorage,2)
                    drawnow
                    pause(0.01)
                    plot(hStorage(:,i)',z)
                    xlim([-10,80])
                    ylim([0,80])
                    formatSpec = 'Pressure Head: Current Time Step is %.1f.';
                    title(sprintf(formatSpec,i));
                    ylabel('Z axis (m)');
                end
            case 'WaterContent'
                disp('You have selected to visualise the water content.')
                phi0    = zeros(1,length(z));
                for j = 1:size(hStorage,2)
                    for i = 1:length(z)
                        phi0(i) = phiFunc(hStorage(i,j)',z(i));
                    end
                    drawnow;
                    plot(phi0,z);    
                    xlim([0,0.35]);
                    ylim([0,80]);
                    formatSpec = 'Water Content: Current Time Step is %.1f.';
                    title(sprintf(formatSpec,j));
                    ylabel('Z axis (m)');
                end
            case 'Saturation'
                disp('You have selected to visualise the saturation.')
                SFunc0    = zeros(1,length(z));
                for j = 1:size(hStorage,2)
                    for i = 1:length(z)
                        SFunc0(i) = SFunc(hStorage(i,j)',z(i));
                    end
                    drawnow
                    plot(SFunc0,z);
                    xlim([0,1])
                    formatSpec = 'Saturation: Current Time Step is %.1f.';
                    title(sprintf(formatSpec,j));
                    ylabel('Z axis (m)');
                end

            case 'Permeability'
                disp('You have selected to visualise the permability.')
                kFunc0    = zeros(1,length(z));
                for j = 1:size(hStorage,2)
                    for i = 1:length(z)
                        kFunc0(i) = kFunc(hStorage(i,j)',z(i));
                    end
                    drawnow
                    plot(kFunc0,z);
                    xlim([0,1])
                    pause(0.01)
                    formatSpec = 'Permability: Current Time Step is %.1f.';
                    title(sprintf(formatSpec,j));
                    ylabel('Z axis (m)');
                end
        end
        
    case 'Statistics'
        disp('You have selected to perform statistics on the model output')
        check = inf;
        zAtNum = zeros(2,1);
        zAtNum(1) = 1;
        % Estimate where the water table is
        while abs(70 - z(zAtNum(1))) < check
            check = 70 - z(zAtNum(1));
            zAtNum(1) = zAtNum(1)+1;
        end
        zAtNum(2) = z(zAtNum(1));
        
        % Count the number of times the water table is about 70.
        counter = 0;
        for i = 1:length(hStorage)
            if hStorage(zAtNum(1),i) > 0 
                counter = counter + 1;
            end
        end
        
        
        % Calculate the average error in the water budget over the time
        % interval
        period = tInterval(2)-tInterval(1);
        epsBar = 0;
        for i = 1:length(espis)
            epsBar = epsBar + timeSteps(i) * espis(i)^2;
        end
        epsBar = epsBar/period;
end
        
 % Output all of the statistics.
fprintf('The number of F(x) evaluations was %.0f \n',sum(FEvals))
fprintf('The number of F(x) evaluations per time iterations was %.3f \n',sum(FEvals)/sum(timeSteps))
fprintf('The water table stayed above z = 70m %.4f percent of the time \n',counter/(tInterval(2)))
fprintf('The water budget residual is %.10f \n',epsBar)
fprintf('The model runtime was %.3f seconds \n',runTime)

end


