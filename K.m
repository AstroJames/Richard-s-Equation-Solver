function k = K(x,z,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K retrieves the conductivity of material 'r' at the coordinate (x,z)
% within the unit cells.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% x         =   the x coordinate of the material on the unit cell.
% z         =   the z coordinate of the material on the unit cell.
% r         =   the type of the material of interest, where
%               - 'a' is the Alluvium Aquifer,
%               - 'c' is the Walloon Coal Measure,
%               - 't' is the benchmark test case.
% Output:
% k         =   the conductivity of material r at the coordinate (x,z).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check domain
if x > 1 || x < 0 || z > 1 || z < 0
    error('Coordinate lies outside of defined domain')
end

switch r
    case 'a' % Alluvium Aquifer
        if x < 0.75 && z >= 0.75 % Soil A above
            k = 6.07;
        elseif x >= 0.25 && (z >= 0.25 && z < 0.5) % Soil A below
            k = 6.07;
        else % Soil B
            k = 3.69;
        end
    case 'c' % Walloon Coal Measure
        if x < 0.25 || (x >=0.5 && x < 0.75) % 1st and 3rd column
            if z < 0.25 || (z >=0.5 && z < 0.75) % 1st and 3rd row
                k = 2.32*10^-1; % Soil A
            else % 2nd and 4th row
                k = 8.40*10^-1; % Soil B
            end
        else % 2nd and 4th column
            if z < 0.25 || (z >=0.5 && z < 0.75) % 1st and 3rd row
                k = 8.40*10^-1; % Soil B
            else
                k = 2.32*10^-1; % Soil A
            end
        end
    case 't' % Benchmark Problem
        if z < 0.5 % Soil A
            k = 16; % Can change value here
        else % Soil B
            k = 4;  % Can change value here
        end
        
        % Note: Theoretical Answer:
        % Dzz = (2*K_A*K_B)/(KA + KB)
        % Dxx = (K_A + K_B)/2
    otherwise
        error('Unknown input material type.')
end

end

