function k = kFunc(h,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The k function for the model.
%
% INPUTS:
% h = pressure at the node.
% z = height of the node.
% OUTPUTS:
% k = k evaluated at the node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = ConstChooser('n',z);
m = 1 - 1/n;

if h < 0
    k = sqrt(SFunc(h,z)) * (1 - (1 - SFunc(h,z)^(1/m) )^m )^2;
else
    k = 1;
end

end


function x = ConstChooser(string,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the correct constant for a node in the model
%
% INPUTS:
% z      = the height of the node
% string = what constant is needed, takes:
%           'sat', 'res', 'alpha', and 'n'
% OUTPUTS:
% x      = the constant need.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch string
    case 'sat'
        if z < 50
            x = .1;
        else
            x = .33;
        end
    case 'res'
        if z < 50
            x = 0;
        else
            x = 0;
        end
    case 'alpha'
        if z < 50
            x = 1.43;
        else
            x = 1.43;
        end
    case 'n'
        if z < 50
            x = 1.51;
        else
            x = 1.51;
        end
    otherwise
        error('Incorrect String: Only can input sat, res, slpha, n.   Error @ ConstChooser Function')
end
end

