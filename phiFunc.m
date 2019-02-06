function phi = phiFunc(h,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The moisture content function for the model.
%
% INPUTS:
% h = pressure at the node.
% z = height of the node.
% OUTPUTS:
% phi = phi evaluated at the node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sat = ConstChooser('sat',z);
res = ConstChooser('res',z);

if h < 0
    phi = res + SFunc(h,z) * (sat - res);
else
    phi = sat;
end

end

function S = SFunc(h,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The S function for the model.
%
% INPUTS:
% h = pressure at the node.
% z = height of the node.
% OUTPUTS:
% S = S evaluated at the node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha   = ConstChooser('alpha',z);
n       = ConstChooser('n',z);
m       = 1 - 1/n;

if h < 0
    S  = (1 + (-alpha*h)^n)^-m;
else
    S = 1;
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
            x = 0.1;
        else
            x = 0.33;
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