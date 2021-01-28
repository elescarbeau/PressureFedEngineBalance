%% Component Pressure Drop Calculations Function
% Author: Elysse Lescarbeau, Date: 12/17/2020
% Equations/methodology from Aerospace Fluids Component Handbook
% NOTE: Temperature assumed constant across components

function [P_2,rho_2,dP] = componentdP(mdot,flow_area,CvOrCdA,rho_1,varargin)
%% units

% varargin variables in order: T_1,P_1
% mdot: lbm/s
% flow_area: if CdA: in^2, if Cv: dimensionless
% rho_1: lbm/ft^3
% T_1: F
% P_1: psia

%% universal values
nitrogen_gas_constant = 55.165; % ftlbf/lbmR, specific gas constant, assumed constant
grav_constant = 32.17405; % ft/s^2, gravitational constant
spec_heat_ratio = 1.4; % dimensionless, nitrogen specific heat ratio, assumed constant

%% flow area conversion
% if CvOrCdA == 1, Cv is given for flow_area and needs to be converted to CdA
% if CvOrCdA == 2, CdA is given for flow_area and no conversion is needed

if CvOrCdA == 1
    flow_area = pi/4*0.65*((.227*sqrt(flow_area))^2);
else
    flow_area = flow_area;
end

%% switch for fluid type to calculate dP
% incompressible: n=4
% compressible, not choked: n=6
    
n = nargin;

if n == 4
    dP = ((mdot/flow_area)^2)/(2*((rho_1/12^3)*grav_constant*12)); % incompressible pressure drop 
    P_2 = 0; % not needed
    rho_2 = rho_1; % density assumed constant for incompressible flow
else
    T_1 = varargin{1};
    P_1 = varargin{2};
    convergence = 0;
    P_2 = P_1 - .001;
    while convergence == 0
        mdot_calculated = (flow_area/144)*rho_1*((P_2/P_1)^(1/spec_heat_ratio))...
            *sqrt(2*grav_constant*nitrogen_gas_constant*(T_1+459.67)*(spec_heat_ratio/(spec_heat_ratio-1))...
            *(1-(P_2/P_1)^((spec_heat_ratio-1)/spec_heat_ratio))); 
        if (abs(mdot_calculated - mdot)/mdot)*100 < 1
            convergence = 1 ;
        else
            P_2 = P_2 - 0.001;
        end
    end
    dP = P_1-P_2;
    rho_2 = (P_2*144)/(nitrogen_gas_constant*(T_1+459.67)); % lbm/ft^3
end

end