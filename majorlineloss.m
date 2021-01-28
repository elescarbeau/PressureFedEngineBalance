% Incompressible Major Line Loss  Calculations Function
% Author: Elysse Lescarbeau, Date: 12/17/2020
% Equations/methodology from www.engineersedge.com/fluid_flow/pressure_drop/pressure_drop.htm

function [dP,v_1] = majorlineloss(mdot,tube_length,inner_diameter,tube_roughness,rho_1,viscosity)

% mdot: lbm/s
% tube_length: in
% inner_diameter: in 
% tube_roughness: in
% rho_1: lbm/ft^3
% viscosity: lbf*s/ft^2
% P_2: psia
% dP: psid

grav_constant = 32.17405; % ft/s^2, gravitational constant

Re = 4*(mdot/grav_constant)/(pi*(inner_diameter/12)*viscosity); % Reynold's Number

friction_coeff = 0.0001; % Darcy coefficient
convergence = 0;
while convergence == 0
    term1 = 1/sqrt(friction_coeff);
    term2 = -2*log10(2.51/(Re*sqrt(friction_coeff))+(tube_roughness/inner_diameter)*0.269);
        if abs(term1 - term2) < .01
            convergence = 1;
        else
            friction_coeff = friction_coeff + 0.0001;
        end
end

A_1 = pi/4*(inner_diameter/12)^2; % ft^2
v_1 = mdot/(rho_1*A_1); % ft/s
dP = (friction_coeff*(tube_length/inner_diameter)*(rho_1/2)*v_1^2)/(grav_constant*144); % psia
end
