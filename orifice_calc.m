% Orifice diameter calculator
% Elysse Lescarbeau, 2/5/20
% update log:
% 4/17/20: added if statement for compressible or incompressible, added calculation of density for compressible fluid based on upstream temp andpressure
% 12/30/20: converted script into function
% Reference: Crane Flow of Fluids, 2010 version, pg. 4-6

function [orifice_diameter,dP] = orifice_calc(mdot_desired,tube_inner_diameter,P_1,P_2,rho_1,fluid_type)
%% units/universal values

% mdot_desired: lbm/s
% tube_inner_diameter: in
% P_1, P_2: psia
% rho_1: lbm/ft^3

spec_heat_ratio = 1.4; % for nitrogen, dimensionless, assumed constant
dP = P_1-P_2;

%% switch for fluid type to calculate required orifice diameter

% fluid_type: if 3, compressible; if 4, incompressible

convergence = 0;
orifice_diameter_guess = .000001;% orifice diameter initial guess, lbm/s 

while convergence == 0
    beta = orifice_diameter_guess/tube_inner_diameter;
    Y = 1 - (.351 + .256*(beta^4) + .93*(beta^8)*(1-((P_2/P_1)^(1/spec_heat_ratio))));
    Cd = .65; % coefficient of discharge for square edged orifice
    C = Cd/sqrt(1-(beta^4));
    if fluid_type == 3
        mdot = Y*C*(((orifice_diameter_guess/12)^2)*(pi()/4))*sqrt(2*rho_1*32.2*(dP*144)); % compressible mdot
    else
        delta_p_immediate = dP*((sqrt(1-(beta^4)*(1-(Cd^2)))+Cd*(beta^2))/(sqrt(1-(beta^4)*(1-(Cd^2)))-Cd*(beta^2))); % accounts for recoverable pressure drop
        mdot = C*(((orifice_diameter_guess/12)^2)*(pi()/4))*sqrt(2*rho_incomp*32.2*(delta_p_immediate*144)); % incompressible mdot
    end
    if (abs(mdot_desired-mdot)/mdot)*100 < .1
        convergence = 1 ;
    else
        orifice_diameter_guess=orifice_diameter_guess+.00001;
    end
end

orifice_diameter = orifice_diameter_guess;

end


