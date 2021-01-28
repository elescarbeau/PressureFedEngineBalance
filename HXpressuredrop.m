% Heat Exchanger Pressurant Pressure Drop Calculation Function
% Author: Elysse Lescarbeau, Date: 12/30/20
% Equations/methodology from Heat Exchanger Dimensioning by J. Saari,
% p.84-87 (https://drive.google.com/file/d/1DAyVj1EhW7SGOJYqiDseTSWOlxGQSsTN/view?usp=sharing)

function [P_2,rho_2,dP,friction_coeff] = HXpressuredrop(nominal_pressurant_mdot,num_tubes,tube_length,inner_diameter,outer_diameter,tube_roughness,P_1,rho_1,T_2,viscosity)
%% units/universal values
% mdot: lbm/s
% tube_length: in
% inner_diameter: in 
% tube_roughness: in
% P_1, P_2: psia
% rho_1, rho_2: lbm/ft^3
% T_1, T_2: F
% viscosity: lbm/ft*s, ASSUMED CONSTANT UNTIL COOLPROP DOWNLOADED
% dP: psid

grav_constant = 32.17405; % ft/s^2, gravitational constant
nitrogen_gas_constant = 55.165; % ftlbf/lbmR, specific gas constant, assumed constant

%% pressure drop calculation

Re = 4*(nominal_pressurant_mdot/num_tubes)/(pi*(inner_diameter/12)*viscosity); % Reynold's Number

friction_coeff = 0.0001; % from https://www.engineersedge.com/fluid_flow/pressure_drop/pressure_drop.htm
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

fanning_friction_coeff = friction_coeff/4;
sigma = (pi/4*inner_diameter^2)/(pi/4*outer_diameter^2);
entrance_loss_coeff = .5; % from https://www.engineeringtoolbox.com/minor-loss-coefficients-pipes-d_626.html
exit_loss_coeff = 1; % from https://www.engineeringtoolbox.com/minor-loss-coefficients-pipes-d_626.html

P_2_guess = P_1; % psia
rho_2_guess = (P_2_guess*144)/(nitrogen_gas_constant*(T_2+459.67)); % lbm/ft^3
rho_average_guess = (rho_1+rho_2_guess)/2; % lbm/ft^3
convergence = 0;
while convergence == 0
    dP = (((((nominal_pressurant_mdot/num_tubes)/(pi/4*inner_diameter^2))^2)/(2*rho_1))*(12^3)/grav_constant/12)*...
        ((1-sigma^2+entrance_loss_coeff)+2*(rho_1/rho_2_guess-1)+(rho_1/rho_average_guess)*...
        (fanning_friction_coeff*(4*tube_length/inner_diameter)+(entrance_loss_coeff+exit_loss_coeff))+...
        (rho_1/rho_2_guess)*(exit_loss_coeff+sigma^2-1)); % psid
    P_2 = P_1 - dP;
    if (abs(P_2 - P_2_guess)/P_2_guess)*100 < 1
        convergence = 1;
    else
        P_2_guess = P_2_guess - .01; % psia
        rho_2_guess = (P_2_guess*144)/(nitrogen_gas_constant*(T_2+459.67)); % lbm/ft^3
        rho_average_guess = (rho_1+rho_2_guess)/2; % lbm/ft^3
    end
end
    
rho_2 = rho_2_guess;


end