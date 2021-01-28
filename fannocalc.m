%% Fanno Flow Line Calculations Function
% Author: Elysse Lescarbeau, Date: 12/17/2020
% Equations/methodology from Modern Compressible Flow, 2nd Edition by John
% D. Anderson. Jr., p. 85 (esp Ex 3.12)

function [P_2,rho_2,dP,T_2,v_2,M_2,dT] = fannocalc(mdot,P_1,T_1,rho_1,tube_length,inner_diameter)
%% units/universal values
% mdot: lbm/s
% P_1, P_2: psia
% T_1, T_2: F
% rho_1, rho_2: lbm/ft^3
% tube_length: in
% inner_diameter: in

nitrogen_gas_constant = 55.165; % ftlbf/lbmR, assumed constant
gamma = 1.4; % dimensionless, specific heat ratio, assumed constant
f = .005; % friction coefficient, assumed constant

%% upstream flow conditions
a_1 = (pi/4)*(inner_diameter^2); % in^2, initial tubing inner cross sectional area
v_1 = mdot/(rho_1*(a_1/144)); % ft/s, initial velocity
M_1=v_1/sqrt(gamma*nitrogen_gas_constant*32.2*(T_1+459.67)); % upstream Mach number
charlengthval_L_1_star = (1-M_1^2)/(gamma*M_1^2)+(gamma+1)/(2*gamma)*log(((gamma+1)*M_1^2)/(2+(gamma-1)*M_1^2)); % Reference Equ 3.107
charlengthval_L = (4*f*tube_length)/inner_diameter; 
charlengthval_L_2_star = charlengthval_L_1_star - charlengthval_L;  % see Ex 3.12

%% downstream flow conditions

% downstream Mach number
M_2 = 1; % dimensionless, downstream Mach number, initial guess
convergence = 0;
while convergence == 0
    charlengthval_L_2_star_calculated = (1-M_2^2)/(gamma*M_2^2)+(gamma+1)/(2*gamma)*log(((gamma+1)*M_2^2)/(2+(gamma-1)*M_2^2));
    if (abs(charlengthval_L_2_star_calculated-charlengthval_L_2_star)/charlengthval_L_2_star)*100 < 1
        convergence = 1 ;
    else
        M_2 = M_2 - .0001;
    end
end

% downstream pressure
P_2 =((M_1/M_2)*sqrt((2+(gamma-1)*M_1^2)/(2+(gamma-1)*M_2^2)))*P_1; % psia
dP = P_1-P_2; % psid

% downstream temperature
T_2 = (((2+(gamma-1)*M_1^2)/(2+(gamma-1)*M_2^2))*(T_1+459.67))-459.67;
dT = T_1-T_2; % F

% downstream density
rho_2 = (P_2*144)/(nitrogen_gas_constant*(T_2+459.67)); % lbm/ft^3

% downstream velocity
v_2 = mdot/(rho_2*(a_1/144)); % ft/s, velocity

end