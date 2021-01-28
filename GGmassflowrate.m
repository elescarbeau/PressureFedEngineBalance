% Author: Elysse Lescarbeau, last updated: 1/6/21
% This function is based off the code im "HEx_massflowrate_final.m" but is
% updated to use as a function and to incorporate CoolProp to automatically
% pull fluid statepoints
function [gg_ox_mdot,gg_fuel_mdot,gg_total_mdot] = GGmassflowrate(GG_OF_ratio,...
    engine_ox_mdot,engine_fuel_mdot,ox_density,fuel_density,pressure_ox_tank,...
    pressure_fuel_tank,HX_pressurant_inlet_pressure,HX_pressurant_outlet_pressure,...
    HX_pressurant_inlet_temp,HX_pressurant_outlet_temp,tank_pressurant_temp,...
    HX_EG_inlet_temp,HX_EG_outlet_temp,R_EG_SI,gamma_EG)

% units %
% inputs 
% GG_OF_ratio: dimensionless
% engine_ox_mdot, engine_fuel_mdot: lbm/s
% ox_density, fuel_density: lbm/ft^3
% pressure_ox_tank, pressure_fuel_tank, HX_pressurant_inlet_pressure, HX_pressurant_outlet_pressure: psia
% HX_pressurant_inlet_temp, HX_pressurant_outlet_temp,tank_pressurant_temp, HX_EG_inlet_temp, HX_EG_outlet_temp: F
% R_EG_SI: J/kgK
% gamma_EG: dimensionless
% outputs
% gg_ox_mdot,gg_fuel_mdot,gg_total_mdot: lbm/s

HX_pressurant_inlet_temp_SI = (HX_pressurant_inlet_temp - 32)* 5/9 + 273.15293; % F to K
HX_pressurant_outlet_temp_SI = (HX_pressurant_outlet_temp - 32)* 5/9 + 273.15293; % F to K  
HX_EG_inlet_temp_SI = (HX_EG_inlet_temp - 32)* 5/9 + 273.15293; % F to K
HX_EG_outlet_temp_SI = (HX_EG_outlet_temp - 32)* 5/9 + 273.15293; % F to K
tank_pressurant_temp_SI = (tank_pressurant_temp - 32)* 5/9 + 273.15293; % F to K

cp_N2_inlet = py.CoolProp.CoolProp.PropsSI('CP0MASS','P',HX_pressurant_inlet_pressure*6894.76,...
    'T',HX_pressurant_inlet_temp_SI,'Nitrogen');% J/kgK
cp_N2_outlet = py.CoolProp.CoolProp.PropsSI('CP0MASS','P',HX_pressurant_outlet_pressure*6894.76,...
    'T',HX_pressurant_outlet_temp_SI,'Nitrogen'); % J/kgK
cp_N2 = (cp_N2_inlet + cp_N2_outlet)/2; % J/kgK

R_N2_SI = 296.80; % J/kgK, nitrogen gas constant, assumed constnt
cp_EG_SI = (gamma_EG*R_EG_SI)/(gamma_EG-1) ; % J/kgK

ox_density_SI = ox_density*16.0185; % lbm/ft^3 to kg/m^3
fuel_density_SI = fuel_density*16.0185; % lbm/ft^3 to kg/m^3
pressure_ox_tank_SI = pressure_ox_tank*6894.757; % psia to Pa
pressure_fuel_tank_SI = pressure_fuel_tank*6894.757; % psia to Pa
pressurant_ox_density_SI =  pressure_ox_tank_SI/(R_N2_SI*tank_pressurant_temp_SI);% kg/m^3
pressurant_ox_density = pressurant_ox_density_SI/16.0185; % kg/m^3 to lbm/ft^3
pressurant_fuel_density_SI =  pressure_fuel_tank_SI/(R_N2_SI*tank_pressurant_temp_SI);% kg/m^3 
pressurant_fuel_density = pressurant_fuel_density_SI/16.0185; % kg/m^3 to lbm/ft^3
engine_ox_pressurant_volflowrate = engine_ox_mdot/ox_density; % ft^3/s
engine_fuel_pressurant_volflowrate = engine_fuel_mdot/fuel_density; % ft^3/s
engine_ox_pressurant_mdot = engine_ox_pressurant_volflowrate*pressurant_ox_density; % lbm/s
engine_fuel_pressurant_mdot = engine_fuel_pressurant_volflowrate*pressurant_fuel_density; % lbm/s
engine_pressurant_mdot = engine_ox_pressurant_mdot + engine_fuel_pressurant_mdot; % lbm/s
engine_pressurant_mdot_SI=engine_pressurant_mdot/2.205; %kg/s; 

syms x % x=mdot_gg_fuel_SI
equ = (GG_OF_ratio*x+x)*(cp_EG_SI*(HX_EG_inlet_temp_SI-HX_EG_outlet_temp_SI))/(cp_N2*(HX_pressurant_outlet_temp_SI-HX_pressurant_inlet_temp_SI)) == engine_pressurant_mdot_SI + ((GG_OF_ratio*x)/ox_density_SI)*pressurant_ox_density_SI + (x/fuel_density_SI)*pressurant_fuel_density_SI;
gg_fuel_mdot_SI = solve(equ,x);
gg_fuel_mdot_SI = double(gg_fuel_mdot_SI);
gg_ox_mdot_SI = GG_OF_ratio*gg_fuel_mdot_SI;
gg_total_mdot_SI=gg_fuel_mdot_SI+gg_ox_mdot_SI;
gg_fuel_mdot = gg_fuel_mdot_SI*2.205;
gg_ox_mdot = gg_ox_mdot_SI*2.205;
gg_total_mdot=gg_total_mdot_SI*2.205;

end

