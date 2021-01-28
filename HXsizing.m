function [nominal_pressurant_mdot,overall_heat_transfer_coeff,surface_area_ft2,...
    tube_length,space_btwn_tubes,HX_pressurant_avg_viscosity,overall_heat_transfer_rate,HX_pressurant_thermR_conv,HX_pressurant_thermR_rad,HX_EG_thermR_conv,HX_EG_thermR_rad,HX_EG_avg_outertubing_temp,HX_pressurant_avg_innertubing_temp,tube_length_guess,HX_pressurant_avg_Pr...
    ] = HXsizing(GG_total_mdot,EG_density,...
    EG_avg_viscosity,EG_avg_thermcond,EG_avg_spec_heat,pressurant_mdot,num_tubes,inner_diameter,outer_diameter,...
    HX_pressurant_inlet_pressure,HX_pressurant_outlet_pressure_initial,...
    HX_pressurant_inlet_temp,HX_pressurant_outlet_temp,tank_pressurant_temp,...
    HX_EG_inlet_temp,HX_EG_outlet_temp)




%% units
% inputs
% GG_total_mdot, pressurant_mdot: lbm/s
% EG_density: lbm/s
% EG_avg_viscosity: lbm/ft*s
% EG_avg_thermcond: Btu/sftF
% EG_avg_spec_heat: Btu/lbmF
% inner_diameter, outer_diameter, tube_length: in
% HX_pressurant_inlet_pressure, HX_pressurant_outlet_pressure_initial: psia
% HX_pressurant_inlet_temp, HX_pressurant_outlet_temp, tank_pressurant_temp,
% HX_EG_inlet_temp, HX_EG_outlet_temp: F
% outputs
% nominal_pressurant_mdot: lbm/s
% overall_heat_transfer_coeff: btu/sft^2F
% surface_area_ft2: ft^2
% space_btwn_tubes: in

%% set up convergence values
convergence1 = 0;
convergence2 = 0;
convergence3 = 0;
convergence4 = 0;

tube_length_guess = 50; % in
HX_EG_avg_outertubing_temp_initial = 650; % F, initial guess for inner surface temp, cannot finalize this value until initial HX sizing is found 
HX_pressurant_avg_innertubing_temp_initial = 600; % F, initial guess for inner surface temp, cannot finalize this value until initial HX sizing is found 

%% Log Mean Temperature Difference
LMTD = ((HX_EG_inlet_temp-HX_pressurant_outlet_temp)-(HX_EG_outlet_temp-HX_pressurant_inlet_temp))...
    /log((HX_EG_inlet_temp-HX_pressurant_outlet_temp)/(HX_EG_outlet_temp-HX_pressurant_inlet_temp)); % LMTD for counterflow HX

%% solve for nominal HX pressurant mass flow rate to use bypass valve to bring HX pressurant outlet temp to correct temp
pressurant_mdot_bypassHX_fraction_initial = 0;
pressurant_mdot_enteringHX_fraction_initial = 1 - pressurant_mdot_bypassHX_fraction_initial;
while convergence1 == 0
    HX_pressurant_outlet_temp_mixed = (pressurant_mdot*pressurant_mdot_enteringHX_fraction_initial...
        *HX_pressurant_outlet_temp + pressurant_mdot*pressurant_mdot_bypassHX_fraction_initial...
        *HX_pressurant_inlet_temp)/pressurant_mdot; % HX outlet temperature after being mixed with bypass line
if (abs(HX_pressurant_outlet_temp_mixed - tank_pressurant_temp)/tank_pressurant_temp)*100 < 1
    convergence1 = 1;
else
    pressurant_mdot_bypassHX_fraction_initial = pressurant_mdot_bypassHX_fraction_initial + .001;
    pressurant_mdot_enteringHX_fraction_initial = 1 - pressurant_mdot_bypassHX_fraction_initial;
end
end

pressurant_mdot_bypassHX_fraction = pressurant_mdot_bypassHX_fraction_initial;
pressurant_mdot_enteringHX_fraction = pressurant_mdot_enteringHX_fraction_initial;
nominal_pressurant_mdot = pressurant_mdot*pressurant_mdot_enteringHX_fraction;

%% exhaust gas side convective heat transfer
% solve for space between tubes
velocity_btwn_tubes = 90; % ft/s, recommended max velocity 
csarea_btwn_tubes_ft2 = GG_total_mdot/(EG_density*velocity_btwn_tubes); % cross sectional area between stack of tubes
csarea_btwn_tubes = csarea_btwn_tubes_ft2*144; % in^2
space_btwn_tubes = ((csarea_btwn_tubes/num_tubes)-((outer_diameter^2)-(pi/4*(outer_diameter^2))))/outer_diameter; % in

% solve for thermal resistance
% set up overall loop

while (convergence2 == 0 || convergence3 == 0 || convergence4 == 0)
tube_length_guess_ft = tube_length_guess/12; % in
HX_EG_Re = (EG_density*velocity_btwn_tubes*tube_length_guess_ft)/EG_avg_viscosity; % dimensionless, Reynold's Number
HX_EG_Pr = (EG_avg_viscosity*EG_avg_spec_heat)/EG_avg_thermcond; % dimensionless, Prandtl Number
HX_EG_Nu = 0.0308*(HX_EG_Re^(4/5))*(HX_EG_Pr^(1/3)); % dimensionless, Nusselt Number HT relation for turbulent flow over flat plate w/ constant heat flux (Fund. Heat & Mass Flow, p.446 https://drive.google.com/file/d/1ygQNoKYAoQdCP6pjLTkWYne-Icz-8YX8/view?usp=sharing)
HX_EG_h_conv = (HX_EG_Nu*EG_avg_thermcond)/tube_length_guess_ft; % Btu/sft^2F, convective coefficient
HX_EG_thermR_conv = 1/HX_EG_h_conv; % sft^2F/Btu, convective thermal resistance

%% exhaust gas side radiative heat transfer

SS_emissivity = 0.585; % from https://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
SB_constant = 5.6703 * 10^(-8); % W/m^2K^4, Stefan-Boltzmann constant
HX_EG_avg_temp = (HX_EG_inlet_temp + HX_EG_outlet_temp)/2;
HX_EG_avg_temp_SI = (HX_EG_avg_temp-32)*(5/9) + 273.15 ; % F to K 
HX_EG_avg_outertubing_temp_initial_SI = (HX_EG_avg_outertubing_temp_initial-32)*(5/9) + 273.15 ; % F to K 
HX_EG_h_rad_SI = SS_emissivity*SB_constant*((HX_EG_avg_temp_SI^2)+(HX_EG_avg_outertubing_temp_initial_SI^2))...
    *(HX_EG_avg_temp_SI+HX_EG_avg_outertubing_temp_initial_SI); % W/m^2K, radiation coefficient from https://celsiainc.com/heat-sink-blog/fundamentals-of-thermal-resistance/
HX_EG_h_rad = HX_EG_h_rad_SI/5.678/3600; % Btu/sft^2F, radiation coefficient
HX_EG_thermR_rad = 1/HX_EG_h_rad; % sft^2F/Btu, radiative thermal resistance

%% pressurant side convective heat transfer
nominal_pressurant_mdot_pertube = nominal_pressurant_mdot/num_tubes;
inner_diameter_ft = inner_diameter/12;
outer_diameter_ft = outer_diameter/12;

% imperial to SI conversion
HX_pressurant_inlet_pressure_SI = HX_pressurant_inlet_pressure*6894.76; % psia to Pa
HX_pressurant_inlet_temp_SI = (HX_pressurant_inlet_temp-32)*(5/9) + 273.15; % F to K 
HX_pressurant_outlet_pressure_initial_SI = HX_pressurant_outlet_pressure_initial*6894.76; % psia to Pa
HX_pressurant_outlet_temp_SI = (HX_pressurant_outlet_temp-32)*(5/9) + 273.15; % F to K 
% inlet values
HX_pressurant_inlet_viscosity_SI = py.CoolProp.CoolProp.PropsSI('V','P',HX_pressurant_inlet_pressure_SI,'T',HX_pressurant_inlet_temp_SI,'Nitrogen'); % Pa*s
HX_pressurant_inlet_viscosity = HX_pressurant_inlet_viscosity_SI*0.671968994813; % Pa*s to lbm/ft*s
HX_pressurant_inlet_thermcond_SI = py.CoolProp.CoolProp.PropsSI('L','P',HX_pressurant_inlet_pressure_SI,'T',HX_pressurant_inlet_temp_SI,'Nitrogen'); % W/mK
HX_pressurant_inlet_thermcond = (HX_pressurant_inlet_thermcond_SI*0.5781759824)/3600; % W/mK to Btu/sftF
% HX_pressurant_inlet_Pr = py.CoolProp.CoolProp.PropsSI('Prandtl','P',HX_pressurant_inlet_pressure_SI,'T',HX_pressurant_inlet_temp_SI,'Nitrogen'); % dimensionless, Prandtl Number
HX_pressurant_inlet_spec_heat_SI = py.CoolProp.CoolProp.PropsSI('Cp0mass','P',HX_pressurant_inlet_pressure_SI,'T',HX_pressurant_inlet_temp_SI,'Nitrogen'); % J/kgK
HX_pressurant_inlet_spec_heat = HX_pressurant_inlet_spec_heat_SI*0.00023884589662749592; % J/kgK to Btu/lbm*F
HX_pressurant_inlet_Pr = (HX_pressurant_inlet_spec_heat_SI*HX_pressurant_inlet_viscosity_SI)/HX_pressurant_inlet_thermcond_SI; % dimensionless

% outlet values
HX_pressurant_outlet_viscosity_SI =py.CoolProp.CoolProp.PropsSI('V','P',HX_pressurant_outlet_pressure_initial_SI,'T',HX_pressurant_outlet_temp_SI,'Nitrogen'); % Pa*s
HX_pressurant_outlet_viscosity = HX_pressurant_outlet_viscosity_SI*0.671968994813; % Pa*s to lbm/ft*s
HX_pressurant_outlet_thermcond_SI = py.CoolProp.CoolProp.PropsSI('L','P',HX_pressurant_outlet_pressure_initial_SI,'T',HX_pressurant_outlet_temp_SI,'Nitrogen'); % W/mK
HX_pressurant_outlet_thermcond = (HX_pressurant_outlet_thermcond_SI*0.5781759824)/3600; % W/mK to Btu/sftF
% HX_pressurant_outlet_Pr = py.CoolProp.CoolProp.PropsSI('Prandtl','P',HX_pressurant_outlet_pressure_initial_SI,'T',HX_pressurant_outlet_temp_SI,'Nitrogen'); % dimensionless, Prandtl Number
HX_pressurant_outlet_spec_heat_SI = py.CoolProp.CoolProp.PropsSI('Cp0mass','P',HX_pressurant_outlet_pressure_initial_SI,'T',HX_pressurant_outlet_temp_SI,'Nitrogen'); % J/kgK
HX_pressurant_outlet_spec_heat = HX_pressurant_outlet_spec_heat_SI*0.00023884589662749592; % J/kgK to Btu/lbm*F
HX_pressurant_outlet_Pr = (HX_pressurant_outlet_spec_heat_SI*HX_pressurant_outlet_viscosity_SI)/HX_pressurant_outlet_thermcond_SI; % dimensionless

%average values 
HX_pressurant_avg_viscosity = (HX_pressurant_inlet_viscosity + HX_pressurant_outlet_viscosity)/2;
HX_pressurant_avg_thermcond = (HX_pressurant_inlet_thermcond + HX_pressurant_outlet_thermcond)/2;
HX_pressurant_avg_Pr = (HX_pressurant_inlet_Pr + HX_pressurant_outlet_Pr)/2;
HX_pressurant_avg_spec_heat = (HX_pressurant_inlet_spec_heat + HX_pressurant_outlet_spec_heat)/2;
% solve for thermal resistance
HX_pressurant_Re = (4*nominal_pressurant_mdot_pertube)/(pi*inner_diameter_ft*HX_pressurant_avg_viscosity); % dimensionless, Reynold's Number
HX_pressurant_Nu = 0.023*(HX_pressurant_Re^.8)*(HX_pressurant_avg_Pr^(1/3)); % dimensionless, Nusselt Number HT relation for internal flow through a tube from HX Dimensioning, Saari, pg. 63 ((https://drive.google.com/file/d/1DAyVj1EhW7SGOJYqiDseTSWOlxGQSsTN/view?usp=sharing)
HX_pressurant_h_conv = (HX_pressurant_Nu*HX_pressurant_avg_thermcond)/inner_diameter_ft; % Btu/sft^2F, convective coefficient
HX_pressurant_thermR_conv = 1/HX_pressurant_h_conv; % sft^2F/Btu, convective thermal resistance

%% pressurant side radiative heat transfer

HX_pressurant_avg_temp = (HX_pressurant_inlet_temp + HX_pressurant_outlet_temp)/2;
HX_pressurant_avg_temp_SI = (HX_pressurant_avg_temp-32)*(5/9) + 273.15 ; % F to K 
HX_pressurant_avg_innertubing_temp_initial_SI = (HX_pressurant_avg_innertubing_temp_initial-32)*(5/9) + 273.15; % F to K 
HX_pressurant_h_rad_SI = SS_emissivity*SB_constant*((HX_pressurant_avg_temp_SI^2)+(HX_pressurant_avg_innertubing_temp_initial_SI^2))...
    *(HX_pressurant_avg_temp_SI+HX_pressurant_avg_innertubing_temp_initial_SI); % W/m^2K, radiation coefficient from https://celsiainc.com/heat-sink-blog/fundamentals-of-thermal-resistance/
HX_pressurant_h_rad = HX_pressurant_h_rad_SI/5.678/3600; % Btu/sft^2F, radiation coefficient
HX_pressurant_thermR_rad = 1/HX_pressurant_h_rad; % sft^2F/Btu, radiative thermal resistance

%% overall heat transfer coefficient

overall_heat_transfer_coeff = 1/((((1/HX_pressurant_thermR_conv)+(1/HX_pressurant_thermR_rad))^-1)+(((1/HX_EG_thermR_conv)+(1/HX_EG_thermR_rad))^-1)); % btu/sft^2F, overall heat transfer coefficient
overall_heat_transfer_rate = min((nominal_pressurant_mdot_pertube*HX_pressurant_avg_spec_heat)...
    *(HX_pressurant_outlet_temp - HX_pressurant_inlet_temp),(GG_total_mdot*EG_avg_spec_heat)...
    *(HX_EG_inlet_temp - HX_EG_outlet_temp)); % BTU/s
surface_area_ft2 = overall_heat_transfer_rate/(overall_heat_transfer_coeff*LMTD); % ft^2
tube_length_ft = surface_area_ft2/(pi*((outer_diameter_ft+inner_diameter_ft)/2)); % ft
tube_length = tube_length_ft*12; % in

%% pressurant side tubing surface temperature

HX_pressurant_surface_area = pi*inner_diameter_ft*tube_length_ft; % ft^2
HX_pressurant_avg_innertubing_temp =(overall_heat_transfer_rate*(((1/HX_pressurant_thermR_conv+1/HX_pressurant_thermR_rad)...
    ^-1)/HX_pressurant_surface_area)+HX_pressurant_avg_temp); % F
if (abs(HX_pressurant_avg_innertubing_temp - HX_pressurant_avg_innertubing_temp_initial)/HX_pressurant_avg_innertubing_temp_initial)*100 < 1
    convergence2 = 1;
else
    convergence2 = 0;
    HX_pressurant_avg_innertubing_temp_initial = HX_pressurant_avg_innertubing_temp;
end

%% exhaust gas side tubing surface temperature

HX_EG_surface_area = pi*outer_diameter_ft*tube_length_ft; % ft^2
HX_EG_avg_outertubing_temp = HX_EG_avg_temp-overall_heat_transfer_rate*(((1/(HX_EG_thermR_conv/HX_EG_surface_area))+(1/(HX_EG_thermR_rad/HX_EG_surface_area)))^-1);
if (abs(HX_EG_avg_outertubing_temp - HX_EG_avg_outertubing_temp_initial)/HX_EG_avg_outertubing_temp_initial)*100 < 1
    convergence3 = 1;
else
    convergence3 = 0;
    HX_EG_avg_outertubing_temp_initial = HX_EG_avg_outertubing_temp;
end

%% convergence to final tube length

if (abs(tube_length - tube_length_guess)/tube_length_guess)*100 < 1
    convergence4 = 1;
else
    convergence4 = 0;
    tube_length_guess = tube_length;
end
end
end
