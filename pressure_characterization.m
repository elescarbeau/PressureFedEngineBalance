clear all
%% Pursuit Rocket Pressure Characterization
% Elysse Lescarbeau, 1/8/21
% main script to determine all major values related to fluids across system
% Note: current code does not include pressure drops due to tube bends or
% fittings (negligible dP influence)
% *Code currently takes ~1.5 minutes to run due to multiple convergences across main script and functions*

filename = 'Pressure_Characterization_Data.xlsx';
%% fluid properties

% ox and fuel fluid properties 
engine_ox_pressure = 600 + 14.7; % psia, required engine pressure into injector
engine_fuel_pressure = 750 + 14.7; % psia, required engine pressure into regen channels
engine_ox_mdot = 9.52; % lbm/s, required ox mass flow rate into engine
engine_fuel_mdot = 2.38; % lbm/s, required fuel mass flow rate into engine
ox_density = 64.5; % lbm/ft^3, assumed constant
fuel_density = 49.1; % lbm/ft^3, assumed constant 
ox_temp = -22; % F, assumed constant
fuel_temp = 70; % F, assumed constant
ox_viscosity_omega_term = (309.57-5.24)/(243.15-5.24); % from http://edge.rit.edu/edge/P07106/public/Nox.pdf
ox_viscosity_SI =  0.0000293423*exp(1.6089*((ox_viscosity_omega_term-1)^(1/3))+2.0439*((ox_viscosity_omega_term-1)^(4/3))); % N s/m2, from http://edge.rit.edu/edge/P07106/public/Nox.pdf
ox_viscosity = ox_viscosity_SI*0.02088543423315013; % lbf*s/ft^2, assumed constant, from https://unitconverterapp.com/quantity/dynamic%20viscosity
fuel_viscosity = 2.4*2.088543423315013e-05; % mPa*s to lbf*s/ft^2, assumed constant, from https://www.shell.com/business-customers/chemicals/our-products/solvents-chemical/alcohols/_jcr_content/par/tabbedcontent/tab/textimage.stream/1460023044637/b0b7ef356c65f2fdd262453bd3eb009011c1bfd0/ipa-s1111-eu.pdf and https://unitconverterapp.com/quantity/dynamic%20viscosity 

% pressurant (nitrogen) fluid properties
nitrogen_spec_heat_ratio = 1.4; % dimensionless, assumed constant
nitrogen_gas_constant = 55.17; % ftlbf/lbmR, assumed constant
pressurant_temp = 200; % F, temperature of pressurant in propellant tanks 
grav_constant = 32.17405; % ft/s^2, gravitational constant

% tank properties
empty_ox_dome_volume = 144.99; % in^3, from SW
empty_ox_tank_volume = 16827.42; % in^3, from SW
empty_ox_commondome_volume = 127.69; % in^3, from SW
empty_fuel_commondome_volume = 143.29; % in^3, from SW
empty_fuel_tank_volume = 7084.4; % in^3, from SW
empty_fuel_dome_volume = 149.72; % in^3, from SW
ullage_percent = 5;
fuel_residual_safety_factor = 3; % percent, ensures engine does not run ox rich at end of burn
tank_radius = 5.795; % in

% tubing properties
SS_tube_roughness = 0.00059; % in, from https://www.enggcyclopedia.com/2011/09/absolute-roughness/

% COPV properties
COPVs_pressure = 4500 + 14.7; % psia, pressurant pressure in COPVs
COPVs_internal_volume = 0.318; % ft^3, COPVS internal volume, from http://airtanksforsale.com/
COPVs_external_volume = 0.454; % ft^3, COPVs external volume, measured

% HX properties 
HX_outer_diameter = .375; % in
HX_inner_diameter = .305; % in
HX_num_tubes = 8;
HX_pressurant_outlet_temp = 400; % F
HX_EG_inlet_temp = 1300; % F
HX_EG_outlet_temp = 500; % F
EG_gas_constant_SI = 601.237; % J/kgK, exhaust gas value from NASA CEA 
EG_spec_heat_ratio = 1.1321; % dimensionless, exhaust gas value from NASA CEA
EG_avg_spec_heat = 2.611; %Btu/lbmF, from NASA CEA, average between chamber and "throat"
EG_density = 0.0933; % lbm/ft^3, from NASA CEA, average between chamber and "throat"
EG_avg_viscosity = 1.46*(10^-5); % lbm/ft*s, used weighted average of individual gas output products & NIST (see "HEX HT Coefficient" tab here:https://docs.google.com/spreadsheets/d/1Q-QfuCvjTrlBRowdYNJNm_TELQPWODFdH_oxyiCvEUc/edit?usp=sharing) 
EG_avg_thermcond = 7.47*(10^-5); % Btu/sftF, used weighted average of individual gas output products & NIST (see "HEX HT Coefficient" tab here:https://docs.google.com/spreadsheets/d/1Q-QfuCvjTrlBRowdYNJNm_TELQPWODFdH_oxyiCvEUc/edit?usp=sharing) 

%% read in data needed for calculations

% post tank data %
Ox_PostTank_Component_Info_Table = readtable(filename,'Sheet',1);
Fuel_PostTank_Component_Info_Table = readtable(filename,'Sheet',2);
Ox_PostTank_Line_Info_Table = readtable(filename,'Sheet',3);
Fuel_PostTank_Line_Info_Table = readtable(filename,'Sheet',4);
size_mat_ox_posttank_comp = size(Ox_PostTank_Component_Info_Table);
size_mat_fuel_posttank_comp = size(Fuel_PostTank_Component_Info_Table);
size_mat_ox_posttank_line = size(Ox_PostTank_Line_Info_Table);
size_mat_fuel_posttank_line = size(Fuel_PostTank_Line_Info_Table);
Ox_PostTank_Component_Info_Matrix = table2array(Ox_PostTank_Component_Info_Table(1:size_mat_ox_posttank_comp(1),2:size_mat_ox_posttank_comp(2)));
Fuel_PostTank_Component_Info_Matrix = table2array(Fuel_PostTank_Component_Info_Table(1:size_mat_fuel_posttank_comp(1),2:size_mat_fuel_posttank_comp(2)));
Ox_PostTank_Line_Info_Matrix = table2array(Ox_PostTank_Line_Info_Table(1:size_mat_ox_posttank_line(1),2:size_mat_ox_posttank_line(2)));
Fuel_PostTank_Line_Info_Matrix = table2array(Fuel_PostTank_Line_Info_Table(1:size_mat_fuel_posttank_line(1),2:size_mat_fuel_posttank_line(2)));

% pre tank data
PreTank_Component_Info_Table = readtable('Pressure_Characterization_Data.xlsx','Sheet',5);
PreTank_Line_Info_Table = readtable('Pressure_Characterization_Data.xlsx','Sheet',6);
Ox_PreTank_Component_Info_Table = readtable('Pressure_Characterization_Data.xlsx','Sheet',7);
Fuel_PreTank_Component_Info_Table = readtable('Pressure_Characterization_Data.xlsx','Sheet',8);
Ox_PreTank_Line_Info_Table = readtable('Pressure_Characterization_Data.xlsx','Sheet',9);
Fuel_PreTank_Line_Info_Table = readtable('Pressure_Characterization_Data.xlsx','Sheet',10);
size_mat_pretank_comp = size(PreTank_Component_Info_Table);
size_mat_pretank_line = size(PreTank_Line_Info_Table);
size_mat_ox_pretank_comp = size(Ox_PreTank_Component_Info_Table);
size_mat_fuel_pretank_comp = size(Fuel_PreTank_Component_Info_Table);
size_mat_ox_pretank_line = size(Ox_PreTank_Line_Info_Table);
size_mat_fuel_pretank_line = size(Fuel_PreTank_Line_Info_Table);
PreTank_Component_Info_Matrix = table2array(PreTank_Component_Info_Table(1:size_mat_pretank_comp(1),2:size_mat_pretank_comp(2)));
PreTank_Line_Info_Matrix = table2array(PreTank_Line_Info_Table(1:size_mat_pretank_line(1),2:size_mat_pretank_line(2)));
Ox_PreTank_Component_Info_Matrix = table2array(Ox_PreTank_Component_Info_Table(1:size_mat_ox_pretank_comp(1),2:size_mat_ox_pretank_comp(2)));
Fuel_PreTank_Component_Info_Matrix = table2array(Fuel_PreTank_Component_Info_Table(1:size_mat_fuel_pretank_comp(1),2:size_mat_fuel_pretank_comp(2)));
Ox_PreTank_Line_Info_Matrix = table2array(Ox_PreTank_Line_Info_Table(1:size_mat_ox_pretank_line(1),2:size_mat_ox_pretank_line(2)));
Fuel_PreTank_Line_Info_Matrix = table2array(Fuel_PreTank_Line_Info_Table(1:size_mat_fuel_pretank_line(1),2:size_mat_fuel_pretank_line(2)));

% final empty ox and fuel volumes
empty_ox_volume_total = empty_ox_dome_volume + empty_ox_tank_volume + empty_ox_commondome_volume + (pi/4*(Ox_PostTank_Line_Info_Matrix(1,2)^2))*Ox_PostTank_Line_Info_Matrix(1,1); % in^3
empty_fuel_volume_total = empty_fuel_commondome_volume + empty_fuel_tank_volume + empty_fuel_dome_volume; % in^3

%% find pressurant mass flow rate

% determine tank pressure & nitrogen mass flow rate via required pressure
% into engine and losses up to tank %

% max hydrostatic pressure %
max_ox_hydrostatic_pressure_guess = 20.32; % psid
max_fuel_hydrostatic_pressure_guess = 4.75; % psid
convergence_hydrostaticpress = 0;
while convergence_hydrostaticpress == 0

% post-tank ox side pressure losses %
ox_posttank_component_losses = zeros(size_mat_ox_posttank_comp(1),1);
ox_posttank_line_losses = zeros(size_mat_ox_posttank_line(1),1);

for j = 1:size_mat_ox_posttank_comp(1)
    [~,~,ox_posttank_component_losses(j,1)] = componentdP(engine_ox_mdot,Ox_PostTank_Component_Info_Matrix(j,1),Ox_PostTank_Component_Info_Matrix(j,2),ox_density); % psid
end

for j= 1:size_mat_ox_posttank_line(1)
    [ox_posttank_line_losses(j,1)] = majorlineloss(engine_ox_mdot,Ox_PostTank_Line_Info_Matrix(j,1),Ox_PostTank_Line_Info_Matrix(j,2),SS_tube_roughness,ox_density,ox_viscosity); % psid
end

total_ox_posttank_losses=sum(ox_posttank_component_losses)+sum(ox_posttank_line_losses); % psid

% post-tank fuel side pressure losses
fuel_posttank_component_losses = zeros(size_mat_fuel_posttank_comp(1),1);
fuel_posttank_line_losses = zeros(size_mat_fuel_posttank_line(1),1);

for j= 1:size_mat_fuel_posttank_comp(1)
    [~,~,fuel_posttank_component_losses(j,1)] = componentdP(engine_fuel_mdot,Fuel_PostTank_Component_Info_Matrix(j,1),Fuel_PostTank_Component_Info_Matrix(j,2),fuel_density); % psid
end

for j= 1:size_mat_fuel_posttank_line(1)
    [fuel_posttank_line_losses(j,1)] = majorlineloss(engine_fuel_mdot,Fuel_PostTank_Line_Info_Matrix(j,1),Fuel_PostTank_Line_Info_Matrix(j,2),SS_tube_roughness,fuel_density,fuel_viscosity); % psid
end

total_fuel_posttank_losses=sum(fuel_posttank_component_losses)+sum(fuel_posttank_line_losses); % psid

% initial COPV temperature guess 
COPVs_pressurant_temp_guess = -6.67 ; % F, initial guess at COPVs outlet temp
convergence_COPVsoutletemp = 0; 
while convergence_COPVsoutletemp == 0

COPVs_pressurant_density = (COPVs_pressure*144)/(nitrogen_gas_constant*(COPVs_pressurant_temp_guess+459.67));

% initial GG values
GG_total_mdot_guess = .08; % lbm/s, required mass flow rate into GG initial guess
GG_OF_ratio = .4;
convergence_GGmdot = 0;
while convergence_GGmdot == 0
GG_ox_mdot_guess = (GG_OF_ratio*GG_total_mdot_guess)/(GG_OF_ratio+1); % lbm/s, required ox mass flow rate into GG
GG_fuel_mdot_guess = GG_total_mdot_guess/(GG_OF_ratio+1); % lbm/s, required fuel mass flow rate into GG

% determine ox side pressurant mass flow rate
ox_required_tank_pressure = engine_ox_pressure + total_ox_posttank_losses - max_ox_hydrostatic_pressure_guess; % psia
ox_vol_flow_rate_total = (engine_ox_mdot+GG_ox_mdot_guess)/ox_density; % ft^3/s
ox_pressurant_density_tank = (ox_required_tank_pressure*144)/(nitrogen_gas_constant*(pressurant_temp+459.67)); % lbm/ft^3
ox_pressurant_mdot = ox_vol_flow_rate_total*ox_pressurant_density_tank; % lbm/s

% determine fuel side pressurant mass flow rate
fuel_required_tank_pressure = engine_fuel_pressure + total_fuel_posttank_losses - max_fuel_hydrostatic_pressure_guess; % psia
fuel_vol_flow_rate_total = (engine_fuel_mdot+GG_fuel_mdot_guess)/fuel_density; % ft^3/s
fuel_pressurant_density_tank = (fuel_required_tank_pressure*144)/(nitrogen_gas_constant*(pressurant_temp+459.67)); % lbm/ft^3
fuel_pressurant_mdot = fuel_vol_flow_rate_total*fuel_pressurant_density_tank; % lbm/s

% determine total pressurant mass flow rate
total_pressurant_mdot = ox_pressurant_mdot+fuel_pressurant_mdot; % lbm/s

%% determine pre-tank pressure losses, HX sizing, reg set pressure 

% preallocate matrices %
% for line losses, matrices output columns are:
% 1) P_2 (psia)
% 2)rho_2 (lbm/ft^3)
% 3)dP (psid)
% 4)T_2 (F)
% 5) v_2 (ft/s)
% 6) M_2 (-)
% 7) dT (F)
% for component losses, matrices output columns are: 
% 1) P_2 (psia)
% 2) rho_2 (lbm/ft^3)
% 3) dP (psid) 
% (additional column values in matrices that are unused will remain 0)

pretank_losses=zeros(size_mat_pretank_comp(1)+size_mat_pretank_line(1)+1,7); % adds additional row for HX loss
size_mat_pretank_losses = size(pretank_losses);
ox_pretank_losses=zeros(size_mat_ox_pretank_comp(1)+size_mat_ox_pretank_line(1),7);
size_mat_ox_pretank_losses = size(ox_pretank_losses);
fuel_pretank_losses=zeros(size_mat_fuel_pretank_comp(1)+size_mat_fuel_pretank_line(1),7);
size_mat_fuel_pretank_losses = size(fuel_pretank_losses);

% pre-regulator pressure drop due to line loss %
[pretank_losses(1,1),pretank_losses(1,2),pretank_losses(1,3),pretank_losses(1,4),pretank_losses(1,5),pretank_losses(1,6),pretank_losses(1,7)] = fannocalc(total_pressurant_mdot,COPVs_pressure,COPVs_pressurant_temp_guess,COPVs_pressurant_density,PreTank_Line_Info_Matrix(1,1),PreTank_Line_Info_Matrix(1,2));

% post-regulator pressure drops %
% 1) iterate through reg set pressures and calculate all pressure drops due to components and line lengths (including HX) until reg set pressure less pressure drops is
% equal to required tank pressure for either ox or fuel tanks 
% 2) calculate trim orifice size needed to reach required tank pressure for whichever tank pressure is not met by component and line pressure drops
reg_set_pressure_guess = engine_fuel_pressure*1.5; % ensures initial guess is not smaller than actual required pressure so functions work
reg_outlet_temp_guess = (((reg_set_pressure_guess/pretank_losses(1,1))^(1-(1/nitrogen_spec_heat_ratio)))*(pretank_losses(1,4)+459.67))-459.67; % F, temperature drop to due isentropic expansion across regulator
reg_outlet_density_guess = (reg_set_pressure_guess*144)/(nitrogen_gas_constant*(reg_outlet_temp_guess+459.67)); % lbm/ft^3
convergence_reg = 0;
while convergence_reg == 0
    % line loss from regulator to first component %
    [pretank_losses(2,1),pretank_losses(2,2),pretank_losses(2,3),pretank_losses(2,4),pretank_losses(2,5),pretank_losses(2,6),pretank_losses(2,7)]... 
        =fannocalc(total_pressurant_mdot,reg_set_pressure_guess,reg_outlet_temp_guess,reg_outlet_density_guess,PreTank_Line_Info_Matrix(2,1),PreTank_Line_Info_Matrix(2,2));

    % loop over components and then line lengths to put losses in correct order
    % losses up to HX %
    count = 3; % count to index into pretank loss matrix
    for i = 1:size_mat_pretank_comp(1)
        [pretank_losses(count,1),pretank_losses(count,2),pretank_losses(count,3)]...
            = componentdP(total_pressurant_mdot,PreTank_Component_Info_Matrix(i,1),PreTank_Component_Info_Matrix(i,2),pretank_losses(count-1,2),pretank_losses(count-1,4),pretank_losses(count-1,1));
        count = count + 1;
        [pretank_losses(count,1),pretank_losses(count,2),pretank_losses(count,3),pretank_losses(count,4),pretank_losses(count,5),pretank_losses(count,6),pretank_losses(count,7)]...
            = fannocalc(total_pressurant_mdot,pretank_losses(count-1,1),pretank_losses(count-2,4),pretank_losses(count-1,2),PreTank_Line_Info_Matrix(i+2,1),PreTank_Line_Info_Matrix(i+2,2));
        HXcount = count; % to be used for indexing into pretank loss matrices for HX dP function to pull pre-HX values
        count = count + 1;
    end
    
    count = HXcount + 1;
    
    % HX sizing %
    % convergence for outlet pressure 
    HX_pressurant_outlet_pressure_initial = pretank_losses(HXcount,1) - .01; 
    convergence_HX = 0;
    while convergence_HX == 0 
    [nominal_pressurant_mdot,overall_heat_transfer_coeff,surface_area_ft2,HX_tube_length,space_btwn_tubes,HX_pressurant_avg_viscosity,overall_heat_transfer_rate,HX_pressurant_thermR_conv,HX_pressurant_thermR_rad,HX_EG_thermR_conv,HX_EG_thermR_rad,HX_EG_avg_outertubing_temp,HX_pressurant_avg_innertubing_temp,tube_length_guess,HX_pressurant_avg_Pr] = HXsizing(GG_total_mdot_guess,EG_density,...
    EG_avg_viscosity,EG_avg_thermcond,EG_avg_spec_heat,total_pressurant_mdot,HX_num_tubes,HX_inner_diameter,HX_outer_diameter,...
    pretank_losses(HXcount,1),HX_pressurant_outlet_pressure_initial,...
    pretank_losses(HXcount,4),HX_pressurant_outlet_temp,pressurant_temp,...
    HX_EG_inlet_temp,HX_EG_outlet_temp);

    % HX loss 
    [pretank_losses(count,1),pretank_losses(count,2),pretank_losses(count,3),friction_coeff] = HXpressuredrop(nominal_pressurant_mdot,...
        HX_num_tubes,HX_tube_length,HX_inner_diameter,HX_outer_diameter,SS_tube_roughness,pretank_losses(HXcount,1),pretank_losses(HXcount,2),pressurant_temp,HX_pressurant_avg_viscosity);
    if (abs(pretank_losses(count,1) - HX_pressurant_outlet_pressure_initial)/HX_pressurant_outlet_pressure_initial)*100 < 1 
        convergence_HX = 1;
    else
        HX_pressurant_outlet_pressure_initial = pretank_losses(count,1);
    end
    end
%     pretank_losses(count,1) = HX_pressurant_outlet_pressure_initial;
    count = count + 1;
    
    % post HX loss up to ox/fuel split %
    [pretank_losses(count,1),pretank_losses(count,2),pretank_losses(count,3),pretank_losses(count,4),pretank_losses(count,5),pretank_losses(count,6),pretank_losses(count,7)]...
            = fannocalc(total_pressurant_mdot,pretank_losses(count-1,1),pressurant_temp,pretank_losses(count-1,2),PreTank_Line_Info_Matrix(size_mat_pretank_line(1),1),PreTank_Line_Info_Matrix(size_mat_pretank_line(1),2));
    count = count + 1;
    
    % ox side losses %
    % initial ox loss
     [ox_pretank_losses(1,1),ox_pretank_losses(1,2),ox_pretank_losses(1,3),ox_pretank_losses(1,4),ox_pretank_losses(1,5),ox_pretank_losses(1,6),ox_pretank_losses(1,7)]...
            = fannocalc(ox_pressurant_mdot,pretank_losses(count-1,1),pressurant_temp,pretank_losses(count-1,2),Ox_PreTank_Line_Info_Matrix(1,1),Ox_PreTank_Line_Info_Matrix(1,2));  
    
    count_OF = 2; % count to index into ox & fuel losses matrix
    for i = 1:size_mat_ox_pretank_comp(1)
        [ox_pretank_losses(count_OF,1),ox_pretank_losses(count_OF,2),ox_pretank_losses(count_OF,3)]...
            = componentdP(ox_pressurant_mdot,Ox_PreTank_Component_Info_Matrix(i,1),Ox_PreTank_Component_Info_Matrix(i,2),ox_pretank_losses(count_OF-1,2),ox_pretank_losses(count_OF-1,4),ox_pretank_losses(count_OF-1,1));
        count_OF = count_OF + 1;
        [ox_pretank_losses(count_OF,1),ox_pretank_losses(count_OF,2),ox_pretank_losses(count_OF,3),ox_pretank_losses(count_OF,4),ox_pretank_losses(count_OF,5),ox_pretank_losses(count_OF,6),ox_pretank_losses(count_OF,7)]...
            = fannocalc(ox_pressurant_mdot,ox_pretank_losses(count_OF-1,1),ox_pretank_losses(count_OF-2,4),ox_pretank_losses(count_OF-1,2),Ox_PreTank_Line_Info_Matrix(i,1),Ox_PreTank_Line_Info_Matrix(i,2));
        count_OF = count_OF + 1;
    end
    ox_tank_inlet_pressure = ox_pretank_losses(count_OF-1,1);
        
    % fuel side losses %
    % initial fuel loss
     [fuel_pretank_losses(1,1),fuel_pretank_losses(1,2),fuel_pretank_losses(1,3),fuel_pretank_losses(1,4),fuel_pretank_losses(1,5),fuel_pretank_losses(1,6),fuel_pretank_losses(1,7)]...
            = fannocalc(fuel_pressurant_mdot,pretank_losses(count-1,1),pressurant_temp,pretank_losses(count-1,2),Fuel_PreTank_Line_Info_Matrix(1,1),Fuel_PreTank_Line_Info_Matrix(1,2));  
    
    count_OF = 2; % count to index into ox & fuel losses matrix
    for i = 1:size_mat_fuel_pretank_comp(1)
        [fuel_pretank_losses(count_OF,1),fuel_pretank_losses(count_OF,2),fuel_pretank_losses(count_OF,3)]...
            = componentdP(fuel_pressurant_mdot,Fuel_PreTank_Component_Info_Matrix(i,1),Fuel_PreTank_Component_Info_Matrix(i,2),fuel_pretank_losses(count_OF-1,2),fuel_pretank_losses(count_OF-1,4),fuel_pretank_losses(count_OF-1,1));
        count_OF = count_OF + 1;
        [fuel_pretank_losses(count_OF,1),fuel_pretank_losses(count_OF,2),fuel_pretank_losses(count_OF,3),fuel_pretank_losses(count_OF,4),fuel_pretank_losses(count_OF,5),fuel_pretank_losses(count_OF,6),fuel_pretank_losses(count_OF,7)]...
            = fannocalc(fuel_pressurant_mdot,fuel_pretank_losses(count_OF-1,1),fuel_pretank_losses(count_OF-2,4),fuel_pretank_losses(count_OF-1,2),Fuel_PreTank_Line_Info_Matrix(i,1),Fuel_PreTank_Line_Info_Matrix(i,2));
        count_OF = count_OF + 1;
    end
    fuel_tank_inlet_pressure = fuel_pretank_losses(count_OF-1,1);
   
    pressure_margin = 1; % psid
    ox_pressure_dP = ox_tank_inlet_pressure-ox_required_tank_pressure;
    fuel_pressure_dP = fuel_tank_inlet_pressure-fuel_required_tank_pressure;
    
    if (ox_pressure_dP < pressure_margin && ox_pressure_dP > 0 && fuel_pressure_dP > 0) || (fuel_pressure_dP < pressure_margin && fuel_pressure_dP > 0 && ox_pressure_dP > 0)
        convergence_reg = 1;
    else
        reg_set_pressure_guess = reg_set_pressure_guess - 0.5;
        reg_outlet_temp_guess = (((reg_set_pressure_guess/pretank_losses(1,1))^(1-(1/nitrogen_spec_heat_ratio)))*(pretank_losses(1,4)+459.67))-459.67; % F, temperature drop to due isentropic expansion across regulator
        reg_outlet_density_guess = (reg_set_pressure_guess*144)/(nitrogen_gas_constant*(reg_outlet_temp_guess+459.67)); % lbm/ft^3
    end
end

reg_set_pressure = reg_set_pressure_guess; % psia
reg_set_pressure_psig = reg_set_pressure-14.7; % psig
reg_outlet_temp = reg_outlet_temp_guess; % F
reg_outlet_density = reg_outlet_density_guess; % lbm/ft^3

% calculate orifice size %
if (ox_pressure_dP < pressure_margin && ox_pressure_dP > 0 && fuel_pressure_dP > 0)
    orifice_side = {'Fuel'};
    [orifice_diameter,orifice_dP] = orifice_calc(fuel_pressurant_mdot,Fuel_PreTank_Line_Info_Matrix(size_mat_fuel_pretank_line(1),2),fuel_tank_inlet_pressure,fuel_required_tank_pressure,fuel_pretank_losses(count_OF-1,2),3);
else
        orifice_side = {'Ox'};
        [orifice_diameter,orifice_dP] = orifice_calc(ox_pressurant_mdot,Ox_PreTank_Line_Info_Matrix(size_mat_ox_pretank_line(1),2),ox_tank_inlet_pressure,ox_required_tank_pressure,ox_pretank_losses(count_OF-1,2),3);
end

%% convergence to find final GG mdot

[~,~,GG_total_mdot] = GGmassflowrate(GG_OF_ratio,engine_ox_mdot,engine_fuel_mdot,...
    ox_density,fuel_density,ox_required_tank_pressure,fuel_required_tank_pressure,pretank_losses(HXcount,1),...
    pretank_losses(HXcount+1,1),pretank_losses(HXcount,4),HX_pressurant_outlet_temp,...
    pressurant_temp,HX_EG_inlet_temp,HX_EG_outlet_temp,EG_gas_constant_SI,EG_spec_heat_ratio); 
% GG mdot will automatically update with the below and be fed into the HX
% sizing function which will then feed into the HX pressure drop function
if (abs(GG_total_mdot - GG_total_mdot_guess)/GG_total_mdot_guess)*100 < 1
    convergence_GGmdot = 1;
else
    GG_total_mdot_guess = GG_total_mdot;
end
end

% final GG mdot values %
GG_ox_mdot = GG_ox_mdot_guess; % lbm/s
GG_fuel_mdot = GG_fuel_mdot_guess; % lbm/s
GG_total_mdot = GG_total_mdot_guess; % lbm/s
    
%% determine minimum blowdown pressure to determine mass of pressurant needed

% find minimum required pre-regulator blowdown pressure during engine burn phase 
% of flight to meet pressurant flow rate requirements (based on min dP across flow area
% of regulator) %

min_blowdown_pressure_guess = pretank_losses(1,1); % psia
reg_flow_area = 1; % regulator Cv
reg_flow_area_CdA = pi/4*0.65*((.227*sqrt(reg_flow_area))^2); % regulator Cv to CdA conversion
convergence_blowdown = 0;
initial_temp = pretank_losses(1,4);
initial_temp_SI = (pretank_losses(1,4)-32)*5/9 + 273.15; % F to K

while convergence_blowdown == 0
    density_blowdown_SI = py.CoolProp.CoolProp.PropsSI('D','P',min_blowdown_pressure_guess*6894.76,'T',initial_temp_SI,'Nitrogen'); % kg/m^3
    density_blowdown = density_blowdown_SI*0.062428; % kg/m^3 to lbm/ft^3
    mdot_calculated_blowdown = (reg_flow_area_CdA/144)*density_blowdown...
    *((reg_set_pressure/min_blowdown_pressure_guess)^(1/nitrogen_spec_heat_ratio))*sqrt(2*grav_constant*nitrogen_gas_constant*(initial_temp+459.67)...
    *(nitrogen_spec_heat_ratio/(nitrogen_spec_heat_ratio-1))*(1-(reg_set_pressure/min_blowdown_pressure_guess)...
    ^((nitrogen_spec_heat_ratio-1)/nitrogen_spec_heat_ratio))); 
    if (abs(mdot_calculated_blowdown - total_pressurant_mdot)/total_pressurant_mdot)*100 < 1
        convergence_blowdown = 1 ;
    else
        min_blowdown_pressure_guess = min_blowdown_pressure_guess - .5;
    end
end

min_blowdown_pressure = min_blowdown_pressure_guess; % psia

% find required pressurant mass for engine using required tank pressures and minimum
% blowdown pressure

%% find total prop volume, amount of pressurant needed, # of COPVs, burn time

ox_prop_volume_guess = 7.25; % ft^3, total volume needed including engine & GG
fuel_prop_volume_guess = 2.63; % ft^3, total volume needed including engine & GG
convergence_vol = 0;
while convergence_vol == 0 
    % total pressurant needed
    ox_pressurant_mass = (ox_prop_volume_guess*144*ox_required_tank_pressure*nitrogen_spec_heat_ratio)...
        /((nitrogen_gas_constant*(pressurant_temp+459.67)*(1-min_blowdown_pressure/COPVs_pressure))); % lbm
    fuel_pressurant_mass = (fuel_prop_volume_guess*144*fuel_required_tank_pressure*nitrogen_spec_heat_ratio)...
        /((nitrogen_gas_constant*(pressurant_temp+459.67)*(1-min_blowdown_pressure/COPVs_pressure))); % lbm
    total_pressurant_mass = ox_pressurant_mass + fuel_pressurant_mass; % lbm
    
    % # of COPVs needed
    initial_ox_num_COPVs = 5; % initial max amount of COPVs in ox tank before storing remaining mass in fuel tank
    max_fuel_num_COPVs = 2; % max amount of COPVs that can fit in fuel tank
    num_COPVs_ox_tank_only =ceil(((total_pressurant_mass*nitrogen_gas_constant*(ox_temp+459.67))/(COPVs_pressure*144))/COPVs_internal_volume); % number of COPVs needed if only stored in ox tank
    if num_COPVs_ox_tank_only <= initial_ox_num_COPVs
        remaining_mass_after_ox = 0; % lbm
        total_num_COPVs_ox_tank = num_COPVs_ox_tank_only;
        total_num_COPVs_fuel_tank = 0;
    else
        remaining_mass_after_ox = total_pressurant_mass-((COPVs_pressure*144)*(COPVs_internal_volume*initial_ox_num_COPVs))/(nitrogen_gas_constant*(ox_temp+459.67)); % lbm
        num_COPVs_fuel_tank = ceil(((remaining_mass_after_ox*nitrogen_gas_constant*(fuel_temp+459.67))/(COPVs_pressure*144))/COPVs_internal_volume);
        if num_COPVs_fuel_tank <= max_fuel_num_COPVs
            remaining_mass_after_fuel = 0; % lbm
            total_num_COPVs_fuel_tank = num_COPVs_fuel_tank;
            total_num_COPVs_ox_tank = initial_ox_num_COPVs;
        else
            remaining_mass_after_fuel = remaining_mass_after_ox-((COPVs_pressure*144)*(COPVs_internal_volume*max_fuel_num_COPVs))/(nitrogen_gas_constant*(fuel_temp+459.67)); % lbm
            num_COPVs_ox_tank_additional = ceil(((remaining_mass_after_fuel*nitrogen_gas_constant*(ox_temp+459.67))/(COPVs_pressure*144))/COPVs_internal_volume);
            total_num_COPVs_ox_tank = initial_ox_num_COPVs + num_COPVs_ox_tank_additional;
            total_num_COPVs_fuel_tank = max_fuel_num_COPVs;
        end
    end
    total_num_COPVs = total_num_COPVs_ox_tank + total_num_COPVs_fuel_tank;
    COPVs_pressurant_temp = (((ox_temp+459.67)*total_num_COPVs_ox_tank+(fuel_temp+459.67)*total_num_COPVs_fuel_tank)/total_num_COPVs)-459.67; % F, weighted average of temperature of COPVs in ox and fuel tanks after reaching equilibrium with propellant temperatures
    
    % find approximate available volume in each tank and burn time
    vol_safety_factor = 5; % to account for additional hardware, tubing, etc. in tanks
    available_ox_volume = (empty_ox_volume_total/(12^3))-((total_num_COPVs_ox_tank*COPVs_external_volume)*(1+vol_safety_factor/100)); % ft^3
    available_fuel_volume = (empty_fuel_volume_total/(12^3))-((total_num_COPVs_fuel_tank*COPVs_external_volume)*(1+vol_safety_factor/100)); % ft^3
    max_possible_ox_mass = (available_ox_volume* (1-ullage_percent/100))*ox_density; % lbm, includes space for ullage
    max_possible_fuel_mass = (available_fuel_volume* (1-ullage_percent/100))*fuel_density; % lbm, includes space for ullage 
    total_ox_mdot = engine_ox_mdot + GG_ox_mdot; % lbm/s
    total_fuel_mdot = engine_fuel_mdot + GG_fuel_mdot; % lbm/s
    burn_time = 10; % seconds, initial guess
    convergence_burntime = 0;
    while convergence_burntime == 0
        total_ox_mass = total_ox_mdot*burn_time; %lbm
        total_fuel_mass = (engine_fuel_mdot*burn_time)*(1+fuel_residual_safety_factor/100)+ GG_fuel_mdot*burn_time; %lbm
        if ((abs(total_ox_mass - max_possible_ox_mass)/max_possible_ox_mass)*100 < 1) && (total_ox_mass <= max_possible_ox_mass) && (total_fuel_mass <= max_possible_fuel_mass)
                convergence_burntime = 1;
            elseif ((abs(total_fuel_mass - max_possible_fuel_mass)/max_possible_fuel_mass)*100 < 1) && (total_ox_mass <= max_possible_ox_mass) && (total_fuel_mass <= max_possible_fuel_mass)
                convergence_burntime = 1;
            else
                burn_time = burn_time + .2;
        end
    end
    ox_prop_volume = total_ox_mass/ox_density; % ft^3, final ox volume
    fuel_prop_volume = total_fuel_mass/fuel_density; % ft^3, final fuel volume
    if (abs(ox_prop_volume - ox_prop_volume_guess)/ox_prop_volume_guess)*100 < 1 && (abs(fuel_prop_volume - fuel_prop_volume_guess)/fuel_prop_volume_guess)*100 < 1
        convergence_vol = 1;
    else
        ox_prop_volume_guess = ox_prop_volume;
        fuel_prop_volume_guess = fuel_prop_volume;
    end
end
ox_prop_volume = ox_prop_volume_guess; % ft^3
fuel_prop_volume = fuel_prop_volume_guess; % ft^3

    % convergence for COPVs outlet temp
    if (abs(COPVs_pressurant_temp - COPVs_pressurant_temp_guess)/COPVs_pressurant_temp_guess)*100 < 1
        convergence_COPVsoutletemp = 1;
    else 
        COPVs_pressurant_temp_guess = COPVs_pressurant_temp;
    end
end
COPVs_pressurant_temp = COPVs_pressurant_temp_guess;

initial_acceleration = 39.7; % ft/s^2
ox_height = (((ox_prop_volume +(total_num_COPVs_ox_tank*COPVs_external_volume))...
    *(1+vol_safety_factor/100))/(pi*((tank_radius/12)^2)))+(Ox_PostTank_Line_Info_Matrix(1,1)/12); % ft, height of ox fluid column
fuel_height = (((fuel_prop_volume +(total_num_COPVs_fuel_tank*COPVs_external_volume))...
    *(1+vol_safety_factor/100))/(pi*((tank_radius/12)^2)))+(Fuel_PostTank_Line_Info_Matrix(1,1)/12); % ft, height of fuel fluid column
max_ox_hydrostatic_pressure = ((ox_density/grav_constant)*(grav_constant+initial_acceleration)*ox_height)/(12^2); % psid
max_fuel_hydrostatic_pressure = ((fuel_density/grav_constant)*(grav_constant+initial_acceleration)*fuel_height)/(12^2); % psid

if ((abs(max_ox_hydrostatic_pressure - max_ox_hydrostatic_pressure_guess)/max_ox_hydrostatic_pressure_guess)*100 < 1)...
        && ((abs(max_fuel_hydrostatic_pressure - max_fuel_hydrostatic_pressure_guess)/max_fuel_hydrostatic_pressure_guess)*100 < 1)
    convergence_hydrostaticpress = 1;
else
    max_ox_hydrostatic_pressure_guess = max_ox_hydrostatic_pressure;
    max_fuel_hydrostatic_pressure_guess = max_fuel_hydrostatic_pressure;
end
end

max_ox_hydrostatic_pressure = max_ox_hydrostatic_pressure_guess; % psid
max_fuel_hydrostatic_pressure = max_fuel_hydrostatic_pressure_guess; % psid

%% output data

% write names to spreadsheet %

pretank_names = cell(size_mat_pretank_losses(1),1);
pretank_names(1:2,1) = PreTank_Line_Info_Table{1:2,1};
pretank_names(3,1) = PreTank_Component_Info_Table {1,1};
pretank_names(4,1) = PreTank_Line_Info_Table{3,1};
pretank_names(5,1) = {'HX'};
pretank_names(6,1) = PreTank_Line_Info_Table{4,1};

pretank_ox_names = cell(size_mat_ox_pretank_losses(1),1);
j = 1;
for i = 1:2:6
    pretank_ox_names(i,1) = Ox_PreTank_Line_Info_Table{j,1};
    pretank_ox_names(i+1,1) = Ox_PreTank_Component_Info_Table{j,1};
    j = j + 1;
end
pretank_ox_names(size_mat_ox_pretank_losses(1),1) = Ox_PreTank_Line_Info_Table{j,1};

pretank_fuel_names = cell(size_mat_fuel_pretank_losses(1),1);
j = 1;
for i = 1:2:6
    pretank_fuel_names(i,1) = Fuel_PreTank_Line_Info_Table{j,1};
    pretank_fuel_names(i+1,1) = Fuel_PreTank_Component_Info_Table{j,1};
    j = j + 1;
end
pretank_fuel_names(size_mat_fuel_pretank_losses(1),1) = Fuel_PreTank_Line_Info_Table{j,1};

posttank_ox_names = cell(size_mat_ox_posttank_comp(1)+size_mat_ox_posttank_line(1),1);
posttank_ox_names(1:size_mat_ox_posttank_comp(1),1) = Ox_PostTank_Component_Info_Table{1:size_mat_ox_posttank_comp(1),1};
posttank_ox_names(1+size_mat_ox_posttank_comp(1):size_mat_ox_posttank_comp(1)+size_mat_ox_posttank_line(1),1) = Ox_PostTank_Line_Info_Table{1:size_mat_ox_posttank_line(1),1};

posttank_fuel_names = cell(size_mat_fuel_posttank_comp(1)+size_mat_fuel_posttank_line(1),1);
posttank_fuel_names(1:size_mat_fuel_posttank_comp(1),1) = Fuel_PostTank_Component_Info_Table{1:size_mat_fuel_posttank_comp(1),1};
posttank_fuel_names(1+size_mat_fuel_posttank_comp(1):size_mat_fuel_posttank_comp(1)+size_mat_fuel_posttank_line(1),1) = Fuel_PostTank_Line_Info_Table{1:size_mat_fuel_posttank_line(1),1};

ox_names = [pretank_names; pretank_ox_names; posttank_ox_names];
fuel_names = [pretank_names; pretank_fuel_names; posttank_fuel_names];

writecell(ox_names,filename,'Sheet',11,'Range','A2');
writecell(fuel_names,filename,'Sheet',12,'Range','A2');

% write values to file %

% ox side
size_mat_ox_names = size(ox_names);
ox_losses = zeros(size_mat_ox_names(1),7);

ox_losses(1:size_mat_pretank_losses(1),1:size_mat_pretank_losses(2)) = pretank_losses;
ox_losses(size_mat_pretank_losses(1)+1:size_mat_pretank_losses(1)+size_mat_ox_pretank_losses(1),...
    1:size_mat_ox_pretank_losses(2)) = ox_pretank_losses;
ox_losses(size_mat_pretank_losses(1)+size_mat_ox_pretank_losses(1)+1:size_mat_pretank_losses(1)+size_mat_ox_pretank_losses(1)...
    +size_mat_ox_posttank_comp(1),3) = ox_posttank_component_losses;
ox_losses(size_mat_pretank_losses(1)+size_mat_ox_pretank_losses(1)+size_mat_ox_posttank_comp(1)+1:size_mat_pretank_losses(1)+size_mat_ox_pretank_losses(1)...
    +size_mat_ox_posttank_comp(1)+size_mat_ox_posttank_line(1),3) = ox_posttank_line_losses;

writematrix(ox_losses,filename,'Sheet',11,'Range','B2');

% fuel side 
size_mat_fuel_names = size(fuel_names);
fuel_losses = zeros(size_mat_fuel_names(1),7);

fuel_losses(1:size_mat_pretank_losses(1),1:size_mat_pretank_losses(2)) = pretank_losses;
fuel_losses(size_mat_pretank_losses(1)+1:size_mat_pretank_losses(1)+size_mat_fuel_pretank_losses(1),...
    1:size_mat_fuel_pretank_losses(2)) = fuel_pretank_losses;
fuel_losses(size_mat_pretank_losses(1)+size_mat_fuel_pretank_losses(1)+1:size_mat_pretank_losses(1)+size_mat_fuel_pretank_losses(1)...
    +size_mat_fuel_posttank_comp(1),3) = fuel_posttank_component_losses;
fuel_losses(size_mat_pretank_losses(1)+size_mat_fuel_pretank_losses(1)+size_mat_fuel_posttank_comp(1)+1:size_mat_pretank_losses(1)+size_mat_fuel_pretank_losses(1)...
    +size_mat_fuel_posttank_comp(1)+size_mat_fuel_posttank_line(1),3) = fuel_posttank_line_losses;

writematrix(fuel_losses,filename,'Sheet',12,'Range','B2');

% write other outputs to file %

% ADD REQUIRED OX, FUEL, TOTAL PROP MASS AND TOTAL DRY MASS (COPV TANKS +
% PRESSURANT)
other_outputs = zeros(27,1);
other_outputs(1,1) = pretank_losses(1,4); % F, pre reg temp
other_outputs(2,1) = pretank_losses(1,1); % psia, pre reg pressure
other_outputs(3,1) = pretank_losses(1,1) - 14.7; % psig, pre reg pressure
other_outputs(4,1) = reg_set_pressure; % psia
other_outputs(5,1) = reg_set_pressure - 14.7; % psig
other_outputs(6,1) = min_blowdown_pressure; % psia
other_outputs(7,1) = min_blowdown_pressure - 14.7; % psia
other_outputs(8,1) = ox_pressurant_mdot; % lbm/s, required ox pressurant mass flow rate
other_outputs(9,1) = fuel_pressurant_mdot; % lbm/s, required fuel pressurant mass flow rate
other_outputs(10,1) = total_pressurant_mdot; % lbm/s, required total pressurant mass flow rate
other_outputs(11,1) = nominal_pressurant_mdot; % lbm/s, nominal HX mdot using bypass valve to bring outlet temp from 400F to 200F
other_outputs(12,1) = GG_total_mdot; % lbm/s
other_outputs(13,1) = pretank_losses(HXcount,4); % F, HX inlet temp
other_outputs(14,1) = HX_tube_length; % in, length per tube
other_outputs(15,1) = space_btwn_tubes; % in, radial spacing between tubing coils 
other_outputs(16,1) = ox_required_tank_pressure; % psia
other_outputs(17,1) = ox_required_tank_pressure - 14.7; % psig
other_outputs(18,1) = fuel_required_tank_pressure; % psia
other_outputs(19,1) = fuel_required_tank_pressure - 14.7; % psig
other_outputs(20,1) = ox_pressurant_mass; % lbm
other_outputs(21,1) = fuel_pressurant_mass; % lbm
other_outputs(22,1) = total_pressurant_mass; % lbm
other_outputs(23,1) = total_num_COPVs_ox_tank;
other_outputs(24,1) = total_num_COPVs_fuel_tank;
other_outputs(25,1) = orifice_diameter; % in, trim orifice diameter
other_outputs(26,1) = orifice_dP; % psid, pressure drop required across trim orifice
other_outputs(27,1) = burn_time; % seconds

writematrix(other_outputs,filename,'Sheet',13,'Range','B2');
writecell(orifice_side,filename,'Sheet',13,'Range','B29');






