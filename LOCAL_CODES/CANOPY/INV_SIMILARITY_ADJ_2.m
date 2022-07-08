function [ Conc_can ] = SIMILARITY_ADJ(FORCING, PARAMS, CONSTANTS, VERTSTRUC, H, E, Conc_top, case_var)

%=========================================================================
% Use similarity to calculate water, heat, windspeed at the top of the
% canopy using water, heat, windspeed at the tower height
%
% Rewritten for multi-direction by Meredith Richardson - 09.27.2020
%
% Written By: Esther Lee with help from Meredith Richardson & Leila
% Hernandez
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       FORCING         % FORCING structure
%       VERTSTRUC       % VERTSTRUC structure
%       PARAMS          % PARAMS structure
%------------------------- Output Variables ------------------------------
%       U               % [m/s] Wind speed 
%       km              % Km momentum diffusivity 
%       tau             % parameter momentum model
%       Conc_can        % [vary] Concentration of tracer considered 
%                         at the top of the canopy
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************
    Utop = FORCING.U;                                                      %[m/s] Wind Speed at the top (Forcing from atmospheric measurements)  
    Ustar = FORCING.Ustar;                                                 %[m/s] Friction Velocity
    hcan = PARAMS.CanStruc.hcan;                                           %[m] Height canopy  
%     h_toplayer = VERTSTRUC.znc(20); %Mere
    hobs = PARAMS.CanStruc.hobs;                                           %[m] Height flux tower %added by Esther to use for Utop calculation
    vonk = CONSTANTS.vonk;                                                 % [] von Karman constant 
    %case_var;                                                              %1. for U  2. for T (assume it is the same as potential T)   3. for specific humidity
    z0 = PARAMS.Soil.z0;                                                   % [m] surface roughness length
    rho = CONSTANTS.rho_dry_air;                                           %density of dry air [kg/m3] check if is is the right rho!!!
    cp = CONSTANTS.cp_JkgK;                                                   %specific heat of air at constant pressure [J/kg/K]
    h_layer = VERTSTRUC.znc;                                               %[m] height of center of all layers
    fLAIz = VERTSTRUC.fLAIz;
 
%*************************************************************************

 h_tower=PARAMS.inv_hobs;    %height of alpine/actice flux tower
 ind_sim = find(~fLAIz); ind_canturb = find(fLAIz>0);
 h_new = h_layer';
 
d0 = PARAMS.inv_hcan*(2/3); 
Conc_can = Conc_top;

if case_var==1                  %case 1. for U      added by Esther
    Conc_can(ind_sim) = Conc_top(ind_sim) + Ustar/vonk * log((h_new(ind_sim)-d0) / (h_tower-d0));
%     Conc_can = Conc_top + Ustar/vonk * log((h_toplayer-d0) / (hobs-d0));
    
elseif case_var==2              %case 2. for T (assume it is the same as potential T)       added by Esther  
    Conc_can(ind_sim) = Conc_top(ind_sim) - H/vonk/Ustar/rho/cp * log((h_new(ind_sim)-d0) / (h_tower-d0));
%     Conc_can = Conc_top - H/vonk/Ustar/rho/cp * log((h_toplayer-d0) / (hobs-d0));

elseif case_var==3              %case 3. for specific humidity      added by Esther
%     Conc_can(ind_sim) = Conc_top(ind_sim) - E./vonk./Ustar./rho .* log((h_new(ind_sim)-d0) ./ (h_tower-d0));
    Conc_can(ind_sim) = Conc_top(ind_sim) - E(ind_sim)./vonk./Ustar./rho.* log((h_new(ind_sim)-d0) ./ (h_tower-d0)); %works for alpine trees
%     Conc_can = Conc_top - E/vonk/Ustar/rho * log((h_toplayer-d0) / (hobs-d0));

end

Conc_can(ind_canturb) =  Conc_can(ind_sim(1));
