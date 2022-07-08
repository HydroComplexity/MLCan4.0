function [ Conc_can ] = SIMILARITY_ADJ(FORCING, PARAMS, CONSTANTS, VERTSTRUC, H, E, Conc_top, case_var)

%=========================================================================
% Use similarity to calculate water, heat, windspeed at the top of the
% canopy using water, heat, windspeed at the tower height
%
% Edited for opposite direction by Meredith Richardson - 03.15.2020
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

%*************************************************************************

 hhigh=hcan;   %height of trees %8
 hlow=PARAMS.inv_hobs;    %height of alpine/actice flux tower
d0_2 = hhigh*(2/3); 
d0_1 = PARAMS.inv_hcan*(2/3); 
% d0 = h_toplayer*2/3;                                                                %Esther: zero plane displacement
% Ustar = Utop*vonk/(log((hobs-d0)/z0));                                        %Esther: friction velocity

% if H>0
%     H
% end

if case_var==1                  %case 1. for U      added by Esther
    Conc_can = Conc_top + Ustar/vonk * log((hhigh-d0_2) / (hlow-d0_1));
%     Conc_can = Conc_top + Ustar/vonk * log((h_toplayer-d0) / (hobs-d0));
    
elseif case_var==2              %case 2. for T (assume it is the same as potential T)       added by Esther  
    Conc_can = Conc_top - H/vonk/Ustar/rho/cp * log((hhigh-d0_2) / (hlow-d0_1));
%     Conc_can = Conc_top - H/vonk/Ustar/rho/cp * log((h_toplayer-d0) / (hobs-d0));

elseif case_var==3              %case 3. for specific humidity      added by Esther
    Conc_can = Conc_top - E/vonk/Ustar/rho * log((hhigh-d0_2) / (hlow-d0_1));
%     Conc_can = Conc_top - E/vonk/Ustar/rho * log((h_toplayer-d0) / (hobs-d0));

end
