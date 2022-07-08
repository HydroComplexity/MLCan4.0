function [Ph, An, Ph_limit, Jc, Jj, Js, photo] = ...
    PHOTOSYNTHESIS_C3(VARIABLES, PARAMS, VERTSTRUC, CONSTANTS, SWITCHES, sunlit,cntspecies) 

%=========================================================================
%   This code solves the leaf dynamics 
%
% Written By: Darren Drewry, Modified by Juan Quijano & Meredith Richardson
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       VARIABLES        % VARIABLES structure
%       PARAMS           % VARIABLES structure
%       VERTSTRUC        % PARAMS structure
%       CONSTANTS        % VARIABLES structure
%       sunlit           % [1] Sunlit [0] Shade
%       cntspecies       % Number of species
%------------------------- Output Variables ------------------------------
%       Ph               % [umol CO2/ m^2 leaf / s] Photosynthetic flux from ecosystem 
%       An               % [umol CO2/ m^2 leaf / s] Photosynthetic minus Leaf Respiration flux from ecosystem 
%       Ph_limit         % [] Type of Photosynthesis 1. Rubisco-limited,
%                               2.light-limited, 3. Sucrose Limited, 4. frozen
%       Jc               % [umol/m^2 leaf area/s] Rubisco-Limited Rate
%       Jj               % [umol/m^2 leaf area/s] Light-Limited Rate 
%       Js               % [umol/m^2 leaf area/s] Sucrose-Limited Rate 
%       Vcmax            % [umol/m^2 leaf area/s] Maximum rate of carboxylation 
%       Jmax             % [umol/m^2 leaf area/s] Maximum Electron transport rate 
%========================================================================              
%                          DE-REFERENCE BLOCK
%*************************************************************************
    % VARIABLES
    if (sunlit)  
        Qabs = VARIABLES.CANOPY.PARabs_sun;                                % [umol/m^2 leaf area/s] Absorbed PAR  
         % VARIABLES FOR EACH SPECIES 
         %*****************************************************************
        Tl =  VARIABLES.CANOPY.Tl_sun(:,cntspecies);                       % [C] Leaf Temperature  
        Ci =  VARIABLES.CANOPY.Ci_sun(:,cntspecies);                       % [umol/mol] Atmospheric Concentration of CO2 in the leaf 
         %*****************************************************************        
    else
        Qabs = VARIABLES.CANOPY.PARabs_shade;                              % [umol/m^2 leaf area/s] Absorbed PAR  
         % VARIABLES FOR EACH SPECIES 
         %*****************************************************************
        Tl =  VARIABLES.CANOPY.Tl_shade(:,cntspecies);                     % [C] Leaf Temperature  
        Ci =  VARIABLES.CANOPY.Ci_shade(:,cntspecies);                     % [umol/mol] Atmospheric Concentration of CO2 in the leaf 
         %*****************************************************************        
    end
    timestep = VARIABLES.timestep;
    Vz = VARIABLES.CANOPY.Vz;                                              % [] Vertical distribution of photosynthetic capacity     
  % variables that are unique for each species  
  %************************************************************************
    vec = PARAMS.Photosyn.Vcmax25_C3{cntspecies};                          % [umol/m^2 leaf area/s] Vector of Vcmax at 25 C  
    Vcmax25 = vec(timestep);                                               % [umol/m^2 leaf area/s] Vcmax at timestep 
    vec = PARAMS.Photosyn.Jmax25_C3{cntspecies};                           % [umol/m^2 leaf area/s] Vector of Jmax at 25 C 
    Jmax25 = vec(timestep);                                                % [umol/m^2 leaf area/s] Jmax at timestep  
    vec = PARAMS.Photosyn.Rd25{cntspecies};                                % [umol/m^2 leaf area/s] Vector of Rd at 25 C  
    Rd25 = vec(timestep);                                                  % [umol/m^2 leaf area/s] Rd at timestep  
    beta = PARAMS.Photosyn.beta_ph_C3(cntspecies);                         % [] Beta Parameter C3 Phothosyntesis model
    O = PARAMS.Photosyn.Oi(cntspecies);                                    % [mmol / mol] Intercellular oxygen concentration  
  %************************************************************************
    
  % VERTSTRUC  
    nvinds_all = VERTSTRUC.nvinds_all;                                     % [] Indicator of Layers where LAI is zero of lower for all species 
    nvinds = nvinds_all{cntspecies};                                       % [] Indicator of Layers where LAI is zero of lower  
    
  % CONSTANTS    
    R_J = CONSTANTS.R;                                                     % [J mol^-1 K^-1] R constant 
   
  % COLD REGIONS
  % Turn off Photosynthesis if Tl<....depth of snow > height of canopy - %Mere
%     hcan = PARAMS.CanStruc.hcan;                                           % [m] canopy height %Mere
%     zsn = VARIABLES.SOIL.zsn;                                              % [mm] depth of snow
photo = VARIABLES.CANOPY.photo;
meanTAday = VARIABLES.CANOPY.meanTAday; %calculated 
ph_minT = PARAMS.Photosyn.ph_minT;
ph_restartT = PARAMS.Photosyn.ph_restartT;
Ts_min = min(VARIABLES.SOIL.Ts); 
Ts_top = VARIABLES.SOIL.Ts(1);
% if Ts_top>Ts_min
%     Ts_min
% end
%*************************************************************************
%*************************************************************************
    if (SWITCHES.coldregion) %Check & update photosynthesis switches
        if (photo==0 && meanTAday > ph_restartT && VARIABLES.SOIL.zsn/1000<PARAMS.true_hcan...
                && timestep>PARAMS.photottlimit(1) && timestep<PARAMS.photottlimit(2)...
                && Ts_top>0)%NR!: 5000, 14000
            photo = 1
            fprintf('timestep %d \n',timestep);
        elseif (photo==1 && ((meanTAday < ph_minT && (timestep<PARAMS.photottlimit(3)...
                || timestep>PARAMS.photottlimit(4))) || VARIABLES.SOIL.zsn/1000>PARAMS.true_hcan...
                || Ts_top<-1))%NR1:8500,12000
            photo = 0
            fprintf('timestep %d \n',timestep);
        end
    end
 
    if (Qabs<=0)
        Ph = zeros(size(Tl));
        An = zeros(size(Tl));
        Ph_limit = zeros(size(Tl));
        gamstar = zeros(size(Tl));
        Wc = zeros(size(Tl));
        Wj = zeros(size(Tl));
    end
    
 
    R = R_J/1000;    % [kJ mol^-1 K^-1] 
    
    TlK = Tl + 273.15;  
    
    Vcmax25 = Vz * Vcmax25;
    Jmax25 = Vz * Jmax25;
    
    gamstar = exp(19.02 - 37.83./(R*TlK));
    Ko = exp(20.30 - 36.38./(R*TlK));
    Kc = exp(38.05 - 79.43./(R*TlK));
    
    Rd = Rd25 * exp(18.72 - 46.39./(R*TlK));
                
    Vcmax = Vcmax25 .* exp(26.35 - 65.33./(R*TlK));
     
    phiPSIImax = 0.352 + 0.022*Tl - 0.00034*Tl.^2;
    
    thetaPSII = 0.76 + 0.018*Tl - 0.00037*Tl.^2;
    
        
    Q2 = Qabs .* phiPSIImax * beta;

    
    Jmax = Jmax25 .* exp(17.57 - 43.54./(R*TlK));
    J = (Q2 + Jmax - sqrt((Q2+Jmax).^2 - 4*thetaPSII.*Q2.*Jmax)) ./ (2*thetaPSII);
    
    
        % Check phiPSII and thetaPSII
    
    if (min(phiPSIImax)<0 | min(thetaPSII)<0)
        ind1 = phiPSIImax<0;
        phiPSIImax(ind1)=0;
        ind2 = thetaPSII<0;
        thetaPSII(ind2)=0;
        J(ind1)=0;
        J(ind2)=0;
    end   
    

    
    Wc = (Vcmax .* Ci) ./ (Ci + Kc.*(1+O./Ko));
    Wj = (J .* Ci) ./ (4.5*Ci + 10.5*gamstar);
    
    % Limiting Rates
    Jc = (1 - gamstar./Ci) .* Wc;
    Jj = (1 - gamstar./Ci) .* Wj;
    Js = Vcmax ./ 2;
    
    % Solve quadratics from [Collatz et al, AFM 1991] to account for co-limitation between rates
    tt = 0.98;
    bb = 0.96;
    Jp = ( (Jc+Jj) - sqrt( (-(Jc+Jj)).^2 - 4*tt*Jc.*Jj ) ) ./ (2*tt);
    Ph = ( (Jp+Js) - sqrt( (-(Jp+Js)).^2 - 4*bb*Jp.*Js ) ) ./ (2*bb);
    
    %**********************************************************************
    indph = (Ci<=gamstar & Qabs > 0);
    Ph(indph) = Rd(indph);
    
if (SWITCHES.coldregion && ~photo) %Mere
        Ph = zeros(size(Tl));
        Ph_limit = 4.*ones(size(Tl)); %photolimit is temperature
else 
    Ph_limit = NaN(size(Jc));
    Jcinds = find(Jc<Jj & Jc<Js);
    Ph_limit(Jcinds) = 1;
    Jjinds = find(Jj<=Jc & Jj<=Js);
    Ph_limit(Jjinds) = 2;
    Jsinds = find(isnan(Ph_limit));
    Ph_limit(Jsinds) = 3;
end     

    An =  Ph - Rd;
    
    Ph(isnan(Tl)) = 0;
    An(isnan(Tl)) = 0;
                