function [VERTSTRUC] = SOIL_PROPERTIES(PARAMS, VERTSTRUC, SWITCHES)

%=========================================================================
% This function calculates the profiles of the soil properties for Arctic
% peat beyond the soil/clay percentages
%
% See Oleson et al (2004) (CLM Documentation) for original equation # references
%
% Peat Equations: Zhao et al, 2019 [doi: 10.1016/j.agrformet.2019.04.004]
% Parameterized to Canadian Arctic from Wu et al, 2016 [10.5194/gmd-9-2639-2016]
% & Krogh et al, 2017 [10.1016/j.jhydrol.2017.05.042]
%
% Written by Meredith Richardson (2020) by editing "SOIL_PROPERTIES" written
% by Darren Drewry
%=========================================================================

% Dereference Structure Values
    nl_soil = VERTSTRUC.nl_soil;

    sand = VERTSTRUC.sand*ones(nl_soil,1);
    clay = VERTSTRUC.clay*ones(nl_soil,1);
%    load sand.mat;
%    load clay.mat;
    
    zhs = VERTSTRUC.zhs;
    smpmin = PARAMS.Soil.smpmin;
    scalek = PARAMS.Soil.scalek;
    
    hemic_peat=find(zhs<0.10); %varies by location
    sapric_peat=find(zhs>0.10);
    

% HEAT CAPACITY OF SOIL SOLIDS [J / m^3 / K]
    HC_sol = 2.5 * 10^6.*ones(size(zhs));     %Wu et al, 2016
    
% POROSITY = WATER CONTENT AT SATURATION [-]
    porsl(hemic_peat,1) = 0.825; %Krogh et al 2017
    porsl(sapric_peat,1) = 0.775;
     
% MINIMUM SOIL SUCTION = SOIL POTENTIAL AT SATURATION [mm]
    psi0 = -10 * ( 10.^(1.88-0.0131*sand) );                    % (7.75)
%     psi0(hemic_peat,1) = 10.2;                    % Wu et al, 2016
%     psi0(sapric_peat,1) = 10.1;
     
% Clapp-Hornberger "b" parameter [-]
    bsw(hemic_peat,1) = 6.1;                    % Wu et al, 2016
    bsw(sapric_peat,1) = 12.0;
    
% THERMAL CONDUCTIVITY OF SOIL MINERALS [W / m / K]
    TK_min = 7.7; TK_org = 0.25; 
    f_min(hemic_peat,1)=0.2;f_min(sapric_peat,1)=0.1;
    f_org = 1-f_min;
    TK_sol =  (TK_min.^(f_min)) .* (TK_org.^(f_org)); %[Eqn 3]
     
% BULK DENSITY [kg / m^3]
    rhod = 2700*(1 - porsl);                                    % (before 6.62)
    
% THERMAL CONDUCTIVITY OF DRY SOIL [W / m / K]
    TK_dry = 0.30.*10.^(-0.87.*porsl); %[Zhao et al, 2019]
     
% HYDRAULIC CONDUCTIVITY AT SATURATION [mm / s] 
    HKsat(hemic_peat,1) = 1.04e-03; %top 10cm upper peat
    HKsat(sapric_peat,1) = 9.88e-04; %top 10cm upper peat

% SOIL MOISTURE CONTENT AT COMPLETE DRYNESS
    theta_dry = porsl .* ( psi0./smpmin ).^(1./bsw);            % 

   
% Dongkook: Start
% Soil moisture parameter for vanGenuchten
% Minasny and McBratney, 2007

% For n
%S1=1./(1+exp(-(24.547 - 0.238.*sand - 0.082.*clay)));
%S2=1./(1+exp(-(-3.569 + 0.081.*sand)));
%S3=1./(1+exp(-(0.694 - 0.024.*sand + 0.048.*clay)));

%n=2.18+0.11.*(48.087 - 44.954.*S1 - 1.023.*S2-3.896.*S3);

% For alpha
% parameter related to the inverse of the air entry suction [1/cm] to [1/mm]
% Markus Tuller and Dani Or, 2003
% http://www.mathworks.com/matlabcentral/fileexchange/45468-soil-classification/content/soil_classification.m
sandp=sand./100;
clayp=clay./100;

siltp=1-sandp-clayp;
for i=1:length(clayp)
    if (siltp(i)+1.5*clayp(i))<.15
        %SC{i,:}='SAND';
        alpha(i) = 0.035*(1/10);
        n(i) = 3.19;
        thetar(i) = 0.058;
    elseif ((siltp(i)+1.5*clayp(i))>=.15)&((siltp(i)+2*clayp(i))<.3)
        %SC{i,:}='LOAMY SAND';
        alpha(i) = 0.035*(1/10);
        n(i) = 2.39;
        thetar(i) = 0.074;
    elseif (clayp(i)>=0.07) & (clayp(i)<=0.2) & (sandp(i)>0.52) & ((siltp(i)+2*clayp(i))>=0.3)
        %SC{i,:}='SANDY LOAM';
        alpha(i) = 0.021*(1/10);
        n(i) = 1.61;
        thetar(i) = 0.067;
    elseif (clayp(i)<0.07) & (siltp(i) < 0.5) & ((siltp(i)+2*clayp(i))>=0.3)
        %SC{i,:}='SANDY LOAM';
        alpha(i) = 0.021*(1/10);
        n(i) = 1.61;
        thetar(i) = 0.067;
    elseif (clayp(i)>=0.07) & (clayp(i)<=0.27) & (siltp(i)>=0.28) & (siltp(i)<0.5) & (sandp(i)<=0.52)
        %SC{i,:}='LOAM';
        alpha(i) = 0.036*(1/10);
        n(i) = 1.31;
        thetar(i) = 0.083;
    elseif ((siltp(i)>0.5) & clayp(i)>0.12 & clayp(i)<0.27) | (siltp(i)>=0.5 & siltp(i)<0.8 & clayp(i)<0.12)
        %SC{i,:}='SILT LOAM';
        alpha(i) = 0.025*(1/10);
        n(i) = 1.39;
        thetar(i) = 0.061;
    elseif siltp(i)>=0.8 & clayp(i)<0.12
        %SC{i,:}='SILT';
        alpha(i) = 0.006*(1/10);
        n(i) = 1.53;
        thetar(i) = 0.123;
    elseif clayp(i)>=0.2 & clayp(i)<0.35 & siltp(i)<0.28 & sandp(i)>0.45
        %SC{i,:}='SANDY CLAY LOAM';
        alpha(i) = 0.033*(1/10);
        n(i) = 1.49;
        thetar(i) = 0.086;
    elseif clayp(i)>=0.27 & clayp(i) <0.4 & sandp(i)>0.2 & sandp(i)<=0.45
        %SC{i,:}='CLAY LOAM';
        alpha(i) = 0.030*(1/10);
        n(i) = 1.37;
        thetar(i) = 0.129;
    elseif clayp(i)>=0.27 & clayp(i)<0.4 & sandp(i)<=0.2
        %SC{i,:}='SILTY CLAY LOAM';
        alpha(i) = 0.027*(1/10);
        n(i) = 1.41;
        thetar(i) = 0.098;
    elseif clayp(i)>=0.35 & sandp(i)>45
        %SC{i,:}='SANDY CLAY';
        alpha(i) = 0.025*(1/10);
        n(i) = 1.4;
        thetar(i) = 0.08;
    elseif clayp(i)>=0.4 & siltp(i)>=0.4
        %SC{i,:}='SILTY CLAY';
        alpha(i) = 0.023*(1/10);
        n(i) = 1.39;
        thetar(i) = 0.163;
    elseif clayp(i)>= 0.4 & sandp(i)<=0.45 & siltp(i)<0.4
        %SC{i,:}='CLAY'
        alpha(i) = 0.021*(1/10);
        n(i) = 1.20;
        thetar(i) = 0.102;
    end
end
alpha=alpha(:);
n=n(:);
thetar=thetar(:);

%alpha(:)   = 0.01*(1/10);         % parameter related to the inverse of the air entry suction [1/cm] to [1/mm]
%n(:)       = 2.0;                  % Pore-size distributions [-]
% Dongkook: End

%*************************************************************************
% STORE IN STRUCTURE
%*************************************************************************

    % ASSIGN
        VERTSTRUC.HC_sol = HC_sol;
        VERTSTRUC.porsl = porsl;
        VERTSTRUC.psi0 = psi0;
        VERTSTRUC.bsw = bsw;
        VERTSTRUC.TK_sol = TK_sol;
        VERTSTRUC.TK_dry = TK_dry;
        VERTSTRUC.HKsat = HKsat;
        VERTSTRUC.theta_dry = theta_dry;        
        VERTSTRUC.eff_poros = VERTSTRUC.porsl;
        % Dongkook Woo - Edit 
        VERTSTRUC.rhod= rhod;
        VERTSTRUC.clayColum = clay;
        
        VERTSTRUC.VanGen_n=n;
        VERTSTRUC.VanGen_alpha=alpha;
        VERTSTRUC.VanGen_thetar=thetar;
        % Dongkook Woo - Edit End
    
    
    