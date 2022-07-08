% function DRIVER_CROP_2_0
% MLCan Looper
% Susana Roque-Malo, 8/22/19
% Edited by Meredith Richardson Martin, 2021


% Instructions: First, save this code in the same directory as DRIVER_CROP_2_0.
% Then, open the MLCan Gui and set it up as if you were
% about to run 1 year on MLCan. Close the Gui. This step saves a new temp
% file with the correct filepaths for the crop and forcings you have
% chosen. The loop below will overwrite the temp file with updated initial
% conditions after each run. This uses a 365 day forcing file. If you want 
% to use the same forcing file over and over, leave the code as is. If you 
% want to use a different forcing file for each year of simulation, you
% will need to prepare the forcing files for each year and save them with a
% naming convention that matches the commented section below. Uncomment the
% section and it will call a different forcing file for each year.

clear all
% close all
clc

tempfilename = './Temps/temp_variable.mat';
load(tempfilename)
mtemp = matfile(tempfilename, 'Writable', true); %N
save ./Temps/temp_variable.mat -v7.3
% save(tempfilename,'-v7.3');

numyearstoloop = 2;
j=NaN(numyearstoloop);

for Y = 1:numyearstoloop
    save('Yloop.mat','Y','numyearstoloop')
    % UNCOMMENT THIS SECTION TO CHANGE FORCING FILE WITH EACH LOOP
    forcinglocation = 'D:\Users\Meredith\From Lyra\MLCan_MeredithEdits_working\MLCan4.0\Data\Italy\MBo_Alpine\2yr_run_forcing_trees\'; % add your file location here
    forcingfile = 'MBo_2yr_Forcing_trees_';
    sitestr = '.MBo_trees_';
    sitestr2 = 'MBo_trees_';
%     mtemp.fullpath_forcings = strcat(forcinglocation,forcingfile,num2str(Y));
    Y
    if Y>2
    %     mtemp = matfile(tempfilename, 'Writable', true); %N
        fullpath_forcings_new = strcat(forcinglocation,forcingfile,num2str(Y));
        forcinglocationshort = '.\Data\Italy\MBo_Alpine\2yr_run_forcing_trees\';
        working_forcings_new = strcat(forcinglocationshort,forcingfile,num2str(Y));
        tempfilename = './Temps/temp_variable.mat';
        load(tempfilename)
        fullpath_forcings = fullpath_forcings_new;
        working_forcings = working_forcings_new;
        clear fullpath_forcings_new working_forcings_new
        Y=Y+1;
        save(tempfilename,'-v7.3');
    elseif Y==2
    %     mtemp = matfile(tempfilename, 'Writable', true); %N
        fullpath_forcings_new = strcat(forcinglocation,forcingfile,num2str(Y));
        forcinglocationshort = '.\Data\Italy\MBo_Alpine\2yr_run_forcing_trees\';
        working_forcings_new = strcat(forcinglocationshort,forcingfile,num2str(Y));
        tempfilename = './Temps/temp_variable.mat';
        load(tempfilename)
        fullpath_forcings = fullpath_forcings_new;
        working_forcings = working_forcings_new;
        clear fullpath_forcings_new working_forcings_new
        save(tempfilename,'-v7.3');
        
    else
        mtemp.fullpath_forcings = strcat(forcinglocation,forcingfile,num2str(Y));
    end

% Canopy-Root-Soil-Atmosphere Exchange Model
%
%   Written By : Darren Drewry (dtd2@illinois.edu)
%              : Juan Quijano (quijano2@illinois.edu)
%              : Dongkook Woo (dwoo5@illinois.edu)

% STRUCTURES:
%   SWITCHES --> model conditional switches
%   VERTSTRUC --> variables describing vertical structure of canopy & soil
%   VARIABLES.CANOPY --> canopy variables
%            .SOIL --> soil variables
%            .ROOT --> root variables
%   FORCING --> holds current timestep forcing variables
%   CONSTANTS --> site independent constants, unit conversions
%   PARAMS
%         .CanStruc --> canopy structural parameters
%         .Rad --> radiation parameters
%         .Photosyn --> photosynthesis paramters
%         .StomCond --> stomtatal conductance parameters
%         .Resp --> ecosystem respiration parameters
%         .MicroEnv --> canopy microenvironment parameters
%         .Soil --> soil paramters
%
%**************************************************************************
%                           USER SPECIFICATIONS
%**************************************************************************

load './Temps/temp_variable.mat' ...
    'DOY_start' 'DOY_end'
doys = [DOY_start:DOY_end];

%******************************  SWITCHES   *******************************
load './Temps/temp_variable.mat'...
    'Turbulence' 'set_para_root1'  'vanGen' 'RHC' 'Soil_heat' ...
    'Temp_Elev' 'Temp_Elev_con' 'CO2_Ambient' 'CO2_Elev'  ...
    'CO2_Elev_con' 'Soil_nutrient' 'litter_depth'

SWITCHES.plants = 1;              %  Allow the presence or no of plants. 1. With plants 2. No plants (bare soil, no roots)
SWITCHES.arcticpeat = 0;          %  Arctic soil conditions: permafrost, peated soils
SWITCHES.coldregion = 1;          %  Turn off photosynthesis with frozen soil/deep snow
PARAMS.Photosyn.ph_minT = -7;%-3; %-7
PARAMS.Photosyn.ph_restartT = 3;%5; %3
PARAMS.photottlimit = [0, 18000, 18000, 0]; % photo hibernation limits [restart min, restart max, stop max, stop min]
% default: [0, 18000, 18000, 0]
PARAMS.refl_PAR_can_snow = 0.1;
PARAMS.refl_NIR_can_snow = 0.1;        %  95-mixforest or trans wood/shrub,mineral land;110-mixfor,peat;150-conif,peat
 % 80-broad,peat; 140-conif,mineral; 125-trans wood/shrub,peat; 170-peatbogs; 
 % 175-trans wood/shrub<10% cover; 65-broad,min(assumed); [Heinila et al, 2019; doi:10.1016/j.jag.2018.10.017] 
SWITCHES.Lv_T = 1;                %  Time variant Latent heat of vaporization; Eqn from Henderson-Sellers, 1984; doi: 10.1002/qj.49711046626
SWITCHES.replacing_on = 0; %1;
SWITCHES.invsimilarity = 1;
PARAMS.inv_hobs = 2.5; PARAMS.inv_hcan = 1; %[m]
PARAMS.true_hcan  = 36;
SWITCHES.LT = 0;                  %  Run long term dynamics with stochastically generated data 1. Yes 0. No
SWITCHES.turb_on = 1;% Turbulence;    %  1 = scalar profiles resolved, otherwise not resolved
if litter_depth > 0
    SWITCHES.litter = 1;          %  1 = Include litter dynamics. 0 does not include litter dynamics
else
    SWITCHES.litter = 0;          %  1 = Include litter dynamics. 0 does not include litter dynamics
end
SWITCHES.littertype = 1;          %  [1 LEs parameterized with aerodynamic resistance

%  2 LEs parameterized with vapor litter diffusivity]
SWITCHES.soilevap = 1;            %  1 = Include soil evaporation. 0 does not include soil evaporation
if set_para_root1{4,2} > 0
    SWITCHES.cutroots = 1;        %  0 if no cut at all
                                  %  (nlc) 1 if cut all the time same layers
                                  %  2 if cut with a threshold
                                  %  3 if cut with a embolism curve (PLC)
                                  %  4 combination of 1, 2 and 3
elseif set_para_root1{4,2} == 0
    SWITCHES.cutroots = 0;        %  0 if no cut at all
                                  %  (nlc) 1 if cut all the time same layers
                                  %  2 if cut with a threshold
                                  %  3 if cut with a embolism curve (PLC)
                                  %  4 combination of 1, 2 and 3
end

SWITCHES.rtcond  = 0;             %  0. Juan's Method, 1. Amenu Method
SWITCHES.rhc = RHC;               %  1 = linearly increasing root hydraulic conductivity with depth
SWITCHES.ns=1;                    %  Numerical Scheme. 1 = Use Implicit, 0 = Use explicit

SWITCHES.soilheat_on = Soil_heat; %  Compute the heat equation 1 = Yes, 0 = No
SWITCHES.canstorheat = 1;         %  Include storage of heat in leaves 1 = Yes, 0 = No
SWITCHES.useG_on=1;               %  Use of G measured in the soil energy balance instead of compute it 0 = Off, 1 = On
SWITCHES.useTs_on=0;              %  Use of Ts (soil temperature) [soil energy balance]

SWITCHES.entropy_on = 1;          %  Compute Entropy, 1 = Yes, 0 = No.
SWITCHES.entropymethod = 2;       %  Compute Entropy Method for LW entropy. 1. Landsberg 1979. 2. Smith 2001. 3. Clausius ***MEREDITH: 3 nonexistent in code; 7/26/2017***

SWITCHES.save_on = 1;             %  1 = save stored variables to .mat file, otherwise no save performed
SWITCHES.plots_on = 0;

SWITCHES.fsv_off = 0;             %  1 = hydraulic constraint turned OFF

SWITCHES.temp_change = Temp_Elev; %  0 = +0 C, 1 = 1 C,  2 = 2 C

% All of the CN model switch and parameters are in core_N
SWITCHES.soilCN_on = Soil_nutrient;%  Compute CN dynamics. 1 = Yes, 0 = No.
SWITCHES.Pedofunctions = 0;       %  1. Use pedo transfer functions with organic matter 0. Do not use
SWITCHES.vanGen = vanGen;         %  1. vanGenuchten, 0. Brooks Corey / Soil moisture chracteristic.

%**************************************************************************
% Code Library Paths
addpath('./LOCAL_CODES/CANOPY/');
addpath('./LOCAL_CODES/ENTROPY/');
addpath('./LOCAL_CODES/ROOT_SOIL/');
addpath('./LOCAL_CODES/ROOT_SOIL/IMPLICIT');
addpath('./LOCAL_CODES/ROOT_SOIL/CN_MODEL');
addpath('./LOCAL_CODES/NUMERICAL/');
addpath('./LOCAL_CODES/NUMERICAL/OTHERS');
addpath('./LOCAL_CODES/PLOTTING/');

%************************** LOAD INFORMATION  *****************************
year = nan;                       % Initialize year, otherwise it is recognized as function
LOAD_SITE_INFO;

ybeginds = nan(length(Each_year),1);
yendinds = nan(length(Each_year),1);

for yy = 1:length(Each_year)
    ybeginds(yy) = find(year==Each_year(yy), 1, 'first');
    yendinds(yy) = find(year==Each_year(yy), 1, 'last');
end
%************************** CLIMATE CHANGE ********************************
if (SWITCHES.temp_change == 0)    % Ambient Temperature
    temp_change = 0;
elseif SWITCHES.temp_change == 1  % Elevated Temperature
    temp_change = Temp_Elev_con;
    Ta_in = Ta_in + temp_change;
end

CO2base_elevated = CO2_Elev_con;
CO2base_ambient = str2num(CO2_Ambient);
if (CO2_Elev == 1)                % Elevated CO2
    CO2base = CO2base_elevated;
else                              % Ambient CO2
    CO2base = CO2base_ambient;
end
current_time=clock;
current_year=num2str(current_time(1));
current_month=num2str(current_time(2));
current_day=num2str(current_time(3));
if Y>1
    namefilein = ['./Result' '/Result_', num2str(num_species), 'species', '_Saved_on_' ...
    , current_month, '.', current_day, '.', current_year, sitestr, num2str(Y-1),'.mat'];
else
    
    namefilein = ['./Result/MBo_IC_2011.mat'];
%     namefilein = ['./Result' '/Result_', num2str(num_species), 'species',...
%         '_IC_TVan_', num2str(Y),'.mat'];
end

%************************ CANOPY - SOIL - FUNCTION ************************
CANOPY_SOIL_COUPLER;
%**************************************************************************

%******************************* SAVE FUNCTION ****************************
load './Temps/temp_variable.mat'...
    'num_species'
if isstr(num_species) == 1
    Number_S = num_species;
else
    Number_S = num2str(num_species);
end
current_time=clock;
current_year=num2str(current_time(1));
current_month=num2str(current_time(2));
current_day=num2str(current_time(3));
current_hour=num2str(current_time(4));
current_min=num2str(current_time(5));

% namefileout = ['./Result' '/Result_', Number_S, 'species', '_Saved_on_' ...
%     , current_month, '.', current_day, '.', current_year, '.', current_hour,...
%     '.', current_min,'.mat'];
namefileout = ['./Result' '/Result_', Number_S, 'species', '_Saved_on_' ...
    , current_month, '.', current_day, '.', current_year, sitestr, num2str(Y),'.mat'];
namefileout2 = ['./Result/Treeline Analysis/Italy/', sitestr2, num2str(Y+2011), '.mat'];

% Make dianuual averages
MAKE_AVERAGES;

% Show message %Mere uncommented this section --> line 159
if (length(Each_year)==1)
    timevect  = decdoy;
    timelabel = ['DOY: ', num2str(Each_year)];
else
%     timevect  = [1:length(decyear)];
          timevect  = decyear;
    timelabel = ['Decimal Year'];
end

if Y < numyearstoloop
%     overwrite
    disp('Updating Initial Conditions')
%     mtemp.nutrient_int(:,2) = VARIABLES.Cl(2:end);
%     mtemp.nutrient_int(:,3) = VARIABLES.Ch(2:end);
%     mtemp.nutrient_int(:,4) = VARIABLES.Cb(2:end);
%     mtemp.nutrient_int(:,5) = VARIABLES.Amm(2:end);
%     mtemp.nutrient_int(:,6) = VARIABLES.Nit(2:end);
%     mtemp.nutrient_int(:,7) = VARIABLES.CNl(2:end);
%     mtemp.nutrient_int_litter(:,2) = VARIABLES.Cl(1);
%     mtemp.nutrient_int_litter(:,3) = VARIABLES.Ch(1);
%     mtemp.nutrient_int_litter(:,4) = VARIABLES.Cb(1);
%     mtemp.nutrient_int_litter(:,5) = VARIABLES.Amm(1);
%     mtemp.nutrient_int_litter(:,6) = VARIABLES.Nit(1);
%     mtemp.nutrient_int_litter(:,7) = VARIABLES.CNl(1);
    mtemp.root_init(:,2) = volliq;
    mtemp.root_init(:,3) = Ts;
    mtemp.root_init_litter(:,2) = volliq(1);
    mtemp.root_init_litter(:,3) = Tli_store(end);%susana Ts(1)
%     Y
    disp('Simulation is done. File is saved in the Result folder');
    save (namefileout)
    save (namefileout2)
    clear all
    load('Yloop.mat')
else
    disp('This is the last year!!')
    msgbox(['Simulation is done.', namefileout(10:end), ' is saved in the Result folder'],'MLCan Simulation');
    save (namefileout)  
    save (namefileout2)  
end




end


%**************************************************************************


