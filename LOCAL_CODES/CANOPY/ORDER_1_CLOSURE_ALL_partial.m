function [Ca_out] = ORDER_1_CLOSURE_ALL_partial ( Ca_in, znc, dzc, K_in, SS_in, soilf, VERTSTRUC)

%=========================================================================
% Calls a portion of the ORDER_1_CLOSURE_ALL routine for canopy that exists
% in part of the control volume height
%
% Written By: Meredith Richardson; 9/27/2020
% All rights reserved!
%
% 
%========================================================================

%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************

 fLAIz = VERTSTRUC.fLAIz;
 ind_canturb = find(fLAIz>0);
 if(length(ind_canturb <2))
     ind_canturb2 = [1; 2];
 else
     ind_canturb2 = ind_canturb;
 end
 hcan_real = znc(max([ind_canturb; 2]));
  Ca_out = Ca_in;
 
 [Ca_in(ind_canturb2)] = ORDER_1_CLOSURE_ALL ( Ca_in(ind_canturb2), znc(ind_canturb2), dzc, ...
     K_in(ind_canturb2), SS_in(ind_canturb2),soilf, hcan_real );
 
 Ca_out(ind_canturb) = Ca_in(ind_canturb);
 
 %[Conc] = ORDER_1_CLOSURE_ALL ( Ca_in, z_in, dz, K_in, SS_in, soilf, ctz )