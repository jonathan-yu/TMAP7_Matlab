function [ outPUT ] = Perform_Fit_PTTP_v4( inPUT, vChng, pChng )
% Perform_Fit_Pseudo runs TMAP_PTTP_v4 through multiple cycles
%   inPUT is the structure used in TMAP_PTTP_v4
%   nChng is the number of steps away from central value
%   vChng point to the variables to be changed
%       i.e. 1 == surface peak concentration (due to trap 1:2 only)
%            2 == damage peak concentration (3:6 only heavy ion damage)
%            3 == intrinisic trap concentration
%            4 == E_trap (trapping energy)
%                0.9   1.1   1.4   1.7   1.9   2.1   [eV]
%                425   500   640   730   840   940   [K]
%   pChng is the percent change for each variable in vChng
%   resetFit 1 or 0 (t or f) to determine if the calculation should be redone
%       for the central value (i.e. if nChng = 2 and reset = 0
%       then variations 1,2,4,and 5 performed and 3 would have no
%       variation) Only use reset = 0, if inPUT has not been adjusted in
%       the workspace
%   refineFit 1 or 0 (t or f) to determine if a refined fit will be performed
%

nTraps = inPUT.nTraps;          % Determine number of traps

% Surface___Traps = 1:2;          % index/variable 1 for zone 1
% Damage____Traps = 3:nTraps;     % index/variable 2 for zone 2
% Intrinsic_Traps = 1:nTraps;     % index/variable 3 for zone 3
% Energy____Traps = 1:nTraps;     % index/variable 4

% version 5 to include check if variable is set to zero, if so then ignore
Surface___Traps = 2:3;          % index/variable 1 for zone 1
Damage____Traps = 3:nTraps;     % index/variable 2 for zone 2
Intrinsic_Traps = 1:nTraps;     % index/variable 3 for zone 3
Energy____Traps = 1:nTraps;     % index/variable 4

pChng = pChng/100;              % modify to represent fractional delta change

if inPUT.fitP.loop.resetFit == 1
    inPUT.fitP.tot = 0;         % reset to zero to perform center calculation
end
    
for ii = 1:inPUT.fitP.loop.outter  
 for mm = 1:length(vChng)
  iChng1 = vChng(mm);
  dChng1 = pChng(mm);  
  % Point to traps associated with an index
  if     iChng1 == 1
      tChng1 = Surface___Traps;
  elseif iChng1 == 2 
      tChng1 = Damage____Traps;
  elseif iChng1 == 3
      tChng1 = Intrinsic_Traps;
  elseif iChng1 == 4 
      tChng1 = Energy____Traps;
  end
      
  dChng = dChng1*0.90^(ii-1);       % if kk>1, then outter loop is reduced to 90%
  disp([' | i = ' num2str(ii,'%3d') ...
        ' | Delta Change = ' num2str(dChng*100,'%9.7f') ' %'])  
  inPUT.fitP.fracChng = dChng;      % fractional delta step
  for jj=1:length(iChng1)           % Loop over all indexes to be fit
      jValue =  iChng1(jj);         % Point to the index
      for kk = 1:length(tChng1)     % Loop over all traps in a particular index
          kValue = tChng1(kk);      % Point to the trap
          inPUT.fitP.trapChng = kValue;
          inPUT.fitP.indxChng = jValue;
          inPUT = TMAP_PTTP_v4(inPUT);                
          save('Workspace_Fit', 'inPUT');
          disp('Saved inPUT');
          if inPUT.fitP.loop.refineFit == 1
              inPUT = Refine_Fit(inPUT,1);
          end
      end
    outPUT = inPUT;
  end 
    
end
 
end
end

function [inPUT] = Refine_Fit(in_,nAdjust)

% Save the old values to reset at the end
old.nn_Chng  = in_.fitP.nn_Chng;
old.fracChng = in_.fitP.fracChng;

% Decide how to adjust and refine the fit 
in_.fitP.nnChng = in_.fitP.nn_Chng + nAdjust;
in_.fitP.fracChng = in_.fitP.fracChng/(in_.fitP.nnChng*2);

inPUT = TMAP_PTTP_v4(in_);                
inPUT.fitP.nn_Chng = old.nn_Chng;  
inPUT.fitP.fracChng = old.fracChng;

save('Workspace_Fit', 'inPUT');
disp('Saved inPUT');
end
    