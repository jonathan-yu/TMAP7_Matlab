function [outputStruct] = TMAP_PTTP_v4(inputStruct)
%% TMAP_PTTP_v4
%  For inputStruct see Make_trap_Implant_02
%
%  Writes input files and Runs TMAP7 with a LeapFrog scheme
%  Specify which trap parameter is to be adjusted by a fractional value
%  
%  Program Structure
%       Organize_Input()
%       Import_Depth_Grid()
%       % Perform Loop over adjusted variable to find best fit
%           % Implantation phase of TMAP Simulation
%           Implant_phase()
%               Calculate_traps()
%               Write_Input_File()      % See below for Struct
%               Run_Updated_TMAP() 
%               Separate_traps()
%           % Match sum of D filled traps to NRA spatial profile
%           Plot_traps()
%           % Perform Loop over LeapFrog transitions
%               % TDS phase of TMAP Simulation
%               TDS_phase()
%               Write_Input_File()      % See below for Struct
%               Plot_traps()
%               Run_Updated_TMAP()
%           % Match D surface flux to TDS temperature profile
%           Plot_TDS_TMAP()
%               get_Extracted_Data()    
%
%  Extract trap concentrations from TMAP output files
%       get_Current_traps()
%           load_traps()
%               parse_TMAP_out()
%
%  TMAP Subroutine Structure
%       Write_Input_File()
%           % Copy Sections 1,3 and 9 from TMAP template files
%           Write_2_main_____()
%           Write_4_thermal__()
%           Write_5_diffusion()
%           Write_6_equation_()
%           Write_7_table____()
%           Write_8_control__()
%               

% Initialize inputStruct prior to any calculations
inA = Organize_Input(inputStruct);       % Determine variable to be adjusted
sample = inputStruct.sample;             % Keep sample info separate

if strcmp(inA.Depth.source,'Import')
    [trapD0.depth] = Import_Depth_Grid();       % Get position and thickness grid
elseif strcmp(inA.Depth.source,'Create')
    [trapD0.depth] = Create_Depth_Grid(inA);    % Get position and thickness grid
end

% Empty vector to store goodness of fit for total number of Fits (2nn+1)
gof.pr = zeros(1,inA.nFits);             % NRA D profile fit
gof.rt = zeros(1,inA.nFits);             % TDS D retention fit
gof.sf = zeros(1,inA.nFits);             % TDS D surface flux fit

% Perform 2nn+1 Fits with a single variable adjusted +/- delta
for ii = 1:inA.nFits                                % Runs through the number of varied parameters
        inA.iFit = ii;
        if ii == inA.rptFit && inA.fitP.tot ~= 0    % Do not repeat TMAP run (save some time)
            gof.pr(inA.rptFit) = inA.fitP.prof;
            gof.sf(inA.rptFit) = inA.fitP.surf;
            gof.rt(inA.rptFit) = inA.fitP.retn;
%             disp('repeat value reused');                        
        else
            inA.iLeap = inA.nLeap+1;                                            
            inA.trap.Matrix(inA.fitP.trapChng,inA.fitP.indxChng) = inA.varV(ii);   % adjust variable in trapInfo
            [inA,trapD0] = Implant_phase(inA,trapD0);   
            trapD_ = trapD0;     % trapD0 is due to implant trapD_ releases
            % Write input file and run TMAP for nLeap transitions 
            gof.pr(ii) = Plot_NRA_TMAP(sample,inA,trapD_,30+ii);
            [inA,trapD_] = TDS_phase(sample,inA,trapD_,ii);
            % Plot and determine goodness of fit for full TDS run
            out = Plot_TDS_TMAP(sample,inA,trapD0,50+ii);
            gof.rt(ii) = out.gof_rt;
            gof.sf(ii) = out.gof_sf;
        end
end

TMAP_NRA.traps = trapD0;
TMAP_NRA.depth = trapD0.depth;
TMAP_TDS  = out.TMAP;

outputStruct = Organize_Output(inputStruct,gof,inA.varV);
outputStruct.sample.TMAP_NRA = TMAP_NRA;
outputStruct.sample.TMAP_TDS = TMAP_TDS;

end

%% Various Subroutines

function [outPUT] = Organize_Input(inPUT)
% Here the data from the input Vector is put into easily mainpulated tables
%   Specific sections relevant to writing TMAP input files are outlined

global folder
    
% Initialize variables and store in arrays as needed
outPUT.eqn  = inPUT.eqn;
outPUT.nTraps = length(inPUT.trap.z2_damg);

% Paramters to create a spatial grid focused on damage zone
outPUT.Depth = inPUT.Depth;

% Trap Profiles for zone 1 and 2 used in Calculate_traps
outPUT.trap.Profile_z1 = inPUT.trap.Profile_z1;
outPUT.trap.Profile_z2 = inPUT.trap.Profile_z2;

outPUT.control = inPUT.control;
outPUT.plot    = inPUT.plot;

% Parameters to be adjusted to fit NRA and TDS profiles
outPUT.trap.Matrix = zeros(outPUT.nTraps,4);    
outPUT.trap.Matrix(:,1) = inPUT.trap.z1_surf;    
outPUT.trap.Matrix(:,2) = inPUT.trap.z2_damg;    
outPUT.trap.Matrix(:,3) = inPUT.trap.z3_intr;    
outPUT.trap.Matrix(:,4) = inPUT.trap.E_trap;

outPUT.trap.Matrix(outPUT.nTraps+1,:) = 0;
implantFlux = inPUT.implant.Fluence/inPUT.implant.time;
outPUT.implantFlux = [                 0, implantFlux;...
                      inPUT.implant.time, implantFlux;...
                                       0,           0;...
                      inPUT.implant.down,           0];

outPUT.implantGauss.center = inPUT.implant.Gauss.center;
outPUT.implantGauss.width_ = inPUT.implant.Gauss.width_;
outPUT.implantGauss.offset = inPUT.implant.Gauss.offset;
outPUT.implantGauss.trans_ = inPUT.implant.Gauss.trans_;
                  
% Number of transitions in LeapFrog scheme
outPUT.nLeap = length(inPUT.TDS_LF.Temp_LF)+1;
outPUT.tstep = inPUT.control.tstep_TDS*ones(outPUT.nLeap,1);  
outPUT.tstep = [outPUT.tstep', inPUT.control.tstep_implant];          % Include Implant at end
outPUT.tstep_initial = inPUT.control.tstep_initial;

% Setup time begin, end, and continue
outPUT.Table1.timeBgn = zeros(outPUT.nLeap+1,1);                      % +1 for implantation
delta_Temp = [inPUT.TDS_LF.TempBgn,inPUT.TDS_LF.Temp_LF,inPUT.TDS_LF.TempEnd];
outPUT.Table1.TempBgn = delta_Temp(1:end-1);
outPUT.Table1.TempEnd = delta_Temp(2:end  );
delta_Temp = outPUT.Table1.TempEnd-outPUT.Table1.TempBgn;
outPUT.Table1.TempBgn = [outPUT.Table1.TempBgn, inPUT.implant.Temp];  % Include Implant at end
outPUT.Table1.TempEnd = [outPUT.Table1.TempEnd, inPUT.implant.Temp];  % Include Implant at end
outPUT.Table1.timeEnd = delta_Temp/inPUT.TDS_LF.ramprate;
outPUT.Table1.timeEnd = [outPUT.Table1.timeEnd, inPUT.implant.time];  % Include Implant at end

outPUT.Table1.timeCont = zeros(outPUT.nLeap+1,1);
for ii = 2:outPUT.nLeap+1
   outPUT.Table1.timeCont(ii) = sum(outPUT.Table1.timeEnd(1:ii-1));
end

outPUT.Table1.n_print = floor(outPUT.Table1.timeEnd./outPUT.tstep);

% Setup the ordering scheme outlined in description above
outPUT.order = zeros(outPUT.nLeap+1,3);
for ii = 1:outPUT.nLeap
    outPUT.order(ii,:) = circshift([ii, ii+1,outPUT.nTraps+1],[1,ii-1]);
end
outPUT.order(end,:) = outPUT.order(1,:);                    % Include Implant at end

outPUT.fitP   =    inPUT.fitP;
outPUT.nFits  = 2*outPUT.fitP.nn__Chng+1;               
outPUT.rptFit =   outPUT.fitP.nn__Chng+1;      

% Initialize parameters for keeping track of indicees
outPUT.iLeap = 1;                                           % Restart cycle
% Organize Folder indicees
folder.TMAP = inPUT.folder.TMAP;
folder.Base = inPUT.folder.Base; 
folder.Temp = inPUT.folder.Temp;                 
folder.Updt = inPUT.folder.Updt;               
folder.Work = inPUT.folder.Work;
folder.tmapName  = [inPUT.fileName '_M' ];                  % naming scheme i.e. M21_M
folder.tmapNames = cell(outPUT.nLeap,outPUT.nFits);         % Empty array for TMAP file names
for ii = 1:outPUT.nFits
    for kk = 1:outPUT.nLeap+1                               % +1 for implantation
        kj = mod(kk,outPUT.nLeap+1);                        % Use 00 for implantation phase
%         kj = kk;
        folder.tmapNames{kk,ii} = [folder.tmapName ...
                                   '_' num2str(kj,'%02d') ...
                                   '_' num2str(ii,'%02d')];
    end
end

trapChng = inPUT.fitP.trapChng;                             % which trap     to adjust
indxChng = inPUT.fitP.indxChng;                             % which variable to adjust
crntTrap = outPUT.trap.Matrix(trapChng,indxChng);           % Store current variable  
fracChng = inPUT.fitP.fracChng;                             % fractional change to variable
dltaChng = fracChng*crntTrap;                               % delta Change  
strtTrap  = crntTrap - outPUT.rptFit*dltaChng;              % Start value will be stepped

outPUT.varV = zeros(1,outPUT.nFits);                        % Empty vector to store varied quantity
% Fill array and matrix for varied parameter in fit and tmapNames
for ii = 1:outPUT.nFits
    outPUT.varV(ii) = strtTrap + ii*dltaChng;
end

outPUT.crntTrap =crntTrap;

outPUT.adjustP = 0;
    
outPUT.plot.YorN = inPUT.plot.YorN;

plot_res = 0;
if plot_res == 1
    % Plot Residence Time info
    Temp2 = 293:1293; 
    res_time = zeros(outPUT.nTraps+1,1001);

    for ii = 1:outPUT.nTraps
        res_time(ii,1:1001) = 1/outPUT.eqn.jmpF*exp(outPUT.trap.Matrix(ii,4)./(outPUT.eqn.k_B*Temp2'));
    end
    max1 = max(outPUT.trap.Matrix(3:end,2) , outPUT.trap.Matrix(3:end,3));
    Temp3 = ones(outPUT.nTraps-1,1)*Temp2;
    exp1 = exp(outPUT.trap.Matrix(3:end,4)*ones(1,1001)./(outPUT.eqn.k_B*Temp3));  
    expEp = sum(max1*ones(1,1001).*exp1/sum(max1));
    res_time(ii+1,1:1001) = 1/outPUT.eqn.jmpF*expEp;

    hold all
    legend_n = cell(outPUT.nTraps+1,1);

    for ii = 1:outPUT.nTraps+1
        plot(Temp2',res_time(ii,:),'--','LineWidth',3);
        legend_n{ii} = [ 'E_' num2str(ii,'%1d') ' = ' ...
                         num2str(outPUT.trap.Matrix(ii,4),'%3.1f') ' [eV]'];
    end

    legend(legend_n);

    plot([293,1293],[1,1])
    plot([293,1293],[1e-3,1e-3])
    plot([293,1293],[1e3,1e3])
    plot([383, 383],[1e-5,1e10])

    xlabel('Temperature [K]')
    ylabel('Residence Time [s]')
    xlim([300,1200])
    ylim([1e-5,1e10])
    set(gca,'YScale','log'); 
end

if isfield(inPUT.implant,'use_mod_template')
    outPUT.implant.use_mod_template = inPUT.implant.use_mod_template;
else
    outPUT.implant.use_mod_template = 0;
end

end

function [outPUT] = Organize_Output(inPUT,goodFit,varA)
 %
 %
      
% Store output vector
 outPUT = inPUT;

% Find best parameters
iFit = zeros(1,4);

[~,iFit(1)] = nanmax(goodFit.pr);
[~,iFit(2)] = nanmax(goodFit.sf);
[~,iFit(3)] = nanmax(goodFit.rt);

goodFit.tot = inPUT.fitP.w_pr*goodFit.pr + ...
              inPUT.fitP.w_sf*goodFit.sf + ...
              inPUT.fitP.w_rt*goodFit.rt;
      
[~,iFit(4)] = nanmax(goodFit.tot);

outPUT.nTraps = length(inPUT.trap.z2_damg);
outPUT.trap.Matrix = zeros(outPUT.nTraps,4);
outPUT.trap.Matrix(:,1) = inPUT.trap.z1_surf;
outPUT.trap.Matrix(:,2) = inPUT.trap.z2_damg;
outPUT.trap.Matrix(:,3) = inPUT.trap.z3_intr;
outPUT.trap.Matrix(:,4) = inPUT.trap.E_trap;
outPUT.trap.Matrix(inPUT.fitP.trapChng,inPUT.fitP.indxChng) = varA(iFit(4));

outPUT.trap.z1_surf = outPUT.trap.Matrix(:,1);
outPUT.trap.z2_damg = outPUT.trap.Matrix(:,2);
outPUT.trap.z3_intr = outPUT.trap.Matrix(:,3);
outPUT.trap.E_trap  = outPUT.trap.Matrix(:,4);

outPUT.fitP.prof  = goodFit.pr(  iFit(4));  
outPUT.fitP.surf  = goodFit.sf(  iFit(4));
outPUT.fitP.retn  = goodFit.rt(  iFit(4));
outPUT.fitP.tot   = goodFit.tot( iFit(4));

% pChng = inA.fitting(4)*100;                     % percent change

outPUT.plot.YorN = inPUT.plot.YorN;           % 1 = Yes, 0 = No

disp(['Best Fit for trap # ' ...
      num2str(inPUT.fitP.trapChng,'%3d') ...
             ' and index # ' ...
      num2str(inPUT.fitP.indxChng,'%3d')])
disp('prof  sflux  rtntn  total')
disp(iFit)

% LEGACY may not work correctly
% stringy = ['Prof, SurfFlux, Rentn = ' ...
%             num2str(iFit(1), '%1d') ', ' ...
%             num2str(iFit(2), '%1d') ', ' ...
%             num2str(iFit(3), '%1d') ': ' ...
%            'Choose ' ...
%             num2str(iFit(4), '%1d') ': ' ...
%            'Trap ' 
%             num2str(trapChng, '%1d') ': ' ...
%            'Var ' 
%             num2str(indxChng, '%1d') ': ' ...
%            'Chng ' 
%             num2str(pChng, '%4.3f') ' %'];

% disp(stringy) 

end

function [depth] = Import_Depth_Grid()
%
%   Find file '4_thermal.inp' with TMAP spatial layout and import into 
%   the structure depth. Both depth and thickness, '.X' and '.T'
%

global folder

fullName = [folder.Temp '4_thermal.inp'];

delimiter = {',','*','='};
formatSpec = '%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'Delimiter', delimiter,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

% find depth data between delx and tempd in .inp file
idx1 = find(strcmp(Var1, 'delx '),1,'first');
idx2 = find(strcmp(Var1,'tempd '),1,'first');

depth_t = str2double(Var1(idx1+1:idx2-2));
depth_t = [0;depth_t;0];

kk = length(depth_t);
depth_m = depth_t;
for i = 2:kk
    depth_m(i) = sum(depth_t(1:i));
end


depth.X = depth_m*1e6;      % convert meters to microns
depth.T = depth_t*1e6;

depth.segments = kk;

end

function [depth] = Create_Depth_Grid(inPUT)
%
%

% find the peak of the damage profile
space1 = 0:0.01:2;
gauss_1 = inPUT.trap.Profile_z2.gauss_1;
gauss_2 = inPUT.trap.Profile_z2.gauss_2;
Damage = gauss_1(1)*exp(-((gauss_1(2)-space1)/gauss_1(3)).^2) + ...
         gauss_2(1)*exp(-((gauss_2(2)-space1)/gauss_2(3)).^2);
[~,dmg_pk_] = max(Damage);
dmg_pk_ = space1(dmg_pk_)*1e-6;

start = 6.0e-10;                            % smallest thickness
srfPk = inPUT.Depth.dmg_depth;
damage1 = dmg_pk_*4;                        % 4 times the damage peak depth
thick = inPUT.Depth.tot_thick*0.5;          % 1/2 the total thickness (turn around)

ii_total = inPUT.Depth.n_six*6 - 2;

ii_0 = floor(ii_total* 4/40);                   % plasma implant zone
ii_1 = floor(ii_total* 4/40);                   % transition to srim
ii_2 = floor(ii_total*10/40);                   % srim damage zone
ii_3 = floor(ii_total*10/40);                   % transition to intrinsic
ii_4 = floor(ii_total* 9/40);                   % intrinsic zone to center
ii_6 = ii_total - (ii_0+ii_1+ii_2+ii_3+ii_4);

ii_0 = ii_0 + 1;                % include +1 for repeated value
ii_1 = ii_1 + 1;            
ii_2 = ii_2 + 1;
ii_3 = ii_3 + 1;
ii_4 = ii_4 + 1;
ii_6 = ii_6 + 1;

zoom0 = srfPk;
zoom1 = srfPk*(ii_1);
zoom2 = damage1;
zoom3 = zoom2*25;

nn_0 = logspace(log10(start),log10(zoom0),ii_0);      % implantation stress

nn_1 = stpspace(zoom0,zoom1,1.05,ii_1);               % transition
nn_2 = stpspace(zoom1,zoom2,1.10,ii_2);               % srim profile
nn_3 = stpspace(zoom2,zoom3,1.25,ii_3);               % transition
nn_4 = stpspace(zoom3,thick,1.30,ii_4);               % intrinsic to center
nn_6 = logspace(log10(start),log10(thick),ii_6);      % center to back surface

% Cut and paste the grid together to define thickness slices
nn_1 = [nn_0,nn_1(2:end),nn_2(2:end),nn_3(2:end),nn_4(2:end)];
mm_1 = nn_1;
for ii = 2:(ii_0+ii_1+ii_2+ii_3+ii_4-5)
    mm_1(ii) = nn_1(ii+1)-nn_1(ii);
end
mm_1 = mm_1(1:end-1);

mm_6 = nn_6;
for ii = 2:ii_6-1
    mm_6(ii) = nn_6(ii+1)-nn_6(ii);
end
mm_6 = mm_6(1:end-1);

depth_t = [0,mm_1, fliplr(mm_6),0]';

kk = length(depth_t);
depth_m = depth_t;
for i = 2:kk
    depth_m(i) = sum(depth_t(1:i));
end

depth_c = depth_m - depth_t/2;

depth.X = depth_m*1e6;      % convert meters to microns
depth.T = depth_t*1e6;
depth.C = depth_c*1e6;

depth.segments = kk;

end

function [spaceX] = stpspace(min_x,max_x,stepsize,Nsteps)
%
%
delta_x = max_x - min_x;
Nsteps = Nsteps-1;
sum_N = 0;
for ii = 1:Nsteps
    sum_N = stepsize^(ii-1) + sum_N;
end

delta_N = delta_x/sum_N;

spaceT = zeros(1,Nsteps);
for ii = 1:Nsteps
    spaceT(ii) = delta_N*stepsize^(ii-1);
end

spaceX = spaceT;
for ii = 2:Nsteps
    spaceX(ii) = sum(spaceT(1:ii));
end
spaceX = [0,spaceX] + min_x;
end

function [trapD0] = Calculate_traps(inA,trapD0)
% Calculations begin
%

% define the depth steps
depthX = trapD0.depth.X;
depthC = trapD0.depth.C;

% Here we assume that the lowest detrapping energy has a trap profile
% [~,iPk1] = min(abs(depthX-nWidth*inA.Depth.dmg_depth*1e6));     % Find depth 
% Plasma  = 1*exp(-((inA.Depth.dmg_depth*1e6/2   - depthC)/ (inA.Depth.dmg_depth*1e6/8)).^2);
% Plasma1 = 1*exp(-((inA.implantGauss.center*1e6 - depthC)/ (inA.Depth.dmg_depth*1e6/4)).^2);

% calculate plasma induced stress zone 1
exp_1 = inA.trap.Profile_z1.exp_1;
Plasma = 1*exp(-(exp_1)*depthC);

% calculate heavy ion damage profile zone 2
gauss_1 = inA.trap.Profile_z2.gauss_1;
gauss_2 = inA.trap.Profile_z2.gauss_2;

Damage = gauss_1(1)*exp(-((gauss_1(2)-depthC)/gauss_1(3)).^2) + ...
         gauss_2(1)*exp(-((gauss_2(2)-depthC)/gauss_2(3)).^2);

zero1 = zeros(length(Damage),1);
ones1 =  ones(length(Damage),1);

% store trap profiles for all traps
zeroMatrix = zeros(inA.nTraps+1,length(depthX));
trapC  = zeroMatrix;                    % trap concentration (at. fraction)
concF  = zeroMatrix;                    % concentration filled (fraction) at t=0
concE  = zeroMatrix;                    % concentration filled when trap turned off

for ii = 1:inA.nTraps
    % trap concentration
    surf_peak = inA.trap.Matrix(ii,1)*Plasma; 
    displcmnt = inA.trap.Matrix(ii,2)*Damage;
    intrinsic = inA.trap.Matrix(ii,3)*ones1;
    trapC(ii,:) = max(surf_peak,max(displcmnt,intrinsic));
end

trapC(:,  1) = 0;                       % Ensure endpoints are zero
trapC(:,end) = 0;

trapD0.trapC = trapC;
trapD0.concF = concF;
trapD0.concF = concF;
trapD0.concE = concE;

trapD0.trapC(inA.nTraps+1,:) = sum(trapD0.trapC(3:inA.nTraps,:),1);    % Here traps 3 to end (e.g. 6) are combined 

plotTraps = 0;

if plotTraps == 1
    figure(99)
    clf
    hold on
    for ii = 1:inA.nTraps
        stairs([0;depthX],[trapC(ii,:),0]*100,'LineWidth',2)
    end
    
    trap_total = sum(trapD0.trapC(1:inA.nTraps,:),1);
    stairs([0;depthX],[trap_total,0]*100,'LineWidth',3)
    
%     trap_total = sum(trapD0.trapC(2:inA.nTraps,:),1);
%     stairs([0;depthX],[trap_total,0]*100,'LineWidth',3)
%     
%     trap_total = sum(trapD0.trapC(3:inA.nTraps,:),1);
%     stairs([0;depthX],[trap_total,0]*100,'LineWidth',3)
%     
%     trap_total = sum(trapD0.trapC(4:inA.nTraps,:),1);
%     stairs([0;depthX],[trap_total,0]*100,'LineWidth',3)
%     
%     trap_total = sum(trapD0.trapC(5:inA.nTraps,:),1);
%     stairs([0;depthX],[trap_total,0]*100,'LineWidth',3)
    
    xlabel('Depth [\mum]')
    ylabel('Trap Concentration [at. %]')
    xlim([0 2.5])
end


end

function [inA,trapD] = Implant_phase(inA,trapD)
% Simulate the Implantation phase
%   SRIM provides Trap Profile used in Calculate_traps
%   Initially write and run TMAP file for 1 second with small time step
%   Reinitiallize TMAP with larger time step for duration of implant time
%   Run TMAP with temperature ramp down to room temp

% Delete all old files in the Updated folder and the Output folder
Delete_Old_Files();

trapD = Calculate_traps(inA,trapD);            % Setup trap concentrations

% Initial TMAP run with small time step at const. elevated Temp
inA.jRun = 1;                                  % jRun = 1 for small step
[inA,trapD] = get_Current_traps(inA,trapD);
inA.trap.Matrix(inA.nTraps+1,4) = Pseudo_Energy(inA,3);
Write_Input_File(inA,trapD);
Run_Updated_TMAP(inA,inA.iLeap,inA.iFit);   

% Continue TMAP run with large time step at const. elevated Temp
inA.jRun = 2;                                  % jRun = 2 for large step
[inA,trapD] = get_Current_traps(inA,trapD);
Write_Input_File(inA,trapD);
Run_Updated_TMAP(inA,inA.iLeap,inA.iFit);    

% Continue TMAP run with large time step with Temp decreasing to RT
inA.jRun = 3;                                  % jRun = 3 for ramp down
[inA,trapD] = get_Current_traps(inA,trapD);
Write_Input_File(inA,trapD);
Run_Updated_TMAP(inA,inA.iLeap,inA.iFit);

% Load completed Implantation of D and trap filling
inA.jRun = 4;
inA.adjustP = 1;
[inA,trapD] = get_Current_traps(inA,trapD);

end

function [gof_pr] = Plot_NRA_TMAP(sample,inA,trapD,mm)
%
%

nTraps = inA.nTraps;

depth_ = trapD.depth.X;
trapC  = trapD.trapC(1:nTraps,:)*100;              % trap concentration [at. %]
concF  = trapD.concF(1:nTraps,:).*trapC;           % filled concentration 
mob_0  = trapD.Crnt.m0/inA.eqn.N_w*100;            % Mobile conc.

trapT  = sum(trapC,1);                             % Total trap concentration
concT  = sum(concF,1);

trapC_0 = zeros(nTraps,1);
depthT = trapD.depth.T;
iEnd = length(trapD.depth.T);
for ij = 1:nTraps
   trapC_0(ij) = sum(trapC(ij,1:iEnd)/100.*depthT(1:iEnd)')/sum(depthT(1:iEnd));
end

T_total = (trapC*depthT)/100*(inA.eqn.N_w*1e-6);
D_total = (concF*depthT)/100*(inA.eqn.N_w*1e-6);

% Calculate goodness of fit
Depth_NRA = sample.NRA.Depth(2:end);        % Experimental depth
Dconc_NRA = sample.NRA.Dconc(1:end-1);

nFitss = 500;
linOUT = Compare_on_same_grid(Depth_NRA,Dconc_NRA,depth_,concT,nFitss);

% Depth_lin = linOUT.x_out;
Dconc_lin = linOUT.y_out_1;
DTMAP_lin = linOUT.y_out_2;

gof_pr = goodnessOfFit(Dconc_lin(2:end),DTMAP_lin(2:end),'NRMSE');

if inA.plot.YorN == 1
    hFig = figure(mm);
    set(hFig, 'Units', 'normalized', 'Position', [0.45 0.525 0.545 0.40])

    clf
    hold all
    
    pColor = zeros(nTraps,3);
    pColor(1,:) = [0 0 1];
    pColor(2,:) = [0.04 0.52 0.78];
    pColor(3,:) = [0 1 1];
    pColor(4,:) = [0.75 0.75 0];
    pColor(5,:) = [1 0 0];
    pColor(6,:) = [0.48 0.06 0.89];

    dark_green = [0 0.5 0];

    stairs(sample.NRA.Depth+1e-100,sample.NRA.Dconc,'-k','LineWidth',6);
    stairs(depth_,[concT(2:end),concT(end)],'-','Color',dark_green,'LineWidth',3);
    for ii = 1:nTraps
        stairs(depth_,[concF(ii,2:end),concF(ii,end)], '--','Color',pColor(ii,:),'LineWidth',2.25);
    end
    stairs(depth_,[trapT(2:end),trapT(end)],'--m','LineWidth',2.25);
    stairs(depth_,[mob_0(2:end),mob_0(end)]*100,'--g','LineWidth',2.25);

    stairs(sample.NRA.Depth+1e-100,sample.NRA.Dconc,'-k','LineWidth',6);
    stairs(depth_,[concT(2:end),concT(end)],'-','Color',dark_green,'LineWidth',3);

    xlabel('Depth [\mum]')
    ylabel('D Concentration [at %]')
    box on

    legend_n = cell(nTraps+2,1);
    legend_n{1} = [ 'NRA : ' num2str(sample.NRA.D_total,'%7.2e') ' [D/m^2] '];
    legend_n{2} = [ 'TMAP Total: '  num2str(sum(D_total)/sample.NRA.D_total*100,'%6.2f') ' %' ];
    for ii = 1:nTraps
        legend_n{ii+2} = [ 'E_' num2str(ii,'%1d') '= ' ...
                                num2str(inA.trap.Matrix(ii,4),'%4.2f') ' [eV]: ' ...
                                num2str(D_total(ii)/sample.NRA.D_total*100,'%4.2f') ' %' ];
    end

    legend(legend_n);
    xlim([.001 50])
%     xlim([.05 50])
    ylim([9.9e-6 ceil(max(sample.NRA.Dconc)*1.02)])

    % How to scale plot
    X_scale = 'log';
    % Y_scale = 'linear';
    Y_scale = 'log';

    set(gca,'XScale',X_scale);
    set(gca,'YScale',Y_scale); 

    if  inA.iLeap == inA.nLeap + 1
        title(['NRA NRMSE = ' num2str(gof_pr,'%9.7f')]);
    else
        Temp_switch = inA.Table1.TempEnd(inA.iLeap);
        res_time = 1/inA.eqn.jmpF*exp(inA.trap.Matrix(:,4)./(inA.eqn.k_B*Temp_switch'));
        title([ 'Depth Profile up to ' ...
                num2str(Temp_switch,'%5.1f') ...
                ' K : Res. Time (' ...
                num2str(inA.iLeap+0,'%2d') ...
                ' & ' ...
                num2str(inA.iLeap+2,'%2d') ...
                ')= '... 
                num2str(res_time(inA.iLeap+0),'%5.3e') ...
                ' & ' ...
                num2str(res_time(inA.iLeap+2),'%5.3e') ...
                ' s']);
    end

end

end

function [out] = Compare_on_same_grid(x_1,y_1,x_2,y_2,nn)
%
%   

% Define spatial grid
x_out = x_1(end)/(nn-1);
x_out = 0:x_out:x_1(end);

% Setup y_1 grid
y_out_1 = zeros(nn,1);
kLength = length(x_1);

for kk = 1:kLength
    kEnd = (kLength+1) - kk;
    for jj = 1:nn
        if  x_out(jj)   <=  x_1(kEnd)
            y_out_1(jj)  =  y_1(kEnd);
        end
    end
end

% Setup y_2 grid
y_out_2 = zeros(nn,1);
kLength = length(x_2);

for kk = 1:kLength
    kEnd = (kLength+1) - kk;
    for jj = 1:nn
        if  x_out(jj)   <=  x_2(kEnd)
            y_out_2(jj)  =  y_2(kEnd);
        end
    end
end

out.x_out = x_out;
out.y_out_1 = y_out_1;
out.y_out_2 = y_out_2;

end

function [inA,trapD] = TDS_phase(sample,inA,trapD,ii)
%
%

for kk = 1:inA.nLeap     
     inA.iLeap = kk;
     [inA,trapD] = TDS_segment(inA,trapD);
     Plot_NRA_TMAP(sample,inA,trapD,30+ii);
end

end

function [inA,trapD] = TDS_segment(inA,trapD)
% Simulate a TDS segment and perform PTTP transition
%

% Initial TMAP run with small time step at const. elevated Temp
inA.jRun = 1;                                  % jRun = 1 for small step
trapD.order.crnt = inA.order(inA.iLeap,:);
trapD = load_traps(inA,trapD);
Write_Input_File(inA,trapD);
Run_Updated_TMAP(inA,inA.iLeap,inA.iFit);

% Continue TMAP run with large time step at const. elevated Temp
inA.jRun = 2;                                  % jRun = 2 for large step
[inA,trapD] = get_Current_traps(inA,trapD);
Write_Input_File(inA,trapD);
Run_Updated_TMAP(inA,inA.iLeap,inA.iFit);

% Load current traps and adjust pseudo trap
inA.jRun = 3;
inA.adjustP = 1;
[inA,trapD] = get_Current_traps(inA,trapD);

end

function [E_pseudo] = Pseudo_Energy(inA,jBgn)
%  
% E_pseudo = -log(  sum( inA.trap.Matrix(jBgn:end,2)...
%                 .*exp(-inA.trap.Matrix(jBgn:end,4)))...
%                 / sum( inA.trap.Matrix(jBgn:end,2)));
    
max1 = max(inA.trap.Matrix(jBgn:end,2) , inA.trap.Matrix(jBgn:end,3));
exp1 = exp(-inA.trap.Matrix(jBgn:end,4)/(inA.eqn.k_B*inA.Table1.TempBgn(inA.iLeap)) );            
E_pseudo = -log(  sum( max1.*exp1 )/ sum( max1 ) )* ...
               (inA.eqn.k_B*inA.Table1.TempBgn(inA.iLeap));

% disp(E_pseudo)           

end

function [inA,trapD] = Pseudo_Trap(inA,trapD)
% Decide how to partition traps
% Use following scheme
% traps = [1,2,p1] where p1 = 3:n for Segment 0 and 1
%       = [p2,2,3] where p2 = 4:n for Segment 2
%       = [4,p3,3] where p2 = 5:n for Segment 3
%       = [4,5,p4] where p2 = 6:n for Segment 4
%       etc...

pseuD = trapD.concF(inA.nTraps+1,:).*trapD.trapC(inA.nTraps+1,:); % at. fraction of D in pseudo trap 

if inA.jRun == 4            % End of Implantation phase
    jPseudo = inA.nTraps-2;
    for jj = 1:jPseudo
       ii = inA.nTraps+1 - jj;
       pseuD_ = max(pseuD-trapD.trapC(ii,:),0);       % completely fill highest trap
       trapD.concF(ii,:) = (pseuD - pseuD_)./trapD.trapC(ii,:);
       pseuD = pseuD_;                                
    end
    % Pseudo trap conc stays the same since 1st TDS segement same as Impantation order
    jPseudo = inA.nTraps+1 - jPseudo;
    inA.trap.Matrix(inA.nTraps+1,4) = Pseudo_Energy(inA,jPseudo);
end

if inA.jRun == 3            % End of TDS Segment
    jPseudo = inA.nTraps - inA.iLeap - 1;
    for jj = 1:jPseudo
        ii = inA.nTraps+1 - jj;
        pseuD_ = max(pseuD-trapD.trapC(ii,:),0);       % completely fill highest trap
        trapD.concF(ii,:) = (pseuD - pseuD_)./trapD.trapC(ii,:);
        pseuD = pseuD_;
    end
    jPseudo = inA.nTraps+1 - jPseudo;
    trapD.trapC(inA.nTraps+1,:) = trapD.trapC(inA.nTraps+1,:) - trapD.trapC(jPseudo,:);
    if inA.iLeap ~= inA.nLeap
        inA.trap.Matrix(inA.nTraps+1,4) = Pseudo_Energy(inA,jPseudo+1);
    end
end

% disp(inA.trap.Matrix(inA.nTraps+1,4))

trapD.concF(isnan(trapD.concF)) = 0;
trapD.concF(isinf(trapD.concF)) = 0;

end

function [] = Delete_Old_Files()
global folder

fID1 = fopen([folder.TMAP 'run_Delete_Old_Files.bat'],'w');
changeDir = ['cd ' folder.Updt];
fprintf(fID1,'%s\n',changeDir);    
fprintf(fID1,'%s\n',['del '  '*.inp']);
changeDir = ['cd ' folder.Work];
fprintf(fID1,'%s\n',changeDir);  
fprintf(fID1,'%s\n',['del ' folder.tmapName '*']);
fclose(fID1);

[~,~] = system([folder.TMAP 'run_Delete_Old_Files']);
end

function [] = Run_Updated_TMAP(inA,kk,ii)
%
%
%
global folder

nLeaps = length(folder.tmapNames(:,1));
kkStep = mod(kk,nLeaps);                   % here kkStep=0 for implantation

fID1 = fopen([folder.TMAP 'run_Updated_TMAP.bat'],'w');
changeDir = ['cd ' folder.TMAP];
runTMAP   = ['START "' ...
             num2str(kkStep,'%02d') ' of ' ...
             num2str(nLeaps-1,'%02d')...
             ' cycles | '...
             num2str(ii,'%02d') ' of ' ...
             num2str(length(folder.tmapNames(1,:)),'%02d')...
             ' variations | trap # '...      
             num2str(inA.fitP.trapChng,'%02d')...
             ' and variable # '...
             num2str(inA.fitP.indxChng,'%02d')...
             ' by '...
             num2str(inA.fitP.fracChng*100,'%2.2f')...
             ' %% '...
             '" /wait cmd /c ' ...
             't7 ' folder.tmapNames{kk,ii}];
fprintf(fID1,'%s\n',changeDir);
fprintf(fID1,'%s\n',runTMAP);
fclose(fID1);

fID2 = fopen([folder.TMAP 'run_Updated_Ex_surfF.bat'],'w');
echo1 = ['(echo ' folder.tmapNames{kk,ii} '.plt'];
echo2 = [' echo ' folder.tmapNames{kk,ii} '_sf' ]; 
echo3 = ['echo n) | ' folder.TMAP 'extract'];
fprintf(fID2,'%s\n',changeDir);
fprintf(fID2,'%s\n',  echo1);
fprintf(fID2,'%s\n',  echo2);
fprintf(fID2,'%s\n',' echo d');
fprintf(fID2,'%s\n',' echo sf');
fprintf(fID2,'%s\n',' echo l');
fprintf(fID2,'%s\n',' echo 1');
fprintf(fID2,'%s\n',' echo 1');
fprintf(fID2,'%s\n',  echo3);
fclose(fID2);

fID3 = fopen([folder.TMAP 'run_Updated_Ex_tempK.bat'],'w');
echo1 = ['(echo ' folder.tmapNames{kk,ii} '.plt'];
echo2 = [' echo ' folder.tmapNames{kk,ii} '_tK' ]; 
echo3 = ['echo n) | ' folder.TMAP 'extract'];
fprintf(fID3,'%s\n',changeDir);
fprintf(fID3,'%s\n',  echo1);
fprintf(fID3,'%s\n',  echo2);
fprintf(fID3,'%s\n',' echo d');
fprintf(fID3,'%s\n',' echo st');
fprintf(fID3,'%s\n',' echo l');
fprintf(fID3,'%s\n',' echo 1');
fprintf(fID3,'%s\n',' echo 1');
fprintf(fID3,'%s\n',  echo3);
fclose(fID3);

fID4 = fopen([folder.TMAP 'run_Updated_Ex_trap1.bat'],'w');
echo1 = ['(echo ' folder.tmapNames{kk,ii} '.plt'];
echo2 = [' echo ' folder.tmapNames{kk,ii} '_t1' ]; 
echo3 = ['echo n) | ' folder.TMAP 'extract'];
fprintf(fID4,'%s\n',changeDir);
fprintf(fID4,'%s\n',  echo1);
fprintf(fID4,'%s\n',  echo2);
fprintf(fID4,'%s\n',' echo d');
fprintf(fID4,'%s\n',' echo ti');
fprintf(fID4,'%s\n',' echo 1');
fprintf(fID4,'%s\n',  echo3);
fclose(fID4);

fID5 = fopen([folder.TMAP 'run_Updated_Ex_trap2.bat'],'w');
echo1 = ['(echo ' folder.tmapNames{kk,ii} '.plt'];
echo2 = [' echo ' folder.tmapNames{kk,ii} '_t2' ]; 
echo3 = ['echo n) | ' folder.TMAP 'extract'];
fprintf(fID5,'%s\n',changeDir);
fprintf(fID5,'%s\n',  echo1);
fprintf(fID5,'%s\n',  echo2);
fprintf(fID5,'%s\n',' echo d');
fprintf(fID5,'%s\n',' echo ti');
fprintf(fID5,'%s\n',' echo 2');
fprintf(fID5,'%s\n',  echo3);
fclose(fID5);

fID6 = fopen([folder.TMAP 'run_Updated_Ex_trap3.bat'],'w');
echo1 = ['(echo ' folder.tmapNames{kk,ii} '.plt'];
echo2 = [' echo ' folder.tmapNames{kk,ii} '_t3' ]; 
echo3 = ['echo n) | ' folder.TMAP 'extract'];
fprintf(fID6,'%s\n',changeDir);
fprintf(fID6,'%s\n',  echo1);
fprintf(fID6,'%s\n',  echo2);
fprintf(fID6,'%s\n',' echo d');
fprintf(fID6,'%s\n',' echo ti');
fprintf(fID6,'%s\n',' echo 3');
fprintf(fID6,'%s\n',  echo3);
fclose(fID6);

fID7 = fopen([folder.TMAP 'run_Updated_Ex_mobl0.bat'],'w');
echo1 = ['(echo ' folder.tmapNames{kk,ii} '.plt'];
echo2 = [' echo ' folder.tmapNames{kk,ii} '_m0' ]; 
echo3 = ['echo n) | ' folder.TMAP 'extract'];
fprintf(fID7,'%s\n',changeDir);
fprintf(fID7,'%s\n',  echo1);
fprintf(fID7,'%s\n',  echo2);
fprintf(fID7,'%s\n',' echo d');
fprintf(fID7,'%s\n',' echo mi');
fprintf(fID6,'%s\n',' echo 1');
fprintf(fID6,'%s\n',' echo 1');
fprintf(fID7,'%s\n',  echo3);
fclose(fID7);

fID8 = fopen([folder.TMAP 'run_Updated_Cleanup.bat'],'w');
fprintf(fID8,'%s\n',changeDir);
moveTMAP1  = ['move ' folder.tmapNames{kk,ii} '.inp ' folder.Work];  
moveTMAP2  = ['move ' folder.tmapNames{kk,ii} '.out ' folder.Work];       
moveTMAP3  = ['move ' folder.tmapNames{kk,ii} '.plt ' folder.Work];      
moveTMAP4  = ['move ' folder.tmapNames{kk,ii} '_sf '  folder.Work]; 
moveTMAP5  = ['move ' folder.tmapNames{kk,ii} '_tK '  folder.Work];          
moveTMAP6  = ['move ' folder.tmapNames{kk,ii} '_t1 '  folder.Work];      
moveTMAP7  = ['move ' folder.tmapNames{kk,ii} '_t2 '  folder.Work];      
moveTMAP8  = ['move ' folder.tmapNames{kk,ii} '_t3 '  folder.Work];      
moveTMAP9  = ['move ' folder.tmapNames{kk,ii} '_m0 '  folder.Work];       
fprintf(fID8,'%s\n',moveTMAP1);
fprintf(fID8,'%s\n',moveTMAP2);
fprintf(fID8,'%s\n',moveTMAP3);
fprintf(fID8,'%s\n',moveTMAP4);
fprintf(fID8,'%s\n',moveTMAP5);
fprintf(fID8,'%s\n',moveTMAP6);
fprintf(fID8,'%s\n',moveTMAP7);
fprintf(fID8,'%s\n',moveTMAP8);
fprintf(fID8,'%s\n',moveTMAP9);
fclose(fID8);

[~,~] = system([folder.TMAP 'run_Updated_TMAP']);

[~,~] = system([folder.TMAP 'run_Updated_Ex_surfF']);
[~,~] = system([folder.TMAP 'run_Updated_Ex_tempK']);

[~,~] = system([folder.TMAP 'run_Updated_Ex_trap1']);
[~,~] = system([folder.TMAP 'run_Updated_Ex_trap2']);
[~,~] = system([folder.TMAP 'run_Updated_Ex_trap3']);
[~,~] = system([folder.TMAP 'run_Updated_Ex_mobl0']);

[~,~] = system([folder.TMAP 'run_Updated_Cleanup']);

% gof.sf = zeros(1,iNames);
% gof.rt = zeros(1,iNames);
% 
% for ii = 1:iNames
%     if ~isempty(folder.tmapNames{ii})
%         kk = 40+ii;
%         [outPlot] = Plot_TMAP_surf_flux_00(folder,folder.tmapNames{ii},sample,kk);
%         gof.sf(ii) = outPlot.normRMSE;
%         gof.rt(ii) = outPlot.gof.rt;
%     end
% end

end

function [out] = Plot_TDS_TMAP(sample,inA,trapD,mm)
%
%
%
global folder

Temp_K = sample.TDS.tempK;              % experimental TDS data
D_flux = sample.TDS.Dflux;
D_rtn  = sample.TDS.D_total;

% Temp_K = sample.TDS.temptr;              % experimental TDS data
% D_flux = sample.TDS.D_flux;
% D_rtn  = sample.TDS.Dreten;

time_s  = [];
surflux = [];
tempK   = [];

trap1=[];
trap2=[];
trap3=[];
mobl0=[];

cycleLen = [];

for kk = 1:inA.iLeap
    fullName    = [folder.Work folder.tmapNames{kk,inA.iFit}]; % finds directory
    fullName_sf = [fullName '_sf']; % appends _sf
    fullName_tK = [fullName '_tK'];
    fullName_t1 = [fullName '_t1']; % appends _t1
    fullName_t2 = [fullName '_t2'];
    fullName_t3 = [fullName '_t3'];
    fullName_m0 = [fullName '_m0'];
    [surf]      = get_Extracted_Data(fullName_sf);
    [temp]      = get_Extracted_Data(fullName_tK);
    [tinv1]     = get_Extracted_Data(fullName_t1);
    [tinv2]     = get_Extracted_Data(fullName_t2);
    [tinv3]     = get_Extracted_Data(fullName_t3);
    [minv0]     = get_Extracted_Data(fullName_m0);
    time_       = surf.time_s(1:end);
    cycleLen    = [cycleLen;length(time_)-1];
    time_s      = [time_s;inA.Table1.timeCont(kk) + time_(2:end)];
    surflux     = [surflux;surf.data_(2:end)];
    tempK       = [  tempK;temp.data_(2:end)];
    trap1       = [trap1;tinv1.data_(2:end)];
    trap2       = [trap2;tinv2.data_(2:end)];
    trap3       = [trap3;tinv3.data_(2:end)];  
    mobl0       = [mobl0;minv0.data_(2:end)];
end

Dflux = -surflux;
deltaTime1 = time_s(3)-time_s(2);                                           % TMAP delta time
rampRate = (tempK(end-2)-tempK(1))/(time_s(end-2)-time_s(1));
deltaTime2 = (Temp_K(end-2) - Temp_K(1))/rampRate / length(Temp_K);         % TDS delta time

out.TMAP.time_s = time_s;
out.TMAP.Dflux = Dflux;
out.TMAP.tempK = tempK;
out.TMAP.trap1 = trap1;
out.TMAP.trap2 = trap2;
out.TMAP.trap3 = trap3;
out.TMAP.mobl0 = mobl0;

% figure(2);plot(tempK,Dflux)

Dflux2  = interp1(tempK,Dflux,Temp_K,'linear');
D_Diff  = D_flux(2:end-8) - Dflux2(2:end-8);                            % 2:end-8 is a bandaid

D_total = sum(max(0,Dflux))*deltaTime1;
E_total = sum(abs(D_Diff))*deltaTime2;

gof_sf = goodnessOfFit(D_flux(2:end-8),Dflux2(2:end-8),'NRMSE');

out.D_total = D_total;
out.gof_sf = gof_sf;
out.gof_rt = 1 - E_total/D_rtn;

trapC = trapD.trapC*100;             % trap concentration [at. %]
concF = trapD.concF.*trapC;          % filled concentration 
concE = trapD.concE*100;

depthT = trapD.depth.T;
D_totalF = (concF*depthT)/100*(inA.eqn.N_w*1e-6);
D_totalE = (concE*depthT)/100*(inA.eqn.N_w*1e-6);

D_rem = D_totalE./D_totalF*100;     % remaining trapped D

if inA.plot.YorN == 1
    hFig = figure(mm);
    set(hFig, 'Units', 'normalized', 'Position', [0.45 0.04 0.545 0.40])

    scale_ = 1e17;
    yMax = ceil(max([max(D_flux) max(Dflux2)])/scale_*1.015);
    clf
    hold all
    
    dark_green  = [0 0.5 0];
    dark_orange = [1 0.5 0];
    dark_blue   = [0.08 0.17 0.55];

    plot(Temp_K,abs(D_flux)/scale_,'k-', 'LineWidth',6)
    plot(Temp_K,    Dflux2 /scale_,'-','Color',dark_green,'LineWidth',3)
    plot(Temp_K(2:end-8),abs(D_Diff)/scale_,'r--','LineWidth',2.25)
    for ii = 1:inA.iLeap-1
        xTemp = [inA.Table1.TempEnd(ii) inA.Table1.TempEnd(ii)];
        yMax2 = [0 yMax];
        plot(xTemp,yMax2,'-.','Color',dark_orange,'LineWidth',2.50)
    end   

    title(['TDS NRMSE = ' num2str(gof_sf,'%6.5f') ])
    xlim([350 1050]);
    ylim([0 yMax])
    xlabel('Thermal Desorption Temperature [K]');
    ylabel('Surface Flux [10^{17} D/m^2/s]');
    box on

    legend_n = cell(inA.iLeap-1+3,1);
    legend_n{1} = [       'TDS = ' num2str(        D_rtn    ,'%4.3e')     ];
    legend_n{2} = [      'TMAP = ' num2str(D_total/D_rtn*100,'%4.2f') ' %'];
    legend_n{3} = ['|Residual| = ' num2str(E_total/D_rtn*100,'%4.2f') ' %'];
    for ii = 1:inA.iLeap-1
        legend_n{ii+3} = [ num2str(ii,'%1d') ' \rightarrow ' num2str(ii+3,'%1d') ' : '...
                           num2str(D_rem(ii),'%3.2f') ' , '  num2str(D_rem(ii+3),'%3.2f') ' % '];
    end
    
    legend(legend_n, 'Location','NorthEast')
    %[~,i_yMax] = max(D_flux);

    % if Temp_K(i_yMax) < ((1150-350)/2 + 350)
    %     legend(legend_n, 'Location','NorthEast')
    % else
    %     legend(legend_n, 'Location','NorthWest')
    % end

end

ShowTrapInv = 0 ;   % Setup only works for 6 traps as written

if ShowTrapInv == 1
% if inA.plot.YorN == 1
    scale_ = 1e20;
    
    hFig = figure(12);
    set(hFig, 'Units', 'normalized', 'Position', [1.003 0.01 0.4 0.90])
    clf
    hold on
    
    % Define the color map 
    cmap = zeros(8,3);
    cmap(1,:) = [0 0 1];            % blue
    cmap(2,:) = [0.04 0.52 0.78];   % turquoise
    cmap(3,:) = [0 1 1];            % cyan
    cmap(4,:) = [0.75 0.75 0];      % gold
    cmap(5,:) = [1 0 0];            % red
    cmap(6,:) = [0.48 0.06 0.89];   % purple
    cmap(7,:) = [0 0.5 0];          % dark green
    cmap(8,:) = [1 0.5 0];          % burnt orange
    
    if inA.nTraps == 6
    TransT = inA.Table1.TempBgn(2:end-1);
    [~,iTrans1] = min(abs(out.TMAP.tempK - TransT(1)));
    [~,iTrans2] = min(abs(out.TMAP.tempK - TransT(2)));
    [~,iTrans3] = min(abs(out.TMAP.tempK - TransT(3)));

    nTemp = length(out.TMAP.tempK);

    jTrap1 = (1:iTrans1-1);
    jTrap2 = (1:iTrans2-1);
    jTrap3 = (iTrans1+1:iTrans3-1);
    jTrap4 = (iTrans2+1:nTemp);
    jTrap5 = (iTrans3+1:nTemp);
    jTrap6 = (iTrans3+1:nTemp);

    jPseudo1 = (1:iTrans1);
    jPseudo2 = (iTrans1+1:iTrans2);
    jPseudo3 = (iTrans2+1:iTrans3);

    total_Inv = (out.TMAP.trap1+out.TMAP.trap2+out.TMAP.trap3);

    hFig = figure(12);
    set(hFig, 'Units', 'normalized', 'Position', [1.003 0.01 0.4 0.90])
    clf
    hold on
    
    plot(out.TMAP.tempK,total_Inv/scale_,'k-','LineWidth',6)
    plot(out.TMAP.tempK(jTrap1),out.TMAP.trap1(jTrap1)/scale_,'--','Color',cmap(1,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap2),out.TMAP.trap2(jTrap2)/scale_,'--','Color',cmap(2,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap3),out.TMAP.trap3(jTrap3)/scale_,'--','Color',cmap(3,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap4),out.TMAP.trap1(jTrap4)/scale_,'--','Color',cmap(4,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap5),out.TMAP.trap2(jTrap5)/scale_,'--','Color',cmap(5,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap6),out.TMAP.trap3(jTrap6)/scale_,'--','Color',cmap(6,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jPseudo1),out.TMAP.trap3(jPseudo1)/scale_,'--','Color',cmap(7,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jPseudo2),out.TMAP.trap1(jPseudo2)/scale_,'--','Color',cmap(7,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jPseudo3),out.TMAP.trap2(jPseudo3)/scale_,'--','Color',cmap(7,:),'LineWidth',2.25)
    maxInv = ceil(max(total_Inv/scale_*2))/2;
    plot([TransT(1),TransT(1)],[0, maxInv],'-.','Color',cmap(8,:),'LineWidth',2.25)
    plot([TransT(2),TransT(2)],[0, maxInv],'-.','Color',cmap(8,:),'LineWidth',2.25)
    plot([TransT(3),TransT(3)],[0, maxInv],'-.','Color',cmap(8,:),'LineWidth',2.25)

    legend('Trap total','Trap 1', 'Trap 2', 'Trap 3', 'Trap 4', 'Trap 5', 'Trap 6', 'Pseudo Trap')
    xlim([350,1050])
    ylim([0,maxInv])
    box on
    xlabel('Thermal Desorption Temperature [K]')
    ylabel('Trap Inventory: I_k = \int C_k dx [10^{20} D/m^2]')
 
    elseif inA.nTraps == 4
    TransT = inA.Table1.TempBgn(2:end-1);
    [~,iTrans1] = min(abs(out.TMAP.tempK - TransT(1)));

    nTemp = length(out.TMAP.tempK);

    jTrap1 = (1:iTrans1-1);
    jTrap2 = (1:nTemp);
    jTrap3 = (iTrans1+1:nTemp);
    jTrap4 = (iTrans1+1:nTemp);

    jPseudo1 = (1:iTrans1);
    
    total_Inv = (out.TMAP.trap1+out.TMAP.trap2+out.TMAP.trap3);

    plot(out.TMAP.tempK,total_Inv/scale_,'k-','LineWidth',6)
    plot(out.TMAP.tempK(jTrap1),out.TMAP.trap1(jTrap1)/scale_,'--','Color',cmap(1,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap2),out.TMAP.trap2(jTrap2)/scale_,'--','Color',cmap(2,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap3),out.TMAP.trap3(jTrap3)/scale_,'--','Color',cmap(3,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jTrap4),out.TMAP.trap1(jTrap4)/scale_,'--','Color',cmap(4,:),'LineWidth',2.25)
    plot(out.TMAP.tempK(jPseudo1),out.TMAP.trap3(jPseudo1)/scale_,'--','Color',cmap(5,:),'LineWidth',2.25)
    maxInv = ceil(max(total_Inv/scale_*2))/2;
    plot([TransT(1),TransT(1)],[0, maxInv],'-.','Color',cmap(8,:),'LineWidth',2.25)
    legend('Trap total','Trap 1', 'Trap 2', 'Trap 3', 'Trap 4', 'Pseudo Trap')

    end
    xlim([350,1050])
    ylim([0,maxInv])
    box on
    xlabel('Thermal Desorption Temperature [K]')
    ylabel('Trap Inventory: I_k = \int C_k dx [10^{20} D/m^2]')
    
 plot2 = 1;
 if plot2 == 1
    scale_ = 1e17;
    hFig = figure(13);
    set(hFig, 'Units', 'normalized', 'Position', [1.403 0.01 0.4 0.90])
    clf
    hold on
    
    if inA.nTraps == 6

%     dTrapT = diff(total_Inv)./diff(out.TMAP.tempK)*0.49;
    dTrap1 = diff(out.TMAP.trap1(jTrap1))./diff(out.TMAP.tempK(jTrap1))*0.49;
    dTrap2 = diff(out.TMAP.trap2(jTrap2))./diff(out.TMAP.tempK(jTrap2))*0.49;
    dTrap3 = diff(out.TMAP.trap3(jTrap3))./diff(out.TMAP.tempK(jTrap3))*0.49;
    dTrap4 = diff(out.TMAP.trap1(jTrap4))./diff(out.TMAP.tempK(jTrap4))*0.49;
    dTrap5 = diff(out.TMAP.trap2(jTrap5))./diff(out.TMAP.tempK(jTrap5))*0.49;
    dTrap6 = diff(out.TMAP.trap3(jTrap6))./diff(out.TMAP.tempK(jTrap6))*0.49;
    dPsdo1 = diff(out.TMAP.trap3(jPseudo1))./diff(out.TMAP.tempK(jPseudo1))*0.49;
    dPsdo2 = diff(out.TMAP.trap1(jPseudo2))./diff(out.TMAP.tempK(jPseudo2))*0.49;
    dPsdo3 = diff(out.TMAP.trap2(jPseudo3))./diff(out.TMAP.tempK(jPseudo3))*0.49;
    
    jTrap1 = jTrap1(1):jTrap1(end-1);
    jTrap2 = jTrap2(1):jTrap2(end-1);
    jTrap3 = jTrap3(1):jTrap3(end-1);
    jTrap4 = jTrap4(1):jTrap4(end-1);
    jTrap5 = jTrap5(1):jTrap5(end-1);
    jTrap6 = jTrap6(1):jTrap6(end-1);
    jPseudo1 = jPseudo1(1):jPseudo1(end-1);
    jPseudo2 = jPseudo2(1):jPseudo2(end-1);
    jPseudo3 = jPseudo3(1):jPseudo3(end-1);
    
%     plot(out.TMAP.tempK(1:end-1),-dTrapT/scale_,'k-','LineWidth',2)
    plot(out.TMAP.tempK(jTrap1),-dTrap1/scale_,'s--','Color',cmap(1,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap2),-dTrap2/scale_,'s--','Color',cmap(2,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap3),-dTrap3/scale_,'s--','Color',cmap(3,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap4),-dTrap4/scale_,'s--','Color',cmap(4,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap5),-dTrap5/scale_,'s--','Color',cmap(5,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap6),-dTrap6/scale_,'s--','Color',cmap(6,:),'LineWidth',2)
    plot(out.TMAP.tempK(jPseudo1),-dPsdo1/scale_,'s--','Color',cmap(7,:),'LineWidth',2)
    plot(out.TMAP.tempK(jPseudo2),-dPsdo2/scale_,'s--','Color',cmap(7,:),'LineWidth',2)
    plot(out.TMAP.tempK(jPseudo3),-dPsdo3/scale_,'s--','Color',cmap(7,:),'LineWidth',2)
    maxOUT = max([max(-dTrap1),max(-dTrap2),max(-dTrap3),...
                  max(-dTrap4),max(-dTrap5),max(-dTrap6),...
                  max(-dPsdo1),max(-dPsdo2),max(-dPsdo3)] );
    minOUT = max([max( dTrap1),max( dTrap2),max( dTrap3),...
                  max( dTrap4),max( dTrap5),max( dTrap6),...
                  max( dPsdo1),max( dPsdo2),max( dPsdo3)] );
    maxPseudo = ceil(maxOUT/scale_/4*10)*4/10;
    minPseudo = ceil(minOUT/scale_/4*10)*4/10;
%     maxPseudo = 9;
%     minPseudo = 1.5;
    plot([TransT(1),TransT(1)],[-minPseudo, maxPseudo],'-.','Color',cmap(8,:),'LineWidth',2.25)
    plot([TransT(2),TransT(2)],[-minPseudo, maxPseudo],'-.','Color',cmap(8,:),'LineWidth',2.25)
    plot([TransT(3),TransT(3)],[-minPseudo, maxPseudo],'-.','Color',cmap(8,:),'LineWidth',2.25)
    legend('Trap 1', 'Trap 2', 'Trap 3', 'Trap 4', 'Trap 5', 'Trap 6', 'Pseudo Trap')
    
    elseif inA.nTraps == 4

%     dTrapT = diff(total_Inv)./diff(out.TMAP.tempK)*0.49;
    dTrap1 = diff(out.TMAP.trap1(jTrap1))./diff(out.TMAP.tempK(jTrap1))*0.49;
    dTrap2 = diff(out.TMAP.trap2(jTrap2))./diff(out.TMAP.tempK(jTrap2))*0.49;
    dTrap3 = diff(out.TMAP.trap3(jTrap3))./diff(out.TMAP.tempK(jTrap3))*0.49;
    dTrap4 = diff(out.TMAP.trap1(jTrap4))./diff(out.TMAP.tempK(jTrap4))*0.49;
    dPsdo1 = diff(out.TMAP.trap3(jPseudo1))./diff(out.TMAP.tempK(jPseudo1))*0.49;
    
    jTrap1 = jTrap1(1):jTrap1(end-1);
    jTrap2 = jTrap2(1):jTrap2(end-1);
    jTrap3 = jTrap3(1):jTrap3(end-1);
    jTrap4 = jTrap4(1):jTrap4(end-1);
    jPseudo1 = jPseudo1(1):jPseudo1(end-1);
    
%     plot(out.TMAP.tempK(1:end-1),-dTrapT/scale_,'k-','LineWidth',6)
    plot(out.TMAP.tempK(jTrap1),-dTrap1/scale_,'s--','Color',cmap(1,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap2),-dTrap2/scale_,'s--','Color',cmap(2,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap3),-dTrap3/scale_,'s--','Color',cmap(3,:),'LineWidth',2)
    plot(out.TMAP.tempK(jTrap4),-dTrap4/scale_,'s--','Color',cmap(4,:),'LineWidth',2)
    plot(out.TMAP.tempK(jPseudo1),-dPsdo1/scale_,'s--','Color',cmap(7,:),'LineWidth',2)
    maxOUT = max([max(-dTrap1),max(-dTrap2),max(-dTrap3),max(-dTrap4),max(-dPsdo1)] );
    minOUT = max([max( dTrap1),max( dTrap2),max( dTrap3),max( dTrap4),max( dPsdo1)] );
    maxPseudo = ceil(maxOUT/scale_/4)*4;
    minPseudo = ceil(minOUT/scale_/4)*4;
%     maxPseudo = 9;
%     minPseudo = 1.5;
    plot([TransT(1),TransT(1)],[-minPseudo, maxPseudo],'-.','Color',cmap(8,:),'LineWidth',2.25)
    legend('Trap 1', 'Trap 2', 'Trap 3', 'Trap 4', 'Pseudo Trap')

    end
    
    plot([0,1200],[0,0],'k--')
    xlim([350,1050])
    ylim([-minPseudo,maxPseudo])
    box on
    xlabel('Thermal Desorption Temperature [K]')
    ylabel('Trap Release: d/dt (I_k) [10^{17} D/m^2/s]')
  
 end

end

end

function [extD] = get_Extracted_Data(fullName)
%
%
%

formatSpec = '%13f%f%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                    'Delimiter', '', 'WhiteSpace',...
                    '', 'EmptyValue' ,NaN,...
                    'ReturnOnError', false);
fclose(fileID);

% Allocate imported array to column variable names
extD.time_s  = dataArray{:, 1};      % time in seconds
extD.data_   = dataArray{:, 2};      % surface flux [D/m^2/s]
                                     % or trap inventory [#]

end

%% Extracting trap concentrations from TMAP output files

function [inA,trapD] = get_Current_traps(inA,trapD)
% Update the trap concentrations
% Determine if in implant phase iLeap == inA.nLeap + 1
%   implant jRun == 1 -> arbDconc for mobile conc., traps point to 1,2,and 7
%   implant jRun == 2 -> mobile and traps cont. from jRun = 1 output file
%                        points to same file M##_M_iLeap_##.out 
%   implant jRun == 3 -> mobile and traps cont. from jRun = 2 output file
% Separate traps from 1,2, and 7 to traps 1:6 and remove 7 altogether
%   gets data directly using parse since load would not account for trap 7 correctly
%
% Determine if in TDS phase
%   TDS iLeap == 1 -> reset mobile to arbDconc and grab traps from M##_M_iLeap_##.out
%   TDS iLeap == 2 -> mobile and traps cont. from iLeap = 1 output file
%                     point to previous file M##_M_1_##.out
%   TDS iLeap == 3 to end continue with previous ouput file
%
% Load the traps at the end
%   

flushDtrapped = 1;      % flush out D still in trap during transition and include in mobile conc? 
trapD.order.crnt = inA.order(inA.iLeap,:);

trapD.order.load = [inA.iLeap,inA.iFit];
trapD.order.last = trapD.order.crnt;             % Something wrong here 
trapD = load_traps(inA,trapD);

if inA.adjustP == 1
   [inA,trapD] = Pseudo_Trap(inA,trapD);
   inA.adjustP = 0;
   if inA.iLeap ~= inA.nLeap + 1
       % trapD.order.last is not pointing to correct order ... 
       if trapD.order.last(1) ~= trapD.order.crnt(1)
          if flushDtrapped == 1
             Dflush      = trapD.trapC(trapD.order.last(1),:).*trapD.concF(trapD.order.last(1),:);
             trapD.Crnt.m0 = trapD.Crnt.m0 + Dflush*inA.eqn.N_w;
          end
       end
       if trapD.order.last(2) ~= trapD.order.crnt(2)
          if flushDtrapped == 1
             Dflush      = trapD.trapC(trapD.order.last(2),:).*trapD.concF(trapD.order.last(2),:);
             trapD.Crnt.m0 = trapD.Crnt.m0 + Dflush*inA.eqn.N_w;
          end
       end    
       if trapD.order.last(3) ~= trapD.order.crnt(3)
          if flushDtrapped == 1            
             Dflush      = trapD.trapC(trapD.order.last(3),:).*trapD.concF(trapD.order.last(3),:);
             trapD.Crnt.m0 = trapD.Crnt.m0 + Dflush*inA.eqn.N_w;
          end
       end
   end
end


end

function [trapD] = load_traps(inA,trapD)
%
%   orderFrom = order of traps from output file
%   orderInto  = order of traps to be inserted into input file

% loadFrom = [inA.iLeap,inA.iFit];      % points to file to be parsed rom
% loadInto = [];

arbDconc = 1.0e-5;
lengthTrap = length(trapD.trapC(1,:));

if inA.jRun == 1 && (inA.iLeap == inA.nLeap+1 || inA.iLeap == 1)
        Crnt.m0(1:lengthTrap) = arbDconc;   % set to arbitrary small D concentration   
else 
    loadFrom  = trapD.order.load;
    parsed    = parse_TMAP_out(loadFrom);
    
    orderFrom = trapD.order.last;
    Crnt.m0 = parsed.m0;
    trapD.concF(orderFrom(1),:) = parsed.c1./trapD.trapC(orderFrom(1),:)/inA.eqn.N_w;
    trapD.concF(orderFrom(2),:) = parsed.c2./trapD.trapC(orderFrom(2),:)/inA.eqn.N_w;
    trapD.concF(orderFrom(3),:) = parsed.c3./trapD.trapC(orderFrom(3),:)/inA.eqn.N_w;

    trapD.concF(isnan(trapD.concF)) = 0 ;       % clean up NANs
    trapD.concF(isinf(trapD.concF)) = 0 ;       % clean up INFs
end

orderInto = trapD.order.crnt;

Crnt.t1 = trapD.trapC(orderInto(1),:); 
Crnt.t2 = trapD.trapC(orderInto(2),:);
Crnt.t3 = trapD.trapC(orderInto(3),:);
Crnt.c1 = trapD.concF(orderInto(1),:);
Crnt.c2 = trapD.concF(orderInto(2),:);
Crnt.c3 = trapD.concF(orderInto(3),:);    

trapD.Crnt = Crnt;

end

function [parsed] = parse_TMAP_out(outputPointer)
%
%
%

global folder

kOut = outputPointer(1);
iOut = outputPointer(2);

% print out for debugging purposes
% disp(folder.tmapNames{kOut,iOut})

fullName = [folder.Work folder.tmapNames{kOut,iOut} '.out']; % finds directory and appends 

% open file and parse first entry of each line
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

idx1 = find(strcmp(Var1,'Mobile Concentration (number/m**3)'),1,'last');
idx2 = find(strcmp(Var1,'Trap  1'                           ),1,'last');
idx3 = find(strcmp(Var1,'Trap  2'                           ),1,'last');
idx4 = find(strcmp(Var1,'Trap  3'                           ),1,'last');

nRows = (idx3-idx2)-1;

Var2 =       Var1(idx1+1:idx1+nRows);
Var2 = [Var2;Var1(idx2+1:idx2+nRows)];
Var2 = [Var2;Var1(idx3+1:idx3+nRows)];
Var2 = [Var2;Var1(idx4+1:idx4+nRows)];

nColumns = textscan(Var2{1},'%f');
nColumns = length(nColumns{1});

tempM = zeros(nRows*4,nColumns);

for jj = 1:nRows*4
    tempV = textscan(Var2{jj},'%f');
    len_V = length(tempV{1});
    dltaV = len_V - 6;                      % TMAP output file has 6 columns max
    if len_V > 6                            % Remove values below 1e-100 reported as 1-100 
        tempC = tempV{1};
        tempC = tempC(1:6 - dltaV);
        zeroC = zeros(dltaV,1);
        tempV{1} = [tempC;zeroC];
    end
    tempM(jj,:) = tempV{1};
end

parsed.m0 = tempM(0*nRows+1:1*nRows,:);
parsed.c1 = tempM(1*nRows+1:2*nRows,:);      % ?trap_conc.? * filled_frac. * atomic density
parsed.c2 = tempM(2*nRows+1:3*nRows,:);
parsed.c3 = tempM(3*nRows+1:4*nRows,:);

parsed.m0 = reshape(parsed.m0',1,[]);       % 
parsed.c1 = reshape(parsed.c1',1,[]);
parsed.c2 = reshape(parsed.c2',1,[]);
parsed.c3 = reshape(parsed.c3',1,[]);

end

%% TMAP Subroutines

function [] = Write_Input_File(inA,trapD)
% TMAP input file separated by section
% 1_title               (unchanged)
% 2_main                (unchanged)
% 3_enclosure           (unchanged)
% 4_thermal             (unchanged)
% 5_diffusion           (changes mobile, traps, and conc.)
% 6_equation            (changes mainly trapping energies 4,5,6)
% 7_table               (changes temperature ramp and start times)
% 8_control             (changes time, timend, and nprint)
% 9_plot                (unchanged)
global folder

Write_2_main_____(    trapD);
Write_4_thermal__(inA,trapD);
Write_5_diffusion(inA,trapD);   %  (changes mobile, inA, and conc.)
Write_6_equation_(inA      );   %  (changes mainly trapping energies 4,5,6)
Write_7_table____(inA      );   %  (changes temperature ramp and start times)
Write_8_control__(inA      );   %  (changes time, timend, and nprint)
Write_9_plot_____(inA      );   %  (changes time, timend, and nprint)

tmapName_k_i = folder.tmapNames{inA.iLeap,inA.iFit};

% utilize the modified control template for control input
if inA.implant.use_mod_template == 1
    if inA.jRun == 1||2 && inA.iLeap == inA.nLeap+1
        %'***** Using modified input file for temperature and flux'
        cmndLine = ['copy ' folder.Temp '7_table_MOD.inp' ' ' folder.Updt '7_table_Up.inp'];
        [~,~] = system(cmndLine);
    end
end

% Create new input file
catName = ['copy ' ...
            folder.Temp '1_title.inp + ' ...
            folder.Updt '2_main_Up.inp + ' ...
            folder.Temp '3_enclosure.inp + ' ...
            folder.Updt '4_thermal_Up.inp + ' ...
            folder.Updt '5_diffusion_Up.inp + ' ...
            folder.Updt '6_equation_Up.inp + ' ...
            folder.Updt '7_table_Up.inp + ' ...     % JHY this overwrites MOD file?
            folder.Updt '8_control_Up.inp + ' ...
            folder.Updt '9_plot_Up.inp ' ...
            folder.TMAP '' tmapName_k_i '.inp'];
[~,~] = system(catName);

end

function [] = Write_2_main_____(trapD)
%
%
%

global folder

fullName = [folder.Temp '2_main.inp']; % finds directory and appends 

% open file and parse first entry of each line
delimiter = {'', ',', '*', '=', '$'};
formatSpec = '%s %[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'Delimiter', delimiter,...
                     'ReturnOnError', false);
Var1 = strtrim(dataArray{:, 1});

fclose(fileID);

idx1 = find(strcmp(Var1,          'espcnme'),1,'first');
idx2 = find(strcmp(Var1,'end of main input'),1,'first');

% open file and parse first entry of each line
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

% % find indicies to cut and create new file
% idx1 = find(strcmp(Var1,'espcnme  =    d2g,end'  ),1,'first');
% idx2 = find(strcmp(Var1,'end of main input'  ),1,'first');

file_2 = [folder.Updt '2_main_Up.inp'];
fileID = fopen(file_2,'w');

segments = ['segnds   =     ' ...
            num2str(trapD.depth.segments,'%3.2d') ...
            ',end'];

format1  = '%s\n';

fprintf(fileID, format1, Var1{1:idx1});
fprintf(fileID, format1, segments);
fprintf(fileID, format1, Var1{idx2-1:end});

fclose(fileID);
end

function [] = Write_4_thermal__(inA,trapD)
% Open and parse template for 4_thermal.inp
%   Save new file with updated depth information
%    in 4_thermal_Up.inp

global folder

fullName = [folder.Temp '4_thermal.inp']; % finds directory and appends 

% open file and parse first entry of each line
delimiter = {'', ',', '*', '=', '$'};
formatSpec = '%s %[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'Delimiter', delimiter,...
                     'ReturnOnError', false);
Var1 = strtrim(dataArray{:, 1});

fclose(fileID);

% find indecies to cut and create new file
idx1 = find(strcmp(Var1,                'delx'),1,'first'); 
idx2 = find(strcmp(Var1,               'tempd'),1,'first'); 
idx3 = find(strcmp(Var1,                'hsrc'),1,'first'); 
idx4 = find(strcmp(Var1,'end of thermal input'),1,'first');

% reopen file and get each line of strings
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

% define thickness slices trapD.depth.T
file_4 = [folder.Updt '4_thermal_Up.inp'];

fileID = fopen(file_4,'w');

format1  = '%s\n';

depthT = cell(trapD.depth.segments-2,1);        % ignoring the zeros at begin and end
for ii = 1:trapD.depth.segments-2
depthT{ii} = num2str(trapD.depth.T(ii+1)*1e-6,'%20.4e');
end

jj = inA.iLeap;

init_temp = ['tempd =  ' ...
            num2str(  trapD.depth.segments,'%3.2d') '*'...
            num2str(inA.Table1.TempBgn(jj),'%6.1f') ',end'];
        
heat_src  = ['hsrc = const,0.,srcpf,' ...
            num2str(  trapD.depth.segments,'%3.2d') '*'...
            '0.,end'];
        
fprintf(fileID, format1, Var1{1:idx1});
fprintf(fileID, format1, depthT{1:end});
fprintf(fileID, format1, Var1{idx2-1});
fprintf(fileID, format1, init_temp);
fprintf(fileID, format1, Var1{idx2+1:idx3-1});
fprintf(fileID, format1, heat_src);
fprintf(fileID, format1, Var1{idx3+1:idx4});

fclose(fileID);
end

function [] = Write_5_diffusion(inA,trapD)
% Open and parse template for 5_diffusion.inp
%  Save new file with updated mobile D, traps, and filled concentrations
%    in 5_diffusion_Up.inp

global folder

fullName = [folder.Temp '5_diffusion.inp']; % finds directory and appends 

% open file and parse first entry of each line
delimiter = {'', ',', '*', '=', '$'};
formatSpec = '%s %[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'Delimiter', delimiter,...
                     'ReturnOnError', false);
Var1 = strtrim(dataArray{:, 1});

fclose(fileID);

% find indecies to cut and create new file
idx0 = find(strcmp(Var1(     1:end),  'nbrden'),1,'first') +    0; % number density
idx1 = find(strcmp(Var1(idx0+1:end),   'concd'),1,'first') + idx0; % mobile d
idx2 = find(strcmp(Var1(idx1+1:end),     'end'),1,'first') + idx1; % mobile d
idx3 = find(strcmp(Var1(idx2+1:end),'trapping'),1,'first') + idx2; % trap1
idx4 = find(strcmp(Var1(idx3+1:end),    'tspc'),1,'first') + idx3; % conc1
idx5 = find(strcmp(Var1(idx4+1:end),    'ttyp'),1,'first') + idx4; % trap2
idx6 = find(strcmp(Var1(idx5+1:end),    'tspc'),1,'first') + idx5; % conc2
idx7 = find(strcmp(Var1(idx6+1:end),    'ttyp'),1,'first') + idx6; % trap3
idx8 = find(strcmp(Var1(idx7+1:end),    'tspc'),1,'first') + idx7; % conc3
idx9 = find(strcmp(Var1(idx8+1:end),     'end'),1,'first') + idx8; % rest of file
idx10= find(strcmp(Var1(idx9+1:end),   'dcoef'),1,'first') + idx9; % rest of file

% reopen file and get each line of strings
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

fileTraps = [folder.Updt '5_diffusion_Up.inp'];

fileID = fopen(fileTraps,'w');

numberDensity = ['nbrden   =    ' num2str(inA.eqn.N_w,'%5.4e') ',end'];
concdBgn = trapD.Crnt.m0(1:end-1);
concdEnd = [num2str(trapD.Crnt.m0(end),'%20.4e') ',end']; 

format1  = '%s\n';
format2  = '%20.4e\n';
% concdEnd = [num2str(trapD.Crnt.m0(end),'%20.4e') ', end']; 

fprintf(fileID, format1, Var1{1:idx0-1 });
fprintf(fileID, format1, numberDensity );
fprintf(fileID, format1, Var1{idx1    });
fprintf(fileID, format2, concdBgn      );
fprintf(fileID, format1, concdEnd      );
fprintf(fileID, format1, Var1{idx3    });
fprintf(fileID, format2, trapD.Crnt.t1);
fprintf(fileID, format1, Var1{idx4    });
fprintf(fileID, format1,'       alpht,equ,3,ctrap');
fprintf(fileID, format2, trapD.Crnt.c1);
fprintf(fileID, format1, Var1{idx5    });
fprintf(fileID, format2, trapD.Crnt.t2);
fprintf(fileID, format1, Var1{idx6    });
fprintf(fileID, format1,'       alpht,equ,3,ctrap');
fprintf(fileID, format2, trapD.Crnt.c2);
fprintf(fileID, format1, Var1{idx7    });
fprintf(fileID, format2, trapD.Crnt.t3);
fprintf(fileID, format1, Var1{idx8    });
fprintf(fileID, format1,'       alpht,equ,3,ctrap');
fprintf(fileID, format2, trapD.Crnt.c3);
fprintf(fileID, format1, Var1{idx9:idx10});
if inA.iLeap == inA.nLeap+1 && (inA.jRun == 1 || inA.jRun == 2)
    source = ['srcsd  = d,tabl,2,srcpf,norm,' ...
              num2str(inA.implantGauss.trans_,'%4.3e') ','... 
              num2str(inA.implantGauss.center,'%4.3e') ','... 
              num2str(inA.implantGauss.width_,'%4.3e') ','... 
              num2str(inA.implantGauss.offset,'%4.3e') ',end'];
    fprintf(fileID, format1, source);
end
fprintf(fileID, format1, Var1{idx10+1:end});



fclose(fileID);
end

function [] = Write_6_equation_(inA)
%
%
%

global folder

fullName = [folder.Temp '6_equation.inp']; % finds directory and appends 

% open file and parse first entry of each line
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

% find indicies to cut and create new file
idx1 = find(strcmp(Var1,'$ (4)  Alphr for trap 1 in tungsten (1/s)'  ),1,'first');
idx2 = find(strcmp(Var1,'$ (5)  Alphr for trap 2 in tungsten (1/s)'  ),1,'first');
idx3 = find(strcmp(Var1,'$ (6)  Alphr for trap 3 in tungsten (1/s)'  ),1,'first');
idx4 = find(strcmp(Var1,'$ (7)  Diffusivity for d in tungsten (m2/s)'),1,'first');

file_6 = [folder.Updt '6_equation_Up.inp'];
fileID = fopen(file_6,'w');

nu = inA.eqn.jmpF;
kb = inA.eqn.k_B;

kPick = inA.order(inA.iLeap,:);
E_t = inA.trap.Matrix(kPick,4);
yTrap = cell(1,3);

% Use all traps in pseudo or only lowest energy for release rate
lowest_E = 0;

if lowest_E == 1
    for ii = 1:3
        if kPick(ii) == inA.nTraps+1    % i.e. the Psuedo trap
            if inA.iLeap == inA.nLeap+1  % i.e. the implantation phase
                pStart = 3;
            else
                pStart = sort(kPick);   % Sort the traps
                pStart = pStart(2)+1;   % Choose highest value other than pseudo trap 
            end
            sumC = inA.trap.Matrix(:,2) + inA.trap.Matrix(:,3);
            sumCt = sum(sumC(pStart:inA.nTraps));
            Energy = inA.trap.Matrix(:,4);
            yPseudoR = [   'y = ' num2str(            nu,'%3.2e') ...
                         '*exp(-' num2str(Energy(pStart),'%5.3f') ...
                              '/' num2str(            kb,'%5.3e') ...
                        '/temp)*' num2str(  sumC(pStart),'%5.3e') ...
                              '/' num2str(         sumCt,'%5.3e') ...
                                                          ',end'];    
        else
            yTrap{ii} = [  'y = ' num2str(      nu,'%3.2e') ...
                         '*exp(-' num2str( E_t(ii),'%5.4f') ...
                              '/' num2str(      kb,'%5.4e') ...
                                              '/temp),end'];
        end
    end

    format1  = '%s\n';

    fprintf(fileID, format1, Var1{1:idx1});
    if isempty(yTrap{1})
        fprintf(fileID, format1, yPseudoR);
    else
        fprintf(fileID, format1, yTrap{1});
    end
    fprintf(fileID, format1, Var1{idx2});
    if isempty(yTrap{2})
        fprintf(fileID, format1, yPseudoR);
    else
        fprintf(fileID, format1, yTrap{2});
    end
    fprintf(fileID, format1, Var1{idx3});
    if isempty(yTrap{3})
        fprintf(fileID, format1, yPseudoR);
    else
        fprintf(fileID, format1, yTrap{3});
    end
    fprintf(fileID, format1, Var1{idx4:end});    
       
else

    for ii = 1:3
        if kPick(ii) == inA.nTraps+1    % i.e. the Psuedo trap
            if inA.iLeap == inA.nLeap+1 % i.e. the implantation phase
                pStart = 3;
            else
                pStart = sort(kPick);   % Sort the traps
                pStart = pStart(2)+1;   % Choose highest value other than pseudo trap 
            end
            pTotal = inA.nTraps-pStart+1;
            sumC = inA.trap.Matrix(:,2) + inA.trap.Matrix(:,3);
            sumCt = sum(sumC(pStart:inA.nTraps));
            Energy = inA.trap.Matrix(:,4);
            yPseudoR = cell(1,pTotal+2);
            yPseudoR{1}   = [ 'y = ' num2str(      nu,'%3.2e') ...
                                 '/' num2str(   sumCt,'%5.4e') '*('];
            yPseudoR{end} = '0),end';
        
            kk = 1;
            for jj = pStart:inA.nTraps
                kk = kk + 1;
                yPseudoR{kk} = [  'exp(-' num2str(Energy(jj),'%5.4f') ...
                                      '/' num2str(        kb,'%5.4e') ...
                                '/temp)*' num2str(  sumC(jj),'%5.4e') ...
                                                              ' + ' ];
            
            end     
        else
            yTrap{ii} = [  'y = ' num2str(      nu,'%3.2e') ...
                     '*exp(-' num2str( E_t(ii),'%5.4f') ...
                          '/' num2str(      kb,'%5.4e') ...
                     '/temp),end'];
        end
    end

    format1  = '%s\n';
    fprintf(fileID, format1, Var1{1:idx1});
    if isempty(yTrap{1})
        for ii = 1:length(yPseudoR)
            fprintf(fileID, format1, yPseudoR{ii});
        end
    else
        fprintf(fileID, format1, yTrap{1});
    end
    fprintf(fileID, format1, Var1{idx2});
    if isempty(yTrap{2})
        for ii = 1:length(yPseudoR)
            fprintf(fileID, format1, yPseudoR{ii});
        end
    else
        fprintf(fileID, format1, yTrap{2});
    end
    fprintf(fileID, format1, Var1{idx3});
    if isempty(yTrap{3})
        for ii = 1:length(yPseudoR)
            fprintf(fileID, format1, yPseudoR{ii});
        end
    else
        fprintf(fileID, format1, yTrap{3});
    end
    fprintf(fileID, format1, Var1{idx4:end});

end


fclose(fileID);

end

function [] = Write_7_table____(inA)
%
%
%

global folder

fullName = [folder.Temp '7_table.inp']; % finds directory and appends 

% open file and parse first entry of each line
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};

fclose(fileID);

% find indecies to cut and create new file
idx1 = find(strcmp(Var1,'$ (1) Temperature history of enclosure 1, plasma chamber (K)'  ),1,'first');
idx2 = find(strcmp(Var1,'$ (2) Implantation flux history (atom/m2/s)'),1,'first');
idx3 = find(strcmp(Var1,'end of table input'),1,'first');

file_7 = [folder.Updt '7_table_Up.inp'];
fileID = fopen(file_7,'w');

jj = inA.iLeap;

timeTemp_0 = ['      ' num2str(inA.Table1.timeBgn(jj),'%6.1f') ...
                  ', ' num2str(inA.Table1.TempBgn(jj),'%6.1f')];
timeTemp_1 = ['      ' num2str(inA.Table1.timeEnd(jj),'%6.1f') ...
                  ', ' num2str(inA.Table1.TempEnd(jj),'%6.1f') ...
              ',  end'];

if inA.jRun == 3
    timeTemp_0 = ['      ' num2str(inA.Table1.timeBgn(jj),'%6.1f') ...
                      ', ' num2str(inA.Table1.TempBgn(jj),'%6.1f')];
    timeTemp_1 = ['      ' num2str(inA.implantFlux(4,1),  '%6.1f') ...
                      ', ' num2str(                   293,'%6.1f') ...
                  ',  end'];   
end

if inA.iLeap == inA.nLeap+1 && inA.jRun == 4
    timeTemp_0 = ['      ' num2str(inA.Table1.timeBgn(jj),'%6.1f') ...
                      ', ' num2str(                   293,'%6.1f')];
    timeTemp_1 = ['      ' num2str(inA.implantFlux(4,1),  '%6.1f') ...
                      ', ' num2str(                   293,'%6.1f') ...
                  ',  end'];
end


format1  = '%s\n';

fprintf(fileID, format1, Var1{1:idx1});
fprintf(fileID, format1, timeTemp_0);
fprintf(fileID, format1, timeTemp_1);
fprintf(fileID, format1, Var1{idx2});
if inA.iLeap == inA.nLeap+1 && (inA.jRun == 1 || inA.jRun == 2)
    % May need to use if else to not include for TDS segments (maybe ignored)          
    timeImpl_0 = ['      ' num2str(inA.implantFlux(1,1),'%6.1f') ...
                      ', ' num2str(inA.implantFlux(1,2),'%6.2e')];
    timeImpl_1 = ['      ' num2str(inA.implantFlux(2,1),'%6.1f') ...
                      ', ' num2str(inA.implantFlux(2,2),'%6.2e') ...
                  ',  end'];
              
    fprintf(fileID, format1, timeImpl_0);
    fprintf(fileID, format1, timeImpl_1);
end
fprintf(fileID, format1, Var1{idx3});

fclose(fileID);
end

function [] = Write_8_control__(inA)
%
%
%

global folder

fullName = [folder.Temp '8_control.inp']; % finds directory and appends 

% open file and parse first entry of each line
delimiter = {'', ',', '*', '=', '$'};
formatSpec = '%s %[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'Delimiter', delimiter,...
                     'ReturnOnError', false);
Var1 = strtrim(dataArray{:, 1});
fclose(fileID);


% find indecies to cut and create new file
idx1 = find(strcmp(Var1,'time'  ),1,'first');
idx2 = find(strcmp(Var1,'bump'),1,'first');

% reopen and save entire line
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};
fclose(fileID);


file_8 = [folder.Updt '8_control_Up.inp'];
fileID = fopen(file_8,'w');

jj = inA.iLeap;

time_   = ['time   = ' num2str(inA.Table1.timeBgn(jj)  ,'%6.1f') ', end'];
tstep_  = ['tstep  = ' num2str(inA.tstep(jj)           ,'%6.3e') ', end'];
timend_ = ['timend = ' num2str(inA.Table1.timeEnd(jj)  ,'%6.3f') ', end'];
nprint_ = ['nprint = ' num2str(inA.Table1.n_print(jj)-1,'%6d'  ) ', end'];
itermx  = ['itermx = ' num2str(inA.control.itermx      ,'%6d'  ) ', end'];
delcmx  = ['delcmx = ' num2str(inA.control.delcmx      ,'%6.3e') ', end'];      % delta convergence max

if inA.jRun == 1
    tstep1  = inA.tstep_initial;
    timeEnd = 0.10;
    nstep1  = ceil(timeEnd/tstep1);                               % 0.10 second / tstep1
    tstep_  = ['tstep  = ' num2str(   tstep1,'%6.3e') ', end'];
    timend_ = ['timend = ' num2str(  timeEnd,'%6.3e') ', end'];
    nprint_ = ['nprint = ' num2str(-1+nstep1,'%6d'  ) ', end'];
elseif inA.jRun == 3 || inA.jRun == 4
    slowStep = inA.tstep(jj)/5;
    tstep_  = ['tstep  = ' num2str(slowStep,'%6.3e') ', end'];
    timeEnd = inA.implantFlux(4,1);
    timend_ = ['timend = ' num2str(timeEnd,'%6.1f')    ', end'];
    nprint_ = ['nprint = ' num2str(floor(timeEnd/slowStep-1),'%6d') ', end'];
end


format1  = '%s\n';
format2  = '  %s\n';

fprintf(fileID, format1, Var1{1:idx1-1});
fprintf(fileID, format2, time_   );
fprintf(fileID, format2, tstep_  );
fprintf(fileID, format2, timend_ );
fprintf(fileID, format2, nprint_ );
fprintf(fileID, format2, itermx  );
fprintf(fileID, format2, delcmx  );

fprintf(fileID, format2, Var1{idx2:end});

fclose(fileID);
end

function [] = Write_9_plot_____(inA)
%
%
%

global folder

fullName = [folder.Temp '9_plot.inp']; % finds directory and appends 

% open file and parse first entry of each line
delimiter = {'', ',', '*', '=', '$'};
formatSpec = '%s %[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'Delimiter', delimiter,...
                     'ReturnOnError', false);
Var1 = strtrim(dataArray{:, 1});
fclose(fileID);


% find indecies to cut and create new file
idx1 = find(strcmp(Var1,'nplot'  ),1,'first');
idx2 = find(strcmp(Var1,'plotseg'),1,'first');

% reopen and save entire line
formatSpec = '%[^\n\r]';
fileID = fopen(fullName,'r');
dataArray = textscan(fileID, formatSpec,...
                     'ReturnOnError', false);
Var1 = dataArray{:, 1};
fclose(fileID);

file_9 = [folder.Updt '9_plot_Up.inp'];
fileID = fopen(file_9,'w');

nPlot = ceil(inA.plot.tplot/inA.control.tstep_TDS);

nplot = ['nplot = ' num2str(nPlot,'%6d') ', end'];

format1  = '%s\n';
format2  = '  %s\n';

fprintf(fileID, format1, Var1{1:idx1-1});
fprintf(fileID, format2, nplot   );
fprintf(fileID, format2, Var1{idx2:end});

fclose(fileID);
end
