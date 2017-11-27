
%First load inPUT file to use as template

%Next, load the TDS data:
load('C:\Users\Jonathan\Documents\TDS_data\WD_melt\WD_00_103_111_112.mat');

%Define the sample name
name = WD112_srs;

%Fill in the inPUT structure fields for Perform_Fit_PTTP_vX
inPUT.sample.TDS.tempK = name.temptr(name.jSpan);
inPUT.sample.TDS.Dflux = name.D_flux(name.jSpan);
inPUT.sample.TDS.D_total = name.Dreten;

%Optional: save the inPUT structure in a .mat file

%Run the script:
%  outPUT = Perform_Fit_PTTP_v4(inPUT,[1,3,4],[10,10,2]);