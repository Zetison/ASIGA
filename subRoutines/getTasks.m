counter = 1;
studies = cell(0,1);

getDefaultTaskValues

if ~exist('studyName','var')
%     studyName = 'getTask_test';
%     studyName = 'getTask_test2';
%     studyName = 'getTask_Simpson';
%     studyName = 'getTask_Simpson_Torus';
%     studyName = 'getTask_Simpson_Torus_gp';
%     studyName = 'getTask_Simpson_PS';
%     studyName = 'getTask_Simpson_sweep';
%     studyName = 'getTask_BeTSSi_IIb';
%     studyName = 'getTask_M5';
%     studyName = 'getTask_M4';
%     studyName = 'getTask_M3';
%     studyName = 'getTask_M3_sweep';
%     studyName = 'getTask_PH';
%     studyName = 'getTask_KarlThomas';
%     studyName = 'getTask_MockShell';
%     studyName = 'getTask_MockShell_P';
%     studyName = 'getTask_MockShell_P_convergence';
%     studyName = 'getTask_MockShell_P_L';
%     studyName = 'getTask_S';
%     studyName = 'getTask_S1_COMSOL';
%     studyName = 'getTask_S1_Demkowicz';
%     studyName = 'getTask_S1_MFS';
%     studyName = 'getTask_S1_sweep';
%     studyName = 'getTask_MA8501';
%     studyName = 'getTask_SspectralConvergence';
%     studyName = 'getTask_S_sweep';
%     studyName = 'getTask_S_sweepBA';
%     studyName = 'getTask_S_MS';
%     studyName = 'getTask_S5';
%     studyName = 'getTask_Fillinger';
%     studyName = 'getTask_testIE';
%     studyName = 'getTask_Hetmaniuk';
%     studyName = 'getTask_articleIGA';
%     studyName = 'getTask_articleIGA_convergenceAnalysis';
%     studyName = 'getTask_articleIGA_Ihlenburg';
%     studyName = 'getTask_articleIGA_Ihlenburg2';
%     studyName = 'getTask_articleIGA_Ihlenburg3';
%     studyName = 'getTask_articleBEM_S1';
%     studyName = 'getTask_articleBEM_S1_convergenceAnalysis';
%     studyName = 'getTask_articleBEM_S1_sweep';
%     studyName = 'getTask_articleBEM_S1parmComp';
%     studyName = 'getTask_articleBEM_S1_qp';
%     studyName = 'getTask_articleBEM_S1_cg';
%     studyName = 'getTask_articleBEM_S1_BA';
%     studyName = 'getTask_articleSBEM_S1';
%     studyName = 'getTask_articleBEM_S1manufactored';
%     studyName = 'getTask_articleBEM_BCA';
%     studyName = 'getTask_articleBEM_BCA_MS';
%     studyName = 'getTask_articleBEM_BCA_P';
%     studyName = 'getTask_articleBEM_BCA_P_sweep';
%     studyName = 'getTask_articleBEM_BCA_BA';
%     studyName = 'getTask_articleBEM_BCA_Pgp';
%     studyName = 'getTask_BC_P_paraview';
%     studyName = 'getTask_BC_P';
%     studyName = 'getTask_BC';
%     studyName = 'getTask_BC_BEM';
%     studyName = 'getTask_BCA';
%     studyName = 'getTask_BCA_P';
%     studyName = 'getTaskGerdes1998tcv';
%     studyName = 'getTaskGerdes1998tcv_P';
end
eval(studyName)
