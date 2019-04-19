#!/bin/bash

#PARAMETERS TO BE SET
################################################################################################
Cent="k010"

cfgFileData="configfiles/config_Ds_data_010.yml"
cfgFileMC="configfiles/config_Ds_MC_010.yml"
cfgFileFit="configfiles/config_Ds_Fit.yml"

accFileName="accfiles/Acceptance_Toy_DsKKpi_yfidPtDep_etaDau09_ptDau100_FONLL5ptshape.root"
predFileName="models/D0DplusDstarPredictions_502TeV_y05_all_021016_BDShapeCorrected.root"
pprefFileName="ppreference/Ds_ppreference_pp5TeV_noyshift_pt_2_3_4_5_6_8_12_16_24_36_50.root"

PtWeightsFileName="ptweights/PtWeigths_LHC19c3a.root"
PtWeightsHistoName="hPtWeightsFONLLtimesTAMUcent"

#assuming cutsets config files starting with "cutset" and are .yml
CutSetsDir="configfiles"
declare -a CutSets=("_010_central_2018")
arraylength=${#CutSets[@]}

OutDirRawyields="outputs/rawyields"
OutDirEfficiency="outputs/efficiency"
OutDirCrossSec="outputs/crosssec"
OutDirRaa="outputs/raa"
################################################################################################

if [ ! -f "${cfgFileData}" ]; then
  echo ERROR: data config file "${cfgFileData}" does not exist!
  exit 2
fi

if [ ! -f "${cfgFileMC}" ]; then
  echo ERROR: MC config file "${cfgFileMC}" does not exist!
  exit 2
fi

if [ ! -f "${cfgFileFit}" ]; then
  echo ERROR: it config file "${cfgFileFit}"does not exist!
  exit 2
fi

if [ ! -f "${accFileName}" ]; then
  echo ERROR: acceptance file "${accFileName}" does not exist!
  exit 2
fi

if [ ! -f "${predFileName}" ]; then
  echo ERROR: FONLL file "${predFileName}" does not exist!
  exit 2
fi

if [ ! -f "${pprefFileName}" ]; then
  echo ERROR: pp reference file "${pprefFileName}" does not exist!
  exit 2
fi

if [ ! -f "${PtWeightsFileName}" ]; then
  echo WARNING: pT-weights file "${PtWeightsFileName}" does not exist!
fi

if [ ! -d "${OutDirRawyields}" ]; then
  mkdir ${OutDirRawyields}
fi

if [ ! -d "${OutDirEfficiency}" ]; then
  mkdir ${OutDirEfficiency}
fi

if [ ! -d "${OutDirCrossSec}" ]; then
  mkdir ${OutDirCrossSec}
fi

if [ ! -d "${OutDirRaa}" ]; then
  mkdir ${OutDirRaa}
fi

#project sparses
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Projecting data sparse 
  python ProjectDplusDsSparse.py ${cfgFileData} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml ${OutDirRawyields}/Distr_Ds_data${CutSets[$iCutSet]}.root
  echo Projecting MC sparses 
  python ProjectDplusDsSparse.py ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_Ds_MC${CutSets[$iCutSet]}.root
done

# #compute raw yields
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Extract raw yields from ${OutDirEfficiency}/Distr_Ds_MC${CutSets[$iCutSet]}.root
  echo '.x GetRawYieldsDplusDs.C+ ('${Cent}',true, "'${OutDirEfficiency}'/Distr_Ds_MC'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYieldsDs_MC'${CutSets[$iCutSet]}'.root")' | root -l
  echo '.q'
  echo Extract raw yields from ${OutDirRawyields}/Distr_Ds_data${CutSets[$iCutSet]}.root
  echo '.x GetRawYieldsDplusDs.C+ ('${Cent}',false, "'${OutDirRawyields}'/Distr_Ds_data'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYieldsDs'${CutSets[$iCutSet]}'.root")' | root -l
  echo '.q'
done

#compute efficiency
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Compute efficiency from ${OutDirEfficiency}/Distr_Ds_MC${CutSets[$iCutSet]}.root
  python ComputeEfficiencyDplusDs.py ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml ${OutDirEfficiency}/Distr_Ds_MC${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Efficiency_Ds${CutSets[$iCutSet]}.root --ptweights ${PtWeightsFileName} ${PtWeightsHistoName}
done

#compute efficiency times acceptance
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Compute efficiency times acceptance
  python CombineAccTimesEff.py ${OutDirEfficiency}/Efficiency_Ds${CutSets[$iCutSet]}.root ${accFileName} ${OutDirEfficiency}/Eff_times_Acc_Ds${CutSets[$iCutSet]}.root
done

#compute HFPtSpectrum
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Compute HFPtspectrum
  echo '.x HFPtSpectrum.C+ (kDsKKpi,"'${predFileName}'","'${OutDirEfficiency}'/Eff_times_Acc_Ds'${CutSets[$iCutSet]}'.root","'${OutDirRawyields}'/RawYieldsDs'${CutSets[$iCutSet]}'.root","hRawYields","hAccEffPrompt","hAccEffFD","hEvForNorm","'${OutDirCrossSec}'/HFPtSpectrumDs'${CutSets[$iCutSet]}'.root",kNb,1.,true,'${Cent}',k2018)' | root -l
  echo '.q'
done

#compute HFPtSpectrumRaa
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Compute HFPtspectrumRaa
  echo '.x HFPtSpectrumRaa.C+ ("'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrumDs'${CutSets[$iCutSet]}'.root","'${OutDirRaa}'/HFPtSpectrumRaaDs'${CutSets[$iCutSet]}'.root",4,1,kNb,'${Cent}',k2018,k5dot023,1./3,3,6,false,1)' | root -l
  echo '.q'
done

#compute yield
for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
do
  echo Compute corrected yield
  echo '.x ComputeDmesonYield.C+ (kDs,'${Cent}',2,1,'${pprefFileName}',"'${OutDirCrossSec}'/HFPtSpectrumDs'${CutSets[$iCutSet]}'.root","","'${OutDirRaa}'/HFPtSpectrumRaaDs'${CutSets[$iCutSet]}'.root","",'${OutDirCrossSec}',1,1./3,3,1)' | root -l
  echo '.q'
done