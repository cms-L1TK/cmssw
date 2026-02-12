# configuration of Demonstrator. This is used to compare FW with SW for the subset fo the chain between LabelIn & LabelOut. FW must be wrapped by EMP & compiled with IPBB.
import FWCore.ParameterSet.Config as cms

from L1Trigger.TrackFindingTracklet.Producer_cfi import TrackFindingTrackletProducer_params
from L1Trigger.TrackFindingTracklet.Analyzer_cfi import TrackFindingTrackletAnalyzer_params

# these parameters are for ModelSim runs of FW
TrackTriggerDemonstrator_params = cms.PSet (

  LabelIn  = TrackFindingTrackletProducer_params.InputLabelTM, #
  LabelOut = TrackFindingTrackletAnalyzer_params.OutputLabelTFP, #
  DirIPBB  = cms.string( "/users/sb19423/integration_test/tf_work/proj/tp_qsim/" ), # path to ipbb proj qsim area
  RunTime  = cms.double( 6. ),                                   # runtime in us

  LinkMappingIn  = cms.vint32(  ),
  LinkMappingOut = cms.vint32( 88, 89 ),

  # For L1T output (only relevant if OutputLabelTFP is used)
  MaxEventsPerL1TFile = cms.int32( 7 ),  # If non-zero: create testvector compatible with L1T from TFP outputs - each tracker region gets its own set of links

  # For HW tests. Following requirements need to be met:
  # 0. Build FW which matches the above in/out labels and the link mapping
  # 1. Apollo herd is running on the board
  # 2. Make sure singularity exists on the machine you're running the demonstrator from
  # 3. If outside CERN network, in a separete terminal tunnel to the board via lxplus e.g. 'ssh -L 3000:CMSAPOLLO214:3000 CERNUSERNAME@lxplus.cern.ch'
  DirTGZ = cms.string( "" ), # path to existing HW tarball, empty if you want to skip HW tests
  # DirTGZ = cms.string( "/users/sb19423/integration_test/tf_work/proj/tp_apollo/package/tp_apollo_hm01_dice_priv_260108_2316.tgz" ), # path to existing HW tarball, empty if you want to skip HW tests
  DockerImage = cms.string( "docker://gitlab-registry.cern.ch/cms-cactus/phase2/pyswatch/centos/algo-test-runner:v0.5.2" ), # Docker image used for submitting HW jobs
  BoardAddress = cms.string( "localhost" ), # DNS name of the board to run on, e.g. "APOLLO214". Use "localhost" if outside CERN network and follow step 3 above.
  BoardType = cms.string( "apollo" ), # `apollo`` or `serenity`
  NameFPGA = cms.string( "F1" ), # name of the FPGA on the board (`F1`` or `F2`` on apollo, `x0` or `x1` on serenity)
  BufferOffset = cms.int32( 611 ), # tx buffer offeset, i.e. when to start capturing the HW output data to match with CMSSW
  SaveAllFiles = cms.bool(True) # save all in/out/diff/pre txt files for each event

)
