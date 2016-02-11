#ifndef __BTAGGINGSCALETOOL_H__
#define __BTAGGINGSCALETOOL_H__

// SFrame include(s):
#include "core/include/SError.h"
#include "plug-ins/include/SToolBase.h"

#include "../NtupleVariables/include/Jet.h"

#include "../include/BTagCalibrationStandalone.h"

class BTaggingScaleTool : public SToolBase {
  
  //
  // Follow examples in https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
  // 
  //
  
 public:
  /// constructor
  BTaggingScaleTool( SCycleBase* parent, const char* name = "BTaggingScaleTool" );

  /// destructor
  ~BTaggingScaleTool();

  /// function booking histograms
  void BeginInputData( const SInputData& id ) throw( SError );
  
  double getScaleFactor( const double& pt, const double& eta, const int& flavour, bool isTagged, const double& sigma = 0.);

  double getScaleFactor( const UZH::Jet& jet, const double& sigma = 0.);
  
  double getPrunedSubjetScaleFactor( const UZH::Jet& jet, const double& sigma = 0.);

  double getScaleFactor( const UZH::JetVec& vJets, const double& sigma = 0.);
  
  double getPrunedSubjetScaleFactor( const UZH::JetVec& vJets, const double& sigma = 0.);
  
  /// function to book histograms for efficiencies
  void bookHistograms();
  
  /// function to fill jet b-tagging efficiencies
  void fillEfficiencies( const UZH::JetVec& vJets );
  
  /// function to fill subjet b-tagging efficiencies
  void fillPrunedSubjetEfficiencies( const UZH::JetVec& vJets );
  
  /// function to convert flavor integer to TString
  TString flavourToString( const int& flavour );
  
  /// helper function to check if jet is b-tagged
  bool isTagged( const UZH::Jet& jet );
  
  /// helper function to check if jet is b-tagged passing CSV value directly
  bool isTagged( const double& csv );
  

 private:

  std::string m_name;                 ///< name of the tool
  std::string m_tagger;
  std::string m_workingPoint;
  std::string m_csvFile;
  std::string m_measurementType_udsg;
  std::string m_measurementType_bc;
  std::string m_effHistDirectory;
  std::string m_effFile;
  
  std::map<std::string, double> wpCuts; // could have a function to set these
  double currentWorkingPointCut;

  std::unique_ptr<BTagCalibrationReader> m_reader;
  std::unique_ptr<BTagCalibrationReader> m_reader_up;
  std::unique_ptr<BTagCalibrationReader> m_reader_down;

};


#endif //  __BTAGGINGSCALETOOL_H__
