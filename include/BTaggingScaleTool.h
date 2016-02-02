#ifndef __BTAGGINGSCALETOOL_H__
#define __BTAGGINGSCALETOOL_H__

// SFrame include(s):
#include "core/include/SError.h"
#include "plug-ins/include/SToolBase.h"

#include "../Common/D3PDVariables/include/Jet.h"

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

  //  double getScaleFactor( const DESY::Jet& jet, double sigma_shift = 0. );
  double getScaleFactor( const DESY::Jet& jet, double sigma = 0.);


  //  double getScaleFactor( const std::vector< DESY::Jet >& vJets, double sigma_shift );
  //  double getScaleFactor( const DESY::JetVec& vJets, double sigma_shift = 0.);
  double getScaleFactor( const DESY::JetVec& vJets, double sigma = 0.);

 private:

  std::string m_name;                 ///< name of the tool
  std::string m_tagger;
  std::string m_workingPoint;
  std::string m_csvFile;
  std::string m_measurementType;

  BTagCalibration*       m_calib;
  BTagCalibrationReader* m_reader;
  BTagCalibrationReader* m_reader_up;
  BTagCalibrationReader* m_reader_down;

};


#endif //  __BTAGGINGSCALETOOL_H__
