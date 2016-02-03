#include "include/BTaggingScaleTool.h"

#include <cstdlib>
#include <limits>

#include "../include/BTagCalibrationStandalone.h"

//
// constructor
//
BTaggingScaleTool::BTaggingScaleTool( SCycleBase* parent, 
                                      const char* name ) : 
  SToolBase( parent ), m_name( name ) {

  SetLogName( name );

  std::string sframe_dir(std::getenv("SFRAME_DIR"));

  m_calib = 0;
  m_reader = 0;
  m_reader_up = 0;
  m_reader_down = 0;

  DeclareProperty( m_name + "_Tagger",    m_tagger = "CSVv2" );
  DeclareProperty( m_name + "_WorkingPoint", m_workingPoint = "Loose" );
  DeclareProperty( m_name + "_CsvFile", m_csvFile = sframe_dir + "/../BTaggingTools/csv/CSVv2.csv" );
  DeclareProperty( m_name + "_MeasurementType", m_measurementType = "comb" );
}

//
// destructor
//
BTaggingScaleTool::~BTaggingScaleTool() {
  delete m_calib;
  delete m_reader;
  delete m_reader_up;
  delete m_reader_down;
}

void BTaggingScaleTool::BeginInputData( const SInputData& ) throw( SError ) {

  m_logger << INFO << "Initializing BTagCalibrationStandalone" << SLogger::endmsg;
  m_logger << INFO << "CSV file:    " << m_csvFile << SLogger::endmsg;
  m_logger << INFO << "Tagger:         " << m_tagger << SLogger::endmsg;
  m_logger << INFO << "WorkingPoint: " << m_workingPoint << SLogger::endmsg;
  m_logger << INFO << "MeasurementType: " << m_measurementType << SLogger::endmsg;
  
  BTagEntry::OperatingPoint wp = BTagEntry::OP_LOOSE;
  if (m_workingPoint == "Loose") {
    wp = BTagEntry::OP_LOOSE;
  }
  else if (m_workingPoint == "Medium") {
    wp = BTagEntry::OP_MEDIUM;
  }
  else if (m_workingPoint == "Tight") {
    wp = BTagEntry::OP_TIGHT;
  }
  else {
    throw SError( ("Unknown working point: " + m_workingPoint).c_str(), SError::SkipCycle );
  }

  m_calib = new BTagCalibration(m_tagger, m_csvFile);
  m_reader = new BTagCalibrationReader(m_calib,    // calibration instance
                             wp,               // working point
                             m_measurementType, // measurement type
                             "central");        // systematics type
  m_reader_up = new BTagCalibrationReader(m_calib, wp, m_measurementType, "up");
  m_reader_down = new BTagCalibrationReader(m_calib, wp, m_measurementType, "down");

  return;

}
                                                                                    
double BTaggingScaleTool::getScaleFactor( const UZH::Jet& jet, double sigma ) {

  // Flavor
  BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;
  if (jet.hadronFlavour() == 5) {
    flavor = BTagEntry::FLAV_B;
  }
  else if (jet.hadronFlavour() == 4) {
    flavor = BTagEntry::FLAV_C;
  }
  else if (jet.hadronFlavour() == 15) {
    flavor = BTagEntry::FLAV_C; // Use C scale factors for tau for now.
  }

  m_logger << DEBUG << "     flavor " << flavor << SLogger::endmsg;

  // range checking, double uncertainty if beyond
  double MaxBJetPt = 670., MaxLJetPt = 1000., MaxEta = 2.4;
  
  if (abs(jet.eta()) > 2.4) {
    return 1.;
  }

  if (jet.hadronFlavour() == BTagEntry::FLAV_B) {
    if (jet.pt() > MaxBJetPt) {
      sigma *= sigma*2;
    }
  }
  else {
    if (jet.pt() > MaxLJetPt) {
      sigma *= sigma*2;
    }
  }

  // Note: this is for b jets, for c jets (light jets) use FLAV_C (FLAV_UDSG)
  double scalefactor = m_reader->eval(flavor, jet.eta(), jet.pt());
  if ((sigma < std::numeric_limits<double>::epsilon()) && (sigma > -std::numeric_limits<double>::epsilon())) {
    if (sigma > 0) {
      double scalefactor_up =  m_reader_up->eval(flavor, jet.eta(), jet.pt()); 
      scalefactor = sigma*(scalefactor_up - scalefactor) + scalefactor;
    }
    else {
      double scalefactor_down =  m_reader_down->eval(flavor, jet.eta(), jet.pt()); 
      scalefactor = abs(sigma)*(scalefactor_down - scalefactor) + scalefactor;
    }
  }
  
  double jetweight = 1.;
  // set effMC to one for now, need to use real value map later
  double effMC = 1.;
  
  if (jet.isTagged()) {
    m_logger << DEBUG << "     Jet is tagged " << SLogger::endmsg;
    jetweight *= scalefactor;
  }
  else {
    m_logger << DEBUG << "     Jet is not tagged " << SLogger::endmsg;
    jetweight *= (1 - (scalefactor * effMC)) / (1 - effMC);
  }
  
  

  m_logger << DEBUG << " jetweight " << jetweight << SLogger::endmsg;

  return jetweight;
}

//
// return scale for Jet collection
//
double BTaggingScaleTool::getScaleFactor( const UZH::JetVec& vJets, double sigma ) {

  double scale = 1.;
  
  m_logger << DEBUG << "BTaggingScaleTool::getScaleFactor" << SLogger::endmsg;

  for (std::vector< UZH::Jet>::const_iterator itJet = vJets.begin(); itJet < vJets.end(); ++itJet) {
    m_logger << DEBUG << "Looking at jet " << itJet - vJets.begin()
	     << ", pT=" << (*itJet).pt() << ", eta=" << (*itJet).eta()
	     << SLogger::endmsg;

    scale *= getScaleFactor(*itJet, sigma);
  }  

  m_logger << DEBUG << "BTaggingScaleTool::getScaleFactor done" << SLogger::endmsg;
  return scale;

}


