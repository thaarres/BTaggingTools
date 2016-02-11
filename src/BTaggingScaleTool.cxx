#include "include/BTaggingScaleTool.h"

#include <cstdlib>
#include <limits>

#include <TH2.h>

//
// constructor
//
BTaggingScaleTool::BTaggingScaleTool( SCycleBase* parent, 
                                      const char* name ) : 
  SToolBase( parent ), m_name( name ) {

  SetLogName( name );

  std::string sframe_dir(std::getenv("SFRAME_DIR"));

  // m_calib = 0;
  // m_reader = 0;
  // m_reader_up = 0;
  // m_reader_down = 0;
  
  wpCuts.clear();
  wpCuts["Loose"] = 0.605;
  wpCuts["Medium"] = 0.89;
  wpCuts["Tight"] = 0.97;
  currentWorkingPointCut = -1;

  DeclareProperty( m_name + "_Tagger",    m_tagger = "CSVv2" );
  DeclareProperty( m_name + "_WorkingPoint", m_workingPoint = "Loose" );
  DeclareProperty( m_name + "_CsvFile", m_csvFile = sframe_dir + "/../BTaggingTools/csv/CSVv2.csv" );
  DeclareProperty( m_name + "_MeasurementType_udsg", m_measurementType_udsg = "comb" );
  DeclareProperty( m_name + "_MeasurementType_bc", m_measurementType_bc = "mujets" );
  DeclareProperty( m_name + "_EffHistDirectory", m_effHistDirectory = "bTagEff" );
  DeclareProperty( m_name + "_EffFile", m_effFile = sframe_dir + "/../BTaggingTools/efficiencies/bTagEff.root" );

}

//
// destructor
//
BTaggingScaleTool::~BTaggingScaleTool() {
  // delete m_calib;
  // delete m_reader;
  // delete m_reader_up;
  // delete m_reader_down;
}

void BTaggingScaleTool::BeginInputData( const SInputData& ) throw( SError ) {

  m_logger << INFO << "Initializing BTagCalibrationStandalone" << SLogger::endmsg;
  m_logger << INFO << "CSV file:    " << m_csvFile << SLogger::endmsg;
  m_logger << INFO << "Tagger:         " << m_tagger << SLogger::endmsg;
  m_logger << INFO << "WorkingPoint: " << m_workingPoint << SLogger::endmsg;
  m_logger << INFO << "MeasurementType udsg: " << m_measurementType_udsg << SLogger::endmsg;
  m_logger << INFO << "MeasurementType bc: " << m_measurementType_bc << SLogger::endmsg;
  m_logger << INFO << "EffHistDirectory: " << m_effHistDirectory << SLogger::endmsg;
  m_logger << INFO << "Efficiency file: " << m_effFile << SLogger::endmsg;
  
  BTagEntry::OperatingPoint wp = BTagEntry::OP_LOOSE;
  if (m_workingPoint == "Loose") {
    wp = BTagEntry::OP_LOOSE;
    currentWorkingPointCut = wpCuts["Loose"];
  }
  else if (m_workingPoint == "Medium") {
    wp = BTagEntry::OP_MEDIUM;
    currentWorkingPointCut = wpCuts["Medium"];
  }
  else if (m_workingPoint == "Tight") {
    wp = BTagEntry::OP_TIGHT;
    currentWorkingPointCut = wpCuts["Tight"];
  }
  else if (m_workingPoint == "Reshaping") {
    wp = BTagEntry::OP_RESHAPING;
    currentWorkingPointCut = wpCuts["Loose"]; //placeholder
  }
  else {
    throw SError( ("Unknown working point: " + m_workingPoint).c_str(), SError::SkipCycle );
  }

  BTagCalibration m_calib(m_tagger, m_csvFile);

  m_reader.reset(new BTagCalibrationReader(wp, "central"));
  m_reader_up.reset(new BTagCalibrationReader(wp, "up"));
  m_reader_down.reset(new BTagCalibrationReader(wp, "down"));
  
  m_reader->load(m_calib, BTagEntry::FLAV_B, m_measurementType_bc);
  m_reader->load(m_calib, BTagEntry::FLAV_C, m_measurementType_bc);
  m_reader->load(m_calib, BTagEntry::FLAV_UDSG, m_measurementType_udsg);
  m_reader_up->load(m_calib, BTagEntry::FLAV_B, m_measurementType_bc);
  m_reader_up->load(m_calib, BTagEntry::FLAV_C, m_measurementType_bc);
  m_reader_up->load(m_calib, BTagEntry::FLAV_UDSG, m_measurementType_udsg);
  m_reader_down->load(m_calib, BTagEntry::FLAV_B, m_measurementType_bc);
  m_reader_down->load(m_calib, BTagEntry::FLAV_C, m_measurementType_bc);
  m_reader_down->load(m_calib, BTagEntry::FLAV_UDSG, m_measurementType_udsg);
  

  return;

}


double BTaggingScaleTool::getScaleFactor( const double& pt, const double& eta, const int& flavour, bool isTagged, const double& sigma ) {

  // Flavor
  BTagEntry::JetFlavor flavorEnum = BTagEntry::FLAV_UDSG;

  double MaxEta = 2.4;
  double abs_eta = fabs(eta);
  if (abs_eta > MaxEta) {
    // outside tracker range
    return 1.;
  }
  
  // range checking, double uncertainty if beyond
  std::pair<float, float> sf_bounds = m_reader->min_max_pt(flavorEnum, abs_eta);
  
  m_logger << DEBUG << "     flavor " << flavorEnum << " - " << sf_bounds.first << " " << sf_bounds.second << SLogger::endmsg;

  float pt_for_eval = pt;
  bool is_out_of_bounds = false;
  if (pt < sf_bounds.first) {
    pt_for_eval = sf_bounds.first + 1e-5;
    is_out_of_bounds = true;
  } else if (pt > sf_bounds.second) {
    pt_for_eval = sf_bounds.second - 0.1;
    is_out_of_bounds = true;
  }
  
  double sigmaScale = sigma;
  // double uncertainty in case jet outside normal kinematics
  if (is_out_of_bounds) {
    m_logger << DEBUG << sf_bounds.first << " - " << sf_bounds.second << SLogger::endmsg;
    m_logger << DEBUG << "out of bounds, using: " << pt_for_eval << " and " << abs_eta << SLogger::endmsg;
    sigmaScale *= 2;
  }
  
  m_logger << DEBUG << "getting scale factor " << SLogger::endmsg;
  double scalefactor = m_reader->eval(flavorEnum, eta, pt_for_eval);
  m_logger << DEBUG << "scale factor: " << scalefactor << SLogger::endmsg;
  if ((sigma > std::numeric_limits<double>::epsilon()) || (sigma < -std::numeric_limits<double>::epsilon())) {
    // m_logger << DEBUG << "limit: " << std::numeric_limits<double>::epsilon() << " value: " << sigma << SLogger::endmsg;
    if (sigma > 0) {
      double scalefactor_up =  m_reader_up->eval(flavorEnum, eta, pt_for_eval);
      scalefactor = sigmaScale*(scalefactor_up - scalefactor) + scalefactor;
    }
    else {
      double scalefactor_down =  m_reader_down->eval(flavorEnum, eta, pt_for_eval);
      scalefactor = fabs(sigmaScale)*(scalefactor_down - scalefactor) + scalefactor;
    }
  }
  if (scalefactor == 0) {
    throw SError( "Scale factor returned is zero!", SError::SkipCycle );
  }
  
  m_logger << DEBUG << "getting final weight " << flavorEnum << SLogger::endmsg;
  
  double jetweight = 1.;
  // set effMC close to one for now, need to use real value map later
  double effMC = .95;
  
  if (isTagged) {
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


double BTaggingScaleTool::getScaleFactor( const UZH::Jet& jet, const double& sigma ) {

  double jetweight = getScaleFactor(jet.pt(), jet.eta(), jet.hadronFlavour(), isTagged(jet), sigma);

  return jetweight;
  
}


double BTaggingScaleTool::getPrunedSubjetScaleFactor( const UZH::Jet& jet, const double& sigma ) {

  double jetweight = 1;
  
  for (int i = 0; i < jet.subjet_pruned_N(); ++i) {
    m_logger << DEBUG << "Looking at pruned subjet " << i
	     << ", pT=" << jet.subjet_pruned_pt()[i] << ", eta=" << jet.subjet_pruned_eta()[i]
	     << SLogger::endmsg;
    jetweight *= getScaleFactor(jet.subjet_pruned_pt()[i], jet.subjet_pruned_eta()[i], jet.subjet_pruned_hadronFlavour()[i], isTagged(jet.subjet_pruned_csv()[i]), sigma);
  }

  return jetweight;
  
}


//
// return scale for Jet collection
//
double BTaggingScaleTool::getScaleFactor( const UZH::JetVec& vJets, const double& sigma ) {

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


//
// return scale for Jet collection
//
double BTaggingScaleTool::getPrunedSubjetScaleFactor( const UZH::JetVec& vJets, const double& sigma ) {

  double scale = 1.;
  
  m_logger << DEBUG << "BTaggingScaleTool::getPrunedSubjetScaleFactor" << SLogger::endmsg;

  for (std::vector< UZH::Jet>::const_iterator itJet = vJets.begin(); itJet < vJets.end(); ++itJet) {
    m_logger << DEBUG << "Looking at jet " << itJet - vJets.begin()
	     << ", pT=" << (*itJet).pt() << ", eta=" << (*itJet).eta()
	     << SLogger::endmsg;

    scale *= getPrunedSubjetScaleFactor(*itJet, sigma);
  }  

  m_logger << DEBUG << "BTaggingScaleTool::getPrunedSubjetScaleFactor done" << SLogger::endmsg;
  return scale;

}


/// function to book histograms for efficiencies
void BTaggingScaleTool::bookHistograms() {
  
  const int nPtBins = 11;
  const int nEtaBins = 4;
  float ptBins[nPtBins+1] = {10, 20, 30, 50, 70, 100, 140, 200, 300, 670, 1000, 1500};
  float etaBins[nEtaBins+1] = {-2.5, -1.5, 0, 1.5, 2.5};
  std::vector<TString> jetCategories = {"jet", "subjet_pruned"};
  
  for (std::vector<TString>::const_iterator jetCat = jetCategories.begin(); jetCat != jetCategories.end(); ++jetCat) {
    Book( TH2F( *jetCat + "_udsg_" + m_workingPoint, *jetCat + "_udsg_" + m_workingPoint, nPtBins, ptBins, nEtaBins, etaBins ), m_effHistDirectory.c_str() );
    Book( TH2F( *jetCat + "_udsg_all", *jetCat + "_udsg_all", nPtBins, ptBins, nEtaBins, etaBins ), m_effHistDirectory.c_str() );
    Book( TH2F( *jetCat + "_c_" + m_workingPoint, *jetCat + "_c_" + m_workingPoint, nPtBins, ptBins, nEtaBins, etaBins ), m_effHistDirectory.c_str() );
    Book( TH2F( *jetCat + "_c_all", *jetCat + "_c_all", nPtBins, ptBins, nEtaBins, etaBins ), m_effHistDirectory.c_str() );
    Book( TH2F( *jetCat + "_b_" + m_workingPoint, *jetCat + "_b_" + m_workingPoint, nPtBins, ptBins, nEtaBins, etaBins ), m_effHistDirectory.c_str() );
    Book( TH2F( *jetCat + "_b_all", *jetCat + "_b_all", nPtBins, ptBins, nEtaBins, etaBins ), m_effHistDirectory.c_str() );
  }
  
}


/// function to fill jet b-tagging efficiencies
void BTaggingScaleTool::fillEfficiencies( const UZH::JetVec& vJets ) {
  
  for (std::vector< UZH::Jet>::const_iterator itJet = vJets.begin(); itJet < vJets.end(); ++itJet) {
    m_logger << DEBUG << "Looking at jet " << itJet - vJets.begin()
	     << ", pT=" << (*itJet).pt() << ", eta=" << (*itJet).eta()
	     << SLogger::endmsg;
    TString flavourString = flavourToString(itJet->hadronFlavour());
    if (isTagged(*itJet)) {
      Hist( "jet_" + flavourString + "_" + m_workingPoint, m_effHistDirectory.c_str() )->Fill( itJet->pt(), itJet->eta() );
    }
    Hist( "jet_" + flavourString + "_all", m_effHistDirectory.c_str() )->Fill( itJet->pt(), itJet->eta() );
  }
  
}

/// function to fill subjet b-tagging efficiencies
void BTaggingScaleTool::fillPrunedSubjetEfficiencies( const UZH::JetVec& vJets ) {
  
  for (std::vector< UZH::Jet>::const_iterator itJet = vJets.begin(); itJet < vJets.end(); ++itJet) {
    m_logger << DEBUG << "Looking at jet " << itJet - vJets.begin()
	     << ", pT=" << (*itJet).pt() << ", eta=" << (*itJet).eta()
	     << SLogger::endmsg;
    for (int i = 0; i < itJet->subjet_pruned_N(); ++i) {
      m_logger << DEBUG << "Looking at pruned subjet " << i
  	     << ", pT=" << itJet->subjet_pruned_pt()[i] << ", eta=" << itJet->subjet_pruned_eta()[i]
  	     << SLogger::endmsg;
      TString flavourString = flavourToString(itJet->subjet_pruned_hadronFlavour()[i]);
      if (isTagged(itJet->subjet_pruned_csv()[i])) {
        Hist( "subjet_pruned_" + flavourString + "_" + m_workingPoint, m_effHistDirectory.c_str() )->Fill( itJet->subjet_pruned_pt()[i], itJet->subjet_pruned_eta()[i] );
      }
      Hist( "subjet_pruned_" + flavourString + "_all", m_effHistDirectory.c_str() )->Fill( itJet->subjet_pruned_pt()[i], itJet->subjet_pruned_eta()[i] );
    }
  }
  
}


TString BTaggingScaleTool::flavourToString( const int& flavour ) {
  
  TString flavourString = "udsg";
  
  if (flavour == 5) {
    flavourString = "b";
  }
  else if (flavour == 4) {
    flavourString = "c";
  }
  else if (flavour == 15) {
    flavourString = "c"; // Use C scale factors for tau for now.
  }
  
  return flavourString;
  
}


bool BTaggingScaleTool::isTagged( const UZH::Jet& jet ) {
  
  if (jet.csv() > currentWorkingPointCut) {
    return true;  
  }
  return false;
  
}


bool BTaggingScaleTool::isTagged( const double& csv ) {
  
  if (csv > currentWorkingPointCut) {
    return true;  
  }
  return false;
  
}

