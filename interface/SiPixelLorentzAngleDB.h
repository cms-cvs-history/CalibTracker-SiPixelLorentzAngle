#ifndef CalibTracker_SiPixelLorentzAngleDB_SiPixelLorentzAngleDB_h
#define CalibTracker_SiPixelLorentzAngleDB_SiPixelLorentzAngleDB_h

#include <map>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/DetId/interface/DetId.h"

// #include "CalibTracker/SiStripLorentzAngle/interface/SiStripLorentzAngleAlgorithm.h"

class SiPixelLorentzAngleDB : public edm::EDAnalyzer
{
 public:
  
  explicit SiPixelLorentzAngleDB(const edm::ParameterSet& conf);
  
  virtual ~SiPixelLorentzAngleDB();
  
  virtual void beginJob(const edm::EventSetup& c);
  
  virtual void endJob(); 
  
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
  
 private:

  std::vector< std::pair<uint32_t, float> > detid_la;
  edm::ParameterSet conf_;
//   double appliedVoltage_;
//   double chargeMobility_;
//   double temperature_;
//   double temperatureerror_;
//   double rhall_;
//   double holeBeta_;
//   double holeSaturationVelocity_;
//   SiStripLorentzAngleAlgorithm *siStripLorentzAngleAlgorithm_;
};


#endif