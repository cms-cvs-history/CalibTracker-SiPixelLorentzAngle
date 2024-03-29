#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CalibTracker/SiPixelLorentzAngle/interface/SiPixelLorentzAngleDB.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

using namespace std;
using namespace edm;

  //Constructor

SiPixelLorentzAngleDB::SiPixelLorentzAngleDB(edm::ParameterSet const& conf) : 
  conf_(conf){
  	magneticField_ = conf_.getParameter<double>("magneticField");
	recordName_ = conf_.getUntrackedParameter<std::string>("record","SiPixelLorentzAngleRcd");
	useFile_ = conf_.getParameter<bool>("useFile");		
	fileName_ = conf_.getParameter<string>("fileName");

        BPixParameters_ = conf_.getUntrackedParameter<Parameters>("BPixParameters");
        FPixParameters_ = conf_.getUntrackedParameter<Parameters>("FPixParameters");
}

  //BeginJob

void SiPixelLorentzAngleDB::beginJob(){
  
}
// Virtual destructor needed.

SiPixelLorentzAngleDB::~SiPixelLorentzAngleDB() {  

}  

// Analyzer: Functions that gets called by framework every event

void SiPixelLorentzAngleDB::analyze(const edm::Event& e, const edm::EventSetup& es)
{

	SiPixelLorentzAngle* LorentzAngle = new SiPixelLorentzAngle();
	   
        int PBDetUnitSize = 0;
        int PFDetUnitSize = 0;
	
	edm::ESHandle<TrackerGeometry> pDD;
	es.get<TrackerDigiGeometryRecord>().get( pDD );
	edm::LogInfo("SiPixelLorentzAngle") <<" There are "<<pDD->detUnits().size() <<" detectors"<<std::endl;


//temporary solution to put the same LA to the given list of det ids
        if(useFile_){

           unsigned int numofrawids = 0;
           unsigned int rawid, aa, bb;
           float la = 0.106;

	   std::ifstream in_file;  // data file pointer
           in_file.open( fileName_.c_str() , std::ios::in ); // in C++
	   if (in_file.bad()) {
      	      edm::LogError("SiPixelLorentzAngleDB")<<"Input file not found"<<std::endl;
  	   } 
           else {

              in_file >> rawid >> aa >> bb;
              while(!in_file.eof()){
                 numofrawids++;
                 cout << rawid << "   " << la << endl;
                 LorentzAngle->putLorentzAngle(rawid,la);
                 in_file >> rawid >> aa >> bb;
              }

              cout << "    DetId input file found with " << numofrawids << " detids " << endl;
           }

        }
//end of temporary solution

	
	for(TrackerGeometry::DetUnitContainer::const_iterator it = pDD->detUnits().begin(); it != pDD->detUnits().end(); it++){
    
	   if( dynamic_cast<PixelGeomDetUnit*>((*it))!=0){
		DetId detid=(*it)->geographicalId();
			
		// fill bpix values for LA 
		if(detid.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) {
				
                   PBDetUnitSize++;
                   PXBDetId pxdetid = PXBDetId(detid);

		   if(!useFile_){
                        cout << " hp:barrel:" << "  layer=" << pxdetid.layer() << "  ladder=" << pxdetid.ladder() << "  module=" << pxdetid.module() << endl;
/*hp
			if ( ! LorentzAngle->putLorentzAngle(detid.rawId(),bPixLorentzAnglePerTesla_) )
			edm::LogError("SiPixelLorentzAngleDB")<<"[SiPixelLorentzAngleDB::analyze] detid already exists"<<std::endl;
*/ //hp

		        for(Parameters::iterator it = BPixParameters_.begin(); it != BPixParameters_.end(); ++it) {
//                           cout << " PSet: " << *it << ", module = " << it->getParameter<unsigned int>("module") << endl;
                           if( it->getParameter<unsigned int>("module") == pxdetid.module() && it->getParameter<unsigned int>("layer") == pxdetid.layer() )
                           {
                              float lorentzangle = (float)it->getParameter<double>("angle");
                              LorentzAngle->putLorentzAngle(detid.rawId(),lorentzangle);
                           }
                        }

		   } else {
			//cout << "    method for reading file implementing (hp)" << endl;
		   }
			
	        // fill fpix values for LA 
		} else if(detid.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap)) {
				
                   PFDetUnitSize++;

                   PXFDetId pxdetid = PXFDetId(detid);

		   if(!useFile_){
                         cout << " hp:endcap:" << "  side=" << pxdetid.side() << "  disk=" << pxdetid.disk() << "  blade=" << pxdetid.blade() << "  panel=" << pxdetid.panel() << "  module=" << pxdetid.module() << endl;
/*hp
		         if ( ! LorentzAngle->putLorentzAngle(detid.rawId(),fPixLorentzAnglePerTesla_) )  edm::LogError("SiPixelLorentzAngleDB")<<"[SiPixelLorentzAngleDB::analyze] detid already exists"<<std::endl;
*/ //hp

		        for(Parameters::iterator it = FPixParameters_.begin(); it != FPixParameters_.end(); ++it) {
//                           cout << " PSet: " << *it << ", module = " << it->getParameter<unsigned int>("module") << endl;
                           if( it->getParameter<unsigned int>("side") == pxdetid.side() && it->getParameter<unsigned int>("disk") == pxdetid.disk() && it->getParameter<unsigned int>("HVgroup") == HVgroup(pxdetid.panel(),pxdetid.module()) )
                           {
                              float lorentzangle = (float)it->getParameter<double>("angle");
                              LorentzAngle->putLorentzAngle(detid.rawId(),lorentzangle);
                           }
                        }

		   } else {
			//cout << "    method for reading file implementing (hp),  rawId = " << detid.rawId() << endl;
                   }

                // neither fpix nor bpix in Pixel
		} else {
			edm::LogError("SiPixelLorentzAngleDB")<<"[SiPixelLorentzAngleDB::analyze] detid is Pixel but neither bpix nor fpix"<<std::endl;
		}
		
	   }
			
	}      
  	
        cout << "  Total number of pixel detector units: " << PBDetUnitSize + PFDetUnitSize << " (barrel: " << PBDetUnitSize << ", forward: " << PFDetUnitSize << ")" << endl;

	edm::Service<cond::service::PoolDBOutputService> mydbservice;
	if( mydbservice.isAvailable() ){
		try{
			if( mydbservice->isNewTagRequest(recordName_) ){
				mydbservice->createNewIOV<SiPixelLorentzAngle>(LorentzAngle,
									       mydbservice->beginOfTime(),
									       mydbservice->endOfTime(),
									       recordName_);
			} else {
				mydbservice->appendSinceTime<SiPixelLorentzAngle>(LorentzAngle,
										  mydbservice->currentTime(),
										  recordName_);
			}
		}catch(const cond::Exception& er){
			edm::LogError("SiPixelLorentzAngleDB")<<er.what()<<std::endl;
		}catch(const std::exception& er){
			edm::LogError("SiPixelLorentzAngleDB")<<"caught std::exception "<<er.what()<<std::endl;
		}catch(...){
			edm::LogError("SiPixelLorentzAngleDB")<<"Funny error"<<std::endl;
		}
	}else{
		edm::LogError("SiPixelLorentzAngleDB")<<"Service is unavailable"<<std::endl;
	}
   

}


unsigned int SiPixelLorentzAngleDB::HVgroup(unsigned int panel, unsigned int module){

   //edm::LogError("SiPixelLorentzAngleDB") << "HVgroup yet to be adjusted to the upgrade" <<std::endl;

/*
   if( 1 == panel && ( 1 == module || 2 == module ))  {
      return 1;
   }
   else if( 1 == panel && ( 3 == module || 4 == module ))  {
      return 2;
   }
   else if( 2 == panel && 1 == module )  {
      return 1;
   }
   else if( 2 == panel && ( 2 == module || 3 == module ))  {
      return 2;
   }
   else {
      // cout << " *** error *** in SiPixelLorentzAngleDB::HVgroup(...), panel = " << panel << ", module = " << module << endl;
      edm::LogError("SiPixelLorentzAngleDB") << "HVgroup panel = " << panel << ", module = " << module << endl;
      return 0;
   }
*/
   return 1;
   
}



void SiPixelLorentzAngleDB::endJob(){


}


