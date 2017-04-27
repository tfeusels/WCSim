#include "WCSimStackingAction.hh"
#include "WCSimDetectorConstruction.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//class WCSimDetectorConstruction;

WCSimStackingAction::WCSimStackingAction(WCSimDetectorConstruction* myDet):DetConstruct(myDet) {
  nWarnings = 0.;

}
WCSimStackingAction::~WCSimStackingAction(){;}


G4ClassificationOfNewTrack WCSimStackingAction::ClassifyNewTrack
(const G4Track* aTrack) 
{
  G4ClassificationOfNewTrack classification    = fWaiting;
  G4ParticleDefinition*      particleType      = aTrack->GetDefinition();
  

  // Make sure it is an optical photon
  if( particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
    {
      // TF: cleaned this up a little: no repetition of code.
      // also don't know why CreatorProcess() == NULL needs to have QE applied.
      // use QE method for ALL.
      if( aTrack->GetCreatorProcess() == NULL ||          // eg. particle gun/gps photons
	  ( aTrack->GetCreatorProcess() != NULL && 
	    ((G4VProcess*)(aTrack->GetCreatorProcess()))->GetProcessType() != fOptical) ) {
	
	G4float photonWavelength = (2.0*M_PI*197.3)/(aTrack->GetTotalEnergy()/eV);
	// MF : translated from skdetsim : better to increase the number of photons
	// than to throw in a global factor  at Digitization time !
	G4float ratio = 1./(1.0-0.25); 
	G4float wavelengthQE = 0;
	if(nWarnings < 10){
	  G4cout << " WARNING: skdetsim fudge factor being applied to QE : " << ratio << G4endl;
	  G4cout << " ====================================================================== " << G4endl;
	  nWarnings++;
	}
	// XQ: get the maximum QE and multiply it by the ratio
	// only work for the range between 240 nm and 660 nm for now 
	// Even with WLS

	G4String * WCIDCollectionNames = DetConstruct->GetIDCollectionNames();
	G4int IDnoTypes = DetConstruct->GetNoIDtypes();

	if (DetConstruct->GetPMT_QE_Method()==1){
	  if(IDnoTypes == 1)
	    wavelengthQE  = DetConstruct->GetPMTQE(WCIDCollectionNames[0],photonWavelength,1,240,700,ratio);
	  else if(IDnoTypes == 2){
	    // For stacking action: need to take the max of both QE curves and then rescale at PMT level (WCSD) 
	    // depending  on ID PMT type
	    G4double wavelengthQE0  = DetConstruct->GetPMTQE(WCIDCollectionNames[0],photonWavelength,1,240,700,ratio);
	    G4double wavelengthQE1  = DetConstruct->GetPMTQE(WCIDCollectionNames[1],photonWavelength,1,240,700,ratio);
	    wavelengthQE = std::max(wavelengthQE0, wavelengthQE1);
	  }
	}else if (DetConstruct->GetPMT_QE_Method()==2){
	  if(IDnoTypes == 1)
	    wavelengthQE  = DetConstruct->GetPMTQE(WCIDCollectionNames[0],photonWavelength,0,240,700,ratio);
	  else if(IDnoTypes == 2){
	    G4double wavelengthQE0 = DetConstruct->GetPMTQE(WCIDCollectionNames[0],photonWavelength,0,240,700,ratio);
	    G4double wavelengthQE1 = DetConstruct->GetPMTQE(WCIDCollectionNames[1],photonWavelength,0,240,700,ratio);
	    wavelengthQE = std::max(wavelengthQE0,wavelengthQE1);
	  }
	}else if (DetConstruct->GetPMT_QE_Method()==3 || DetConstruct->GetPMT_QE_Method() == 4){
	  wavelengthQE = 1.1;
	}
	
	if( G4UniformRand() > wavelengthQE )
	  classification = fKill;
      }
    }
  
  return classification;
}

void WCSimStackingAction::NewStage() {;}
void WCSimStackingAction::PrepareNewEvent() {;}

