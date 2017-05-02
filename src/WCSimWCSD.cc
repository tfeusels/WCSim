#include "WCSimWCSD.hh"
#include "G4ParticleTypes.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4ios.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>

#include "WCSimDetectorConstruction.hh"
#include "WCSimTrackInformation.hh"

#include "WCSimSteppingAction.hh"

WCSimWCSD::WCSimWCSD(G4String CollectionName, G4String name,WCSimDetectorConstruction* myDet)
:G4VSensitiveDetector(name)
{
  // Place the name of this collection on the list.  We can have more than one
  // in principle.  CollectionName is a vector.

  // Note there is some sort of problem here.  If I use the name
  // Which has a "/" in it, I can find this collection later using 
  // GetCollectionID()

  collectionName.insert(CollectionName);
  
  fdet = myDet;
  
  HCID = -1;
}

WCSimWCSD::~WCSimWCSD() {}

void WCSimWCSD::Initialize(G4HCofThisEvent* HCE)
{
  // Make a new hits collection. With the name we set in the constructor
  hitsCollection = new WCSimWCHitsCollection
    (SensitiveDetectorName,collectionName[0]);

  // This is a trick.  We only want to do this once.  When the program
  // starts HCID will equal -1.  Then it will be set to the pointer to
  // this collection.

  
  // Get the Id of the "0th" collection
  if (HCID<0){
    HCID =  GetCollectionID(0); 
  }  
  // Add it to the Hit collection of this event.
  HCE->AddHitsCollection( HCID, hitsCollection );  

  // So far collectionName[0] == SensitiveDetectorName
  // keep track of situations where this is not the case
  if(collectionName[0] != SensitiveDetectorName){
    G4cout << "For collectionID  : " << GetCollectionID(0) << G4endl;
    G4cout << "collName : " << collectionName[0]  << ", and SensName : " << SensitiveDetectorName << G4endl;
    G4cout << "Make sure behaviour is as expected ! " << G4endl;
  }
  // Initialize the Hit map to all tubes not hit.
  PMTHitMap.clear();
  // Trick to access the static maxPE variable.  This will go away with the 
  // variable.

  WCSimWCHit* newHit = new WCSimWCHit();
  newHit->SetMaxPe(0);
  delete newHit;
}

G4bool WCSimWCSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{ 
  G4StepPoint*       preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle  theTouchable = preStepPoint->GetTouchableHandle();
  G4VPhysicalVolume* thePhysical  = theTouchable->GetVolume();


  //XQ 3/30/11 try to get the local position try to add the position and direction
  G4ThreeVector worldPosition = preStepPoint->GetPosition();
  G4ThreeVector localPosition = theTouchable->GetHistory()->GetTopTransform().TransformPoint(worldPosition);
  G4ThreeVector worldDirection = preStepPoint->GetMomentumDirection();
  G4ThreeVector localDirection = theTouchable->GetHistory()->GetTopTransform().TransformAxis(worldDirection);


  WCSimTrackInformation* trackinfo 
    = (WCSimTrackInformation*)(aStep->GetTrack()->GetUserInformation());
  
  G4int primParentID = -1;
  if (trackinfo)
    //Skip secondaries and match to mother process, eg. muon, decay particle, gamma from pi0/nCapture.
    primParentID = trackinfo->GetPrimaryParentID();  //!= ParentID. 
  else // if there is no trackinfo, then it is a primary particle!
    primParentID = aStep->GetTrack()->GetTrackID();    


  G4int    trackID           = aStep->GetTrack()->GetTrackID();
  G4String volumeName        = aStep->GetTrack()->GetVolume()->GetName();
  
  //XQ Add the wavelength there
  G4float  wavelength = (2.0*M_PI*197.3)/( aStep->GetTrack()->GetTotalEnergy()/eV);
  
  G4double energyDeposition  = aStep->GetTotalEnergyDeposit();
  G4double hitTime           = aStep->GetPreStepPoint()->GetGlobalTime();

  G4ParticleDefinition *particleDefinition = 
    aStep->GetTrack()->GetDefinition();

  if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition() 
       && energyDeposition == 0.0) 
    return false;
  // MF : I don't see why other particles should register hits
  // they don't in skdetsim. 
  if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition())
    return false;
 

  // M Fechner : too verbose
  //  if (aStep->GetTrack()->GetTrackStatus() == fAlive)G4cout << "status is fAlive\n";
  if ((aStep->GetTrack()->GetTrackStatus() == fAlive )                                         //// TF: If absorbed in postVol, then it's killed, so will be Hit, if it's alive, then it's not a hit.
      &&(particleDefinition == G4OpticalPhoton::OpticalPhotonDefinition()))
    return false;

  // TF: Problem: photons can go through the sensitive detector (glass)
  // and be killed in the next volume by Absorption, eg. in the BlackSheet, but will then
  // still be counted as hits here. So, either recode as extended/optical/LXe cathode implementation
  // or easier: make sure they cross the cathode, as the logical Boundary "cathode" can't be a
  // sensitive detector (I think).
  
  G4StepPoint        *postStepPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume  *postVol = postStepPoint->GetPhysicalVolume();

  //if (thePhysical)  G4cout << " thePrePV:  " << thePhysical->GetName()  << G4endl;
  //if (postVol) G4cout << " thePostPV: " << postVol->GetName() << G4endl;

  // For the record: volumeName == thePhysical->GetName() 

  //Optical Photon must pass through glass into PMT interior!
  // What about the other way around? TF: current interior won't keep photons alive like in reality
  // Not an issue yet, because then interior needs to be a sensitive detector, when postStepPoint is the glass.
  if(postVol->GetName() != "InteriorWCPMT")
    return false;





  //  if ( particleDefinition ==  G4OpticalPhoton::OpticalPhotonDefinition() ) 
  // G4cout << volumeName << " hit by optical Photon! " << G4endl;
    
  // Make the tubeTag string based on the replica numbers
  // See WCSimDetectorConstruction::DescribeAndRegisterPMT() for matching
  // tag construction.

  std::stringstream tubeTag;

  // Start tubeTag with mother to distinguish different PMT hierarchies
//  G4LogicalVolume *theMother = thePhysical->GetMotherLogical();
//  if (theMother != NULL)
//    tubeTag << theMother->GetName() << ":";

//  tubeTag << thePhysical->GetName(); 
  for (G4int i = theTouchable->GetHistoryDepth()-1 ; i >= 0; i--){
    tubeTag << ":" << theTouchable->GetVolume(i)->GetName();
    tubeTag << "-" << theTouchable->GetCopyNumber(i);
  }
  //  tubeTag << ":" << theTouchable->GetVolume(i)->GetCopyNo(); 

  // Get the tube ID from the tubeTag
  G4int replicaNumber = WCSimDetectorConstruction::GetTubeID(tubeTag.str());

  G4String * WCIDCollectionNames = fdet->GetIDCollectionNames();
  G4int IDnoTypes = fdet->GetNoIDtypes();
  G4int IDno = -1;
  // Figure out which ID PMT
  if(tubeTag.str().find(WCIDCollectionNames[0]) != std::string::npos)
    IDno = 0;
  else {
    IDno = 1;
    if(tubeTag.str().find(WCIDCollectionNames[1]) == std::string::npos)
      G4cerr << "DID not find collection name in " << tubeTag.str() << G4endl;
  }
  G4float theta_angle = 0.;
  G4float effectiveAngularEfficiency = 0.;


  
  G4float ratio = 1.;
  G4float maxQE = 0.;
  G4float photonQE = 0.;
  if (fdet->GetPMT_QE_Method() == 4){
    photonQE = 1.1;
  }else if (fdet->GetPMT_QE_Method() == 1){
    if(IDnoTypes == 1)
      photonQE = 1.1;
    else if(IDnoTypes == 2){
      //// For hybrid ID : scaled in stacking action by max(QE1, QE2, lambda)
      // here rescale based on which PMT
      
      //grab QE_lambda1, QE_lambda2, check which is max and rescale appropriately.
      ratio = 1./(1.-0.25);
      G4double QE_lambda0 = fdet->GetPMTQE(WCIDCollectionNames[0],wavelength,1,240,700, ratio);
      G4double QE_lambda1 = fdet->GetPMTQE(WCIDCollectionNames[1],wavelength,1,240,700, ratio);
      if(QE_lambda0 >= QE_lambda1){
	if(IDno == 0)
	  photonQE = 1.1;
	else if(IDno == 1){
	  photonQE = QE_lambda1/QE_lambda0;
	}
      } else if(QE_lambda0 < QE_lambda1){
	if(IDno == 0)
	  photonQE = QE_lambda0/QE_lambda1;
	else if(IDno == 1){
	  photonQE = 1.1;
	}
      }
    }
  }else if (fdet->GetPMT_QE_Method()==2){
    if(IDnoTypes == 1){
      maxQE = fdet->GetPMTQE(WCIDCollectionNames[IDno],wavelength,0,240,700,ratio);
      photonQE = fdet->GetPMTQE(volumeName, wavelength,1,240,700,ratio);
      photonQE = photonQE/maxQE;
    }else if(IDnoTypes == 2){
      G4double maxQE0 = fdet->GetPMTQE(WCIDCollectionNames[0],wavelength,0,240,700,ratio);
      G4double maxQE1 = fdet->GetPMTQE(WCIDCollectionNames[1],wavelength,0,240,700,ratio);
      photonQE = fdet->GetPMTQE(volumeName, wavelength,1,240,700,ratio);
      photonQE = photonQE/std::max(maxQE0,maxQE1);
    }
  }else if (fdet->GetPMT_QE_Method()==3){
    ratio = 1./(1.-0.25);
    photonQE = fdet->GetPMTQE(volumeName, wavelength,1,240,700,ratio);
  }
  
  

  if (G4UniformRand() <= photonQE){
    
     G4double local_x = localPosition.x();
     G4double local_y = localPosition.y();
     G4double local_z = localPosition.z();
     theta_angle = acos(fabs(local_z)/sqrt(pow(local_x,2)+pow(local_y,2)+pow(local_z,2)))/3.1415926*180.;
     effectiveAngularEfficiency = fdet->GetPMTCollectionEfficiency(theta_angle, volumeName);
     if (G4UniformRand() <= effectiveAngularEfficiency || fdet->UsePMT_Coll_Eff()==0){
       //Retrieve the pointer to the appropriate hit collection. 
       //Since volumeName is the same as the SD name, this works. 
       G4SDManager* SDman = G4SDManager::GetSDMpointer();
       G4RunManager* Runman = G4RunManager::GetRunManager();
       G4int collectionID = SDman->GetCollectionID(volumeName);
       const G4Event* currentEvent = Runman->GetCurrentEvent();
       G4HCofThisEvent* HCofEvent = currentEvent->GetHCofThisEvent();
       hitsCollection = (WCSimWCHitsCollection*)(HCofEvent->GetHC(collectionID));
      
       // note: for hybrid ID: different PMTs, are registered as different collecionID and
       // hence have their own hitsCollection.
       // Below will add hits to that hitsCollection (different for each ID PMT) and
       // register the hit PMT in the PMTHitMap, using tubeID 'replicaNumber' (for all PMTs of both ID).
       // That is just to monitor which PMTs have been hit.

       // If this tube hasn't been hit add it to the collection	 
       if (this->PMTHitMap[replicaNumber] == 0)
       //if (PMTHitMap.find(replicaNumber) == PMTHitMap.end())  TF attempt to fix
	 {
	   WCSimWCHit* newHit = new WCSimWCHit();
	   newHit->SetTubeID(replicaNumber);
	   newHit->SetTrackID(trackID);
	   newHit->SetEdep(energyDeposition); 
	   newHit->SetLogicalVolume(thePhysical->GetLogicalVolume());
	   
	   G4AffineTransform aTrans = theTouchable->GetHistory()->GetTopTransform();
	   newHit->SetRot(aTrans.NetRotation());
	   
	   aTrans.Invert();
	   newHit->SetPos(aTrans.NetTranslation());
	   
	   // Set the hitMap value to the collection hit number
	   PMTHitMap[replicaNumber] = hitsCollection->insert( newHit );
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPe(hitTime);
	   (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddParentID(primParentID);
	   
	   //     if ( particleDefinition != G4OpticalPhoton::OpticalPhotonDefinition() )
	   //       newHit->Print();
	     
	 }
       else {
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddPe(hitTime);
	 (*hitsCollection)[PMTHitMap[replicaNumber]-1]->AddParentID(primParentID);
	 
       }
     }
  }

  return true;
}

void WCSimWCSD::EndOfEvent(G4HCofThisEvent* HCE)
{

 
  if (verboseLevel>0) 
  { 
    //Need to specify which collection in case multiple geometries are built:
    G4String WCIDCollectionName = fdet->GetIDCollectionName(0);
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4int collectionID = SDman->GetCollectionID(WCIDCollectionName);
    hitsCollection = (WCSimWCHitsCollection*)HCE->GetHC(collectionID);
    
    G4int numHits = hitsCollection->entries();

    G4cout << "======================================================================" << G4endl;
    G4cout << "There are " << numHits << " tubes hit in the : " << WCIDCollectionName << G4endl;
    for (G4int i=0; i < numHits; i++) 
      (*hitsCollection)[i]->Print();
    
    G4int IDnoTypes = fdet->GetNoIDtypes();
    if(IDnoTypes == 2){
      WCIDCollectionName = fdet->GetIDCollectionName(1);
      collectionID = SDman->GetCollectionID(WCIDCollectionName);
      hitsCollection = (WCSimWCHitsCollection*)HCE->GetHC(collectionID);
    
      numHits = hitsCollection->entries();
      G4cout << "======================================================================" << G4endl;
      G4cout << "And there are " << numHits << " tubes hit in the : " << WCIDCollectionName << G4endl;
      for (G4int i=0; i < numHits; i++) 
	(*hitsCollection)[i]->Print();
    }
    

      /*
    {
      if(abs((*hitsCollection)[i]->GetTubeID() - 1584)  < 5){
	  std::cout << (*hitsCollection)[i]->GetTubeID() << std::endl;
	  (*hitsCollection)[i]->Print();
      }
      }*/

    /* Detailed debug:
    G4cout << "Through mPMTLV " << WCSimSteppingAction::n_photons_through_mPMTLV << G4endl;
    G4cout << "Through Acrylic " << WCSimSteppingAction::n_photons_through_acrylic << G4endl;
    G4cout << "Through Gel " << WCSimSteppingAction::n_photons_through_gel << G4endl;
    G4cout << "On Blacksheet " << WCSimSteppingAction::n_photons_on_blacksheet << G4endl;
    G4cout << "On small PMT " << WCSimSteppingAction::n_photons_on_smallPMT << G4endl;
    */
  } 
}

