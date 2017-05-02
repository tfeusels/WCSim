#ifndef WCSimRunAction_h
#define WCSimRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4String.hh"

#include "TFile.h"
#include "TTree.h"
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"
#include "WCSimDetectorConstruction.hh"
#include "WCSimEnumerations.hh"        
#include "evNtuple.h"
#include "WCSimRandomParameters.hh"


class G4Run;
class WCSimRunActionMessenger;

class WCSimRunAction : public G4UserRunAction
{
public:
  WCSimRunAction(WCSimDetectorConstruction*, WCSimRandomParameters*);
  ~WCSimRunAction();

public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
  void SetRootFileName(G4String fname) { RootFileName = fname; }
  G4String GetRootFileName() { return RootFileName; }
  void SetOptionalRootFile(G4bool choice) { useDefaultROOTout = choice; }
  G4bool GetRootFileOption() { return useDefaultROOTout; }
  void FillGeoTree();
  TTree* GetTree(){return WCSimTree;}
  TTree* GetFlatMasterTree(){return masterTree;}
  TTree* GetGeoTree(){return geoTree;}
  TTree* GetOptionsTree(){return optionsTree;}
  WCSimRootGeom* GetRootGeom(){return wcsimrootgeom;}
  WCSimRootEvent* GetRootEvent(){return wcsimrootsuperevent;}
  WCSimRootOptions* GetRootOptions(){return wcsimrootoptions;}

  //Need to share with EventAction:
  TTree *GetCherenkovHitsTree(){return cherenkovHitsTree;}
  TTree *GetCherenkovDigiHitsTree(){return cherenkovDigiHitsTree;}
  TTree *GetTracksTree(){return tracksTree;}
  TTree *GetTriggerTree(){return triggerTree;}
  TTree *GetEventInfoTree(){return eventInfoTree;}

  void SetEventHeaderNew(G4int run_id, G4int event_id, G4int subevent_id){
    event = event_id;
    run = run_id;
    subevent = subevent_id;}

  void SetTriggerInfoNew(TriggerType_t trigger_type,int trigger_info, /* std::vector<G4float> trigger_info, */
			 G4double trigger_length, G4double trigger_start){
    trig_type = trigger_type;
    /*
    for(size_t j = 0; j < trigger_info.size(); j++){
      if(j < 10)
	trig_info[j] = trigger_info[j];
      else{
	G4cerr << "Trigger Info too large for array" << G4endl;
	break;
      }
      }*/
    trig_info = trigger_info;
    trig_length = trigger_length;
    trig_start = trigger_start;
    return;
  }

  eventNtuple * GetMyStruct(){return evNtup;}

  void SetTree(TTree* tree){WCSimTree=tree;}
  void SetGeoTree(TTree* tree){geoTree=tree;}
  void SetRootEvent(WCSimRootEvent* revent){wcsimrootsuperevent=revent;}
  void SetRootGeom(WCSimRootGeom* rgeom){wcsimrootgeom=rgeom;}

  int  GetNumberOfEventsGenerated() { return numberOfEventsGenerated;}
  int  GetNtuples(){return ntuples;}

  void incrementEventsGenerated() { numberOfEventsGenerated++;}
  void incrementWaterTubeHits()   { numberOfTimesWaterTubeHit++;} 
  void incrementFVWaterTubeHits() { numberOfTimesFVWaterTubeHit++;} 
  void incrementCatcherHits()     { numberOfTimesCatcherHit++;}
  void SetNtuples(int ntup) {ntuples=ntup;}
  
  //New:
  void FillFlatGeoTree();

private:
  // MFechner : set by the messenger
  std::string RootFileName;
  // Only required for verification scripts and current fiTQun tuning
  // But making initialization very slow due to large TCloneArray init.
  G4bool useDefaultROOTout;

  //
  TTree* WCSimTree;
  TTree* geoTree;
  TTree* optionsTree;
  WCSimRootEvent* wcsimrootsuperevent;
  WCSimRootGeom* wcsimrootgeom;
  WCSimRootOptions* wcsimrootoptions;
  WCSimDetectorConstruction* wcsimdetector;
  WCSimRandomParameters* wcsimrandomparameters;

  int numberOfEventsGenerated;
  int numberOfTimesWaterTubeHit;
  int numberOfTimesFVWaterTubeHit;
  int numberOfTimesCatcherHit;

  WCSimRunActionMessenger* messenger;
  int ntuples;  // 1 for ntuples to be written

  //new: vars for FLAT Trees
  TTree* masterTree;
  TTree* geomTree;
  TTree* cherenkovHitsTree;           //cfr. RecoHitSeriesMap
  TTree* cherenkovDigiHitsTree;       //cfr. RecoHitSeriesMap
  TTree* tracksTree;
  TTree* triggerTree;
  TTree* eventInfoTree;
  TTree* optionsTree_flat;

  // GeoTree
  Char_t geo_type_string[20];   
  double cyl_radius;     
  double cyl_length;     
  Char_t pmt_id0_string[50];
  Char_t pmt_id1_string[50];
  double pmt_radius_id[2];
  Char_t pmt_od_string[50];
  double pmt_radius_od;
  int numPMT_id_tot;
  int numPMT_id[2];
  int numPMT_od;
  int orient; 
  double offset_x;
  double offset_y;
  double offset_z;
  int num_mPMT[2];           // in theory: hybrid 'can' be two types of mPMT, but not urgent.
  int tube_id[5000000];
  int mPMT_id[5000000];
  int mPMT_pmt_id[5000000];
  int collection_id[5000000];
  double tube_x[5000000];
  double tube_y[5000000];
  double tube_z[5000000];
  double dir_x[5000000];
  double dir_y[5000000];
  double dir_z[5000000];
  double theta[5000000];
  double phi[5000000];
  
  CylinderLocation_t cylLocation[5000000];
  
  //EventHeader
  int run;
  int event;
  int subevent;

  //TriggerTree
  TriggerType_t trig_type;
  int trig_info;              //double trig_info[100];
  double trig_length;
  double trig_start;
  
  //Event info: General, Tracks and Hits
  eventNtuple *evNtup;
    
};

#endif
