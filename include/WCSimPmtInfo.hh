#ifndef WCSim_PmtInfo
#define WCSim_PmtInfo

/////////////////////////////////////////////////////////////////
//
// Save the Pmt info 
//
/////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "WCSimEnumerations.hh"
class TDirectory;

class WCSimPmtInfo : public TObject {
private:
  CylinderLocation_t cylocation;
  Double_t trans_x;
  Double_t trans_y;
  Double_t trans_z;
  Double_t orien_x; 
  Double_t orien_y;
  Double_t orien_z;
  //see WCSimConstructGeometryTables.cc for scheme that works for single 20", single mPMT and hybrid.
  Int_t tube;             
  Int_t mPMT;
  Int_t mPMT_pmt;
  Float_t dark_rate;               // is individual pmt property
  Float_t conversion_factor;       // temp?
  Int_t collection_id;             // very useful to have
public: 
  WCSimPmtInfo();
  
  WCSimPmtInfo(CylinderLocation_t i, Double_t t1, Double_t t2, Double_t t3, Double_t o1, Double_t o2, Double_t o3, Int_t tubeid);

  WCSimPmtInfo(CylinderLocation_t i, Double_t t1, Double_t t2, Double_t t3, Double_t o1, Double_t o2, Double_t o3, Int_t tubeid,  Int_t mPMTid, Int_t mPMT_pmtid);
  
  virtual ~WCSimPmtInfo();

  CylinderLocation_t Get_cylocation(){return cylocation;}
  Double_t Get_transx() {return trans_x;}
  Double_t Get_transy() {return trans_y;}
  Double_t Get_transz() {return trans_z;}
  Double_t Get_orienx() {return orien_x;}
  Double_t Get_orieny() {return orien_y;}
  Double_t Get_orienz() {return orien_z;}
  Int_t Get_tubeid(){return tube;}
  Int_t Get_mPMTid(){return mPMT;}
  Int_t Get_mPMT_pmtid(){return mPMT_pmt;}
  Float_t Get_darkrate(){return dark_rate;}
  Float_t Get_conv_factor(){return conversion_factor;}

  void Set_Darkrate(Float_t rate, Float_t conversion){
    dark_rate = rate;
    conversion_factor = conversion;
    return;
  }
  void Set_collectionID(Int_t id){collection_id = id;}
  Int_t Get_collectionID(){return collection_id;}

  ClassDef(WCSimPmtInfo,2);
};

#endif
