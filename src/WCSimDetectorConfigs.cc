//  -*- mode:c++; tab-width:4;  -*-
#include "WCSimDetectorConstruction.hh"
#include "WCSimPMTObject.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

/***********************************************************
 *
 * This file contains the setup functions for various 
 * detector configurations.  These can be set up by 
 * default in WCSimDetectorConstruction.cc or called
 * in mac files by adding them to WCSimDetectorMessenger.cc.
 *
 * Sourcefile for the WCSimDetectorConstruction class
 *
 ***********************************************************/



void WCSimDetectorConstruction::SetSuperKGeometry()
{
  WCDetectorName = "SuperK";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("PMT20inch", WCIDCollectionName[0]);
  WCPMTName[0] = PMT->GetPMTName();
  WCPMTExposeHeight[0] = PMT->GetExposeHeight();
  WCPMTRadius[0] = PMT->GetRadius();

  WCIDDiameter          = 33.6815*m; //16.900*2*cos(2*pi*rad/75)*m; //inner detector diameter
  WCIDHeight            = 36.200*m; //"" "" height
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCBarrelNumPMTHorizontal  = 150; 
  WCBarrelNRings        = 17.;
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3; 
  WCCapPMTSpacing       = 0.707*m; // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}


// Note: the actual coverage is 20.27%
void WCSimDetectorConstruction::SuperK_20inchPMT_20perCent()
{
  WCDetectorName = "SuperK_20inchPMT_20perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("PMT20inch", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 33.6815*m; //16.900*2*cos(2*pi*rad/75)*m; //inner detector diameter
  WCIDHeight            = 36.200*m; //"" "" height
  WCBarrelPMTOffset     = 0.0715*m; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3; 
  WCPMTPercentCoverage  = 20.27;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}


// Note: the actual coverage is 20.27%
void WCSimDetectorConstruction::SuperK_20inchBandL_20perCent()
{
  WCDetectorName = "SuperK_20inchBandL_20perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine20inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 33.6815*m; //16.900*2*cos(2*pi*rad/75)*m; //inner detector diameter
  WCIDHeight            = 36.200*m; //"" "" height
  WCBarrelPMTOffset     = 0.0715*m; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3; 
  WCPMTPercentCoverage  = 20.27;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}


// Note: the actual coverage is 14.59%
void WCSimDetectorConstruction::SuperK_12inchBandL_15perCent()
{
  WCDetectorName = "SuperK_12inchBandL_15perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine12inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 33.6815*m; //16.900*2*cos(2*pi*rad/75)*m; //inner detector diameter
  WCIDHeight            = 36.200*m; //"" "" height
  WCBarrelPMTOffset     = 0.0715*m; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 14.59;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}


// Note: the actual coverage is 13.51%
void WCSimDetectorConstruction::SuperK_20inchBandL_14perCent()
{
  WCDetectorName = "SuperK_20inchBandL_14perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine20inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 33.6815*m; //16.900*2*cos(2*pi*rad/75)*m; //inner detector diameter
  WCIDHeight            = 36.200*m; //"" "" height
  WCBarrelPMTOffset     = 0.0715*m; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 13.51;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}

void WCSimDetectorConstruction::Cylinder_60x74_20inchBandL_14perCent()
{ 
  WCDetectorName = "Cylinder_60x74_20inchBandL_14perCent()";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine20inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 74.0*m;
  WCIDHeight            = 60.0*m;
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 13.51;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}

void WCSimDetectorConstruction::Cylinder_60x74_20inchBandL_40perCent()
{ 
  WCDetectorName = "Cylinder_60x74_20inchBandL_40perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine20inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 74.0*m;
  WCIDHeight            = 60.0*m;
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 40.0;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}

void WCSimDetectorConstruction::Cylinder_60x74_20inchBandL_20perCent()
{ 
  WCDetectorName = "Cylinder_60x74_20inchBandL_20perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine20inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 74.0*m;
  WCIDHeight            = 60.0*m;
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 20.0;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}


void WCSimDetectorConstruction::Cylinder_12inchHPD_15perCent()
{ 
  WCDetectorName = "Cylinder_12inchHPD_15perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  // cylindrical detector with a height of 100m and a diameter of 69m 
  // with 12" HPD and 14.59% photocoverage
  WCSimPMTObject * PMT = CreatePMTObject("HPD12inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 69.0*m;
  WCIDHeight            = 100.0*m;
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 14.59;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}

void WCSimDetectorConstruction::SetHyperKGeometry()
{
  WCDetectorName = "HyperK";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("BoxandLine20inchHQE", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();
  WCIDDiameter          = 70.8*m; // = 74m - 2*(60cm ID wall + 1m OD)
  WCIDHeight            = 54.8*m; // = 60m - 2*(60cm ID wall + 2m OD)
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 40.0;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;
}

void WCSimDetectorConstruction::SetEggShapedHyperKGeometry()
{
  WCDetectorName = "EggShapedHyperK";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCODCollectionName = WCDetectorName + "-glassFaceWCPMT_OD"; 
  WCSimPMTObject * PMT = CreatePMTObject("PMT20inch", WCIDCollectionName[0]);
  WCSimPMTObject * outerPMT = CreatePMTObject("PMT8inch", WCODCollectionName);
  WCPMTName[0] = PMT->GetPMTName();
  innerPMT_Expose = PMT->GetExposeHeight();
  innerPMT_Radius = PMT->GetRadius();
  waterTank_TopR   = 32000.*mm;
  waterTank_BotR   = 30000.*mm;
  waterTank_Height = 48000.*mm;
  waterTank_UpperA =  8000.*mm;
  waterTank_LowerB =  6000.*mm;
  waterTank_Length = 49500.*mm;

  innerPMT_TopR     = 29095.*mm; 
  innerPMT_BotR     = 27095.*mm;
  innerPMT_TopW     = 12038.*mm;
  innerPMT_BotW     = 11004.*mm;
  innerPMT_Height   = 21095.*mm;
  innerPMT_Rpitch   =   990.*mm;
  innerPMT_Apitch   =   990.*mm;

  outerPMT_Name = outerPMT->GetPMTName();
  outerPMT_Expose = outerPMT->GetExposeHeight();
  outerPMT_Radius = outerPMT->GetRadius();
  outerPMT_TopR      = innerPMT_TopR + 900.*mm;
  outerPMT_BotR      = innerPMT_BotR + 900.*mm;
  outerPMT_TopW      = 12394.*mm;
  outerPMT_BotW      = 11319.*mm;
  outerPMT_Height    = innerPMT_Height + 900.*mm;
  outerPMT_TopRpitch = 3. * innerPMT_Rpitch * (outerPMT_TopR/innerPMT_TopR);
  outerPMT_BotRpitch = 3. * innerPMT_Rpitch * (outerPMT_BotR/innerPMT_BotR);
  outerPMT_Apitch    = 2. * innerPMT_Apitch;

  blackSheetThickness = 20.*mm;

  innerPMT_TopN = 0;
  innerPMT_BotN = 0;

  isEggShapedHyperK = true; // Tell DetectorConstruction to build egg-shaped HK geometry

  MatchWCSimAndEggShapedHyperK();
  InitSinglePMT();
}

void WCSimDetectorConstruction::SetEggShapedHyperKGeometry_withHPD()
{
  WCDetectorName = "EggShapedHyperK_withHPD";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCODCollectionName = WCDetectorName + "-glassFaceWCPMT_OD";
  WCSimPMTObject * PMT = CreatePMTObject("HPD20inchHQE", WCIDCollectionName[0]);
  WCSimPMTObject * outerPMT = CreatePMTObject("PMT8inch", WCODCollectionName);
  WCPMTName[0] = PMT->GetPMTName();
  innerPMT_Expose = PMT->GetExposeHeight();
  innerPMT_Radius = PMT->GetRadius();
  waterTank_TopR   = 32000.*mm;
  waterTank_BotR   = 30000.*mm;
  waterTank_Height = 48000.*mm;
  waterTank_UpperA =  8000.*mm;
  waterTank_LowerB =  6000.*mm;
  waterTank_Length = 49500.*mm;


  innerPMT_TopR     = 29095.*mm;
  innerPMT_BotR     = 27095.*mm;
  innerPMT_TopW     = 12038.*mm;
  innerPMT_BotW     = 11004.*mm;
  innerPMT_Height   = 21095.*mm;
  innerPMT_Rpitch   =   990.*mm;
  innerPMT_Apitch   =   990.*mm;


  outerPMT_Expose = outerPMT->GetExposeHeight();
  outerPMT_Radius = outerPMT->GetRadius();
  outerPMT_TopR      = innerPMT_TopR + 900.*mm;
  outerPMT_BotR      = innerPMT_BotR + 900.*mm;
  outerPMT_TopW      = 12394.*mm;
  outerPMT_BotW      = 11319.*mm;
  outerPMT_Height    = innerPMT_Height + 900.*mm;
  outerPMT_TopRpitch = 3. * innerPMT_Rpitch * (outerPMT_TopR/innerPMT_TopR);
  outerPMT_BotRpitch = 3. * innerPMT_Rpitch * (outerPMT_BotR/innerPMT_BotR);
  outerPMT_Apitch    = 2. * innerPMT_Apitch;

  blackSheetThickness = 20.*mm;

  innerPMT_TopN = 0;
  innerPMT_BotN = 0;

  isEggShapedHyperK = true; // Tell DetectorConstruction to build egg-shaped HK geometry

  MatchWCSimAndEggShapedHyperK();
  InitSinglePMT();
}


/**
 * Transfer egg-shaped HK variables needed elsewhere
 * to their generic WC equivalents.
 * This should be included in all egg-shaped HK configurations.
 */
void WCSimDetectorConstruction::MatchWCSimAndEggShapedHyperK()
{
  WCLength = waterTank_Length;
  WCPosition = 0.;
  WCPMTRadius[0] = innerPMT_Radius;
}




// Trevor Towstego: detector with single module
// ToDo: check defaults for mPMT.
void WCSimDetectorConstruction::SetTestSinglemPMTGeometry()
{
  WCDetectorName = "TestSinglemPMT";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";

  mPMT_ID_PMT[0] = "PMT3inchR12199_02"; //"BoxandLine20inchHQE";// (combine with nPMT=1);
  mPMT_OD_PMT[0] = "PMT8inch";          // Only used for the unique string name of mPMT for now!
  
  WCSimPMTObject * PMT = CreatePMTObject(mPMT_ID_PMT[0], WCIDCollectionName[0]);
  WCPMTName[0] = PMT->GetPMTName();
  WCPMTExposeHeight[0] = PMT->GetExposeHeight();
  WCPMTRadius[0] = PMT->GetRadius();
  
  //mPMT params
  vessel_cyl_height[0] = 38.*CLHEP::mm;
  vessel_radius_curv[0] = 342.*CLHEP::mm;
  vessel_radius[0] = 254.*CLHEP::mm;
  dist_pmt_vessel[0] = 2*CLHEP::mm;
  orientation[0] = PERPENDICULAR;
  mPMT_outer_material[0] = "G4_PLEXIGLASS";
  mPMT_inner_material[0] = "SilGel";
  mPMT_material_pmtAssembly[0] = "SilGel";
  mPMT_outer_material_d[0] = 10.*CLHEP::mm;
  // Radius of cone at z=reflectorHeight
  id_reflector_height[0] = 6.53*CLHEP::mm;                //10. > previous 7mm (deprecated number from JINST)
  id_reflector_z_offset[0] = 4.8*CLHEP::mm;            //from KM3Net CAD drawings
  id_reflector_angle[0] = 48*CLHEP::pi/180.*CLHEP::rad; // Based on KM3Net reflector specs
  mPMT_pmt_openingAngle[0] = 8.7*CLHEP::deg;     
  
  // BarrelPMTOffset/WCCapEdgeLimit needs PMT/mPMT height
  G4double vessel_tot_height = vessel_radius[0] + vessel_cyl_height[0];
  // parameters related to filling the ID mPMT
  nID_PMTs[0] = 19; //33;
  config_file[0] = "mPMTconfig_19_nuPrism_3ring.txt";//"mPMTconfig_33_13_1.txt";

  WCIDDiameter             = 10.0*m;
  WCIDHeight               = 30.0*m;
  WCBarrelPMTOffset        = vessel_tot_height;
  WCBarrelNumPMTHorizontal = 3;
  WCBarrelNRings           = 3;
  WCPMTperCellHorizontal   = 1;
  WCPMTperCellVertical     = 1;
  WCPMTPercentCoverage     = 10.;
  WCCapEdgeLimit           = WCIDDiameter/2.0 - vessel_tot_height;
  WCCapPMTSpacing          = WCIDDiameter*2.0;
  WCBlackSheetThickness    = 2.0*cm;
  WCAddGd                  = false;
}




void WCSimDetectorConstruction::Cylinder_60x74_hybrid_40perCent()
{ 
  WCDetectorName = "Cylinder_60x74_hybrid_40perCent";

  IDnoTypes = 2;
  // ID PMT 1
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  mPMT_ID_PMT[0] = "BoxandLine20inchHQE"; 
  WCSimPMTObject * PMT = CreatePMTObject(mPMT_ID_PMT[0], WCIDCollectionName[0]);
  WCPMTName[0]         = PMT->GetPMTName();
  WCPMTRadius[0]       = PMT->GetRadius();

  // ToDo: construct 20" B&L protective cover, now all set to zero
  // potentially even with SilGel
  vessel_cyl_height[0] = 0.*CLHEP::mm;
  vessel_radius_curv[0] = 0.*CLHEP::mm;
  vessel_radius[0] = 0.*CLHEP::mm;
  dist_pmt_vessel[0] = 0*CLHEP::mm;
  orientation[0] = PERPENDICULAR;
  mPMT_outer_material[0] = "Water";
  mPMT_inner_material[0] = "Water";
  mPMT_material_pmtAssembly[0] = "Water";
  mPMT_outer_material_d[0] = 0.*CLHEP::mm;
  // Radius of cone at z=reflectorHeight
  id_reflector_height[0] = 0.*CLHEP::mm;                //10. > previous 7mm (deprecated number from JINST)
  id_reflector_z_offset[0] = 0.*CLHEP::mm;            //from KM3Net CAD drawings
  id_reflector_angle[0] = 0*CLHEP::pi/180.*CLHEP::rad; // Based on KM3Net reflector specs
  mPMT_pmt_openingAngle[0] = 0*CLHEP::deg;     
  
  // BarrelPMTOffset/WCCapEdgeLimit needs PMT/mPMT height
  G4double vessel_tot_height = vessel_radius[0] + vessel_cyl_height[0];
  // parameters related to filling the ID mPMT
  nID_PMTs[0] = 1; 
  config_file[0] = "";


  // ID PMT 2
  WCIDCollectionName[1] = WCDetectorName +"-glassFaceWCmPMT"; 
  mPMT_ID_PMT[1] = "PMT3inchR12199_02"; 
  mPMT_OD_PMT[1] = "PMT3inchR12199_02";          
  
  WCSimPMTObject * PMT2 = CreatePMTObject(mPMT_ID_PMT[1], WCIDCollectionName[1]);
  WCPMTName[1] = PMT2->GetPMTName();
  WCPMTRadius[1] = PMT2->GetRadius();
  
  //mPMT params
  vessel_cyl_height[1] = 38.*CLHEP::mm;
  vessel_radius_curv[1] = 342.*CLHEP::mm;
  vessel_radius[1] = 254.*CLHEP::mm;
  dist_pmt_vessel[1] = 2*CLHEP::mm;
  orientation[1] = PERPENDICULAR;
  mPMT_outer_material[1] = "G4_PLEXIGLASS";
  mPMT_inner_material[1] = "SilGel";
  mPMT_material_pmtAssembly[1] = "SilGel";
  mPMT_outer_material_d[1] = 10.*CLHEP::mm;
  // Radius of cone at z=reflectorHeight
  id_reflector_height[1] = 6.53*CLHEP::mm;                //10. > previous 7mm (deprecated number from JINST)
  id_reflector_z_offset[1] = 4.8*CLHEP::mm;            //from KM3Net CAD drawings
  id_reflector_angle[1] = 48*CLHEP::pi/180.*CLHEP::rad; // Based on KM3Net reflector specs
  mPMT_pmt_openingAngle[1] = 8.7*CLHEP::deg;     
  
  // BarrelPMTOffset/WCCapEdgeLimit needs PMT/mPMT height
  G4double vessel_tot_height2 = vessel_radius[1] + vessel_cyl_height[1];
  // parameters related to filling the ID mPMT
  nID_PMTs[1] = 19; 
  config_file[1] = "mPMTconfig_19_nuPrism_3ring.txt";

  //Pick largest
  G4double largest_vessel_radius = std::max({PMT->GetRadius(), vessel_radius[0],vessel_radius[1]}); //C++11
  G4double largest_offset = std::max({WCPMTRadius[0],vessel_tot_height,vessel_tot_height2});

  WCPMTRadius[0]           = largest_vessel_radius; 
  WCIDDiameter          = 74.0*m;
  WCIDHeight            = 60.0*m;
  WCBarrelPMTOffset     = largest_offset; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 4;    // instead of 3 to make sure barrel is filled alternatingly
  WCPMTPercentCoverage  = 40.0;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*largest_vessel_radius));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - largest_vessel_radius;
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

}






// Note: using vessel_radius (mPMT radius) instead of PMT radius for detector construction. Important!
void WCSimDetectorConstruction::Cylinder_60x74_3inchmPMT_40perCent()
{ 
  WCDetectorName = "Cylinder_60x74_3inchmPMT_40perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  mPMT_ID_PMT[0] = "PMT3inchR12199_02";
  mPMT_OD_PMT[0] = "PMT8inch";
  WCSimPMTObject * PMT = CreatePMTObject(mPMT_ID_PMT[0], WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = PMT->GetRadius();

  //mPMT params:
  vessel_cyl_height[0] = 38.*CLHEP::mm; 
  vessel_radius_curv[0] = 342.*CLHEP::mm;   
  vessel_radius[0] = 254.*CLHEP::mm;   
  dist_pmt_vessel[0] = 2*CLHEP::mm;
  orientation[0] = PERPENDICULAR;
  mPMT_outer_material[0] = "G4_PLEXIGLASS";
  mPMT_inner_material[0] = "SilGel";
  mPMT_material_pmtAssembly[1] = "SilGel";
  mPMT_outer_material_d[0] = 10.*CLHEP::mm;
  // Radius of cone at z=reflectorHeight
  id_reflector_height[0] = 6.53*CLHEP::mm;                //10. > previous 7mm (deprecated number from JINST)
  id_reflector_z_offset[0] = 4.8*CLHEP::mm;            //from KM3Net CAD drawings
  id_reflector_angle[0] = 48*CLHEP::pi/180.*CLHEP::rad; // Based on KM3Net reflector specs
  mPMT_pmt_openingAngle[0] = 8.7*CLHEP::deg;     
  G4double vessel_tot_height = vessel_radius[0] + vessel_cyl_height[0];

  // parameters related to filling the ID mPMT
  nID_PMTs[0] = 19;
  config_file[0] = "mPMTconfig_19_nuPrism_3ring.txt";


  WCIDDiameter          = 74.0*m;
  WCIDHeight            = 60.0*m;
  WCBarrelPMTOffset     = vessel_tot_height; //mPMT cylinder Radius //WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 40.0;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*vessel_radius[0]));
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - vessel_radius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;


}

// Uniform distribution of 3" PMTs
void WCSimDetectorConstruction::Cylinder_60x74_3inch_40perCent()
{ 
  WCDetectorName = "Cylinder_60x74_3inch_40perCent";
  WCIDCollectionName[0] = WCDetectorName +"-glassFaceWCPMT";
  WCSimPMTObject * PMT = CreatePMTObject("PMT3inchR12199_02", WCIDCollectionName[0]);
  WCPMTName[0]           = PMT->GetPMTName();
  WCPMTExposeHeight[0]   = PMT->GetExposeHeight();
  WCPMTRadius[0]         = .04*m; //PMT->GetRadius();
  WCIDDiameter          = 74.0*m;
  WCIDHeight            = 60.0*m;
  WCBarrelPMTOffset     = WCPMTRadius[0]; //offset from vertical
  WCPMTperCellHorizontal= 4;
  WCPMTperCellVertical  = 3;
  WCPMTPercentCoverage  = 40.;
  WCBarrelNumPMTHorizontal = round(WCIDDiameter*sqrt(pi*WCPMTPercentCoverage)/(10.0*WCPMTRadius[0]));  
  WCBarrelNRings           = round(((WCBarrelNumPMTHorizontal*((WCIDHeight-2*WCBarrelPMTOffset)/(pi*WCIDDiameter)))
                                      /WCPMTperCellVertical));
  WCCapPMTSpacing       = (pi*WCIDDiameter/WCBarrelNumPMTHorizontal); // distance between centers of top and bottom pmts
  WCCapEdgeLimit        = WCIDDiameter/2.0 - WCPMTRadius[0];
  WCBlackSheetThickness = 2.0*cm;
  WCAddGd               = false;

  InitSinglePMT();
}



void WCSimDetectorConstruction::InitSinglePMT(){
  
  for(int i = 0; i<2; i++){
	vessel_cyl_height[i] = 0.*mm;
	vessel_radius_curv[i] = 0.1*mm;
	vessel_radius[i] = 0.1*mm;
	dist_pmt_vessel[i] = 0.*mm;
	orientation[i] = PERPENDICULAR;
	mPMT_ID_PMT[i] = "";
	mPMT_OD_PMT[i] = "";
	mPMT_outer_material[i] = "Water";
	mPMT_inner_material[i] = "Water";
	mPMT_material_pmtAssembly[i] = "Water";
	mPMT_outer_material_d[i] = 0.*CLHEP::mm;
	id_reflector_height[i] = 0.*CLHEP::mm;
	id_reflector_z_offset[i] = 0.*CLHEP::mm;
	id_reflector_angle[i] = 0.*CLHEP::rad; 
	nID_PMTs[i] = 1;   
	config_file[i] = "";
	mPMT_pmt_openingAngle[i] = 0*CLHEP::deg;     
  }
}
