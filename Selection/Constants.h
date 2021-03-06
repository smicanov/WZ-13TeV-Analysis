// ### Preselection ###

// ## electrons ##

// Barrel and Endcap
#define ETASCBARREL  1.479 
#define ETASCENDCAP  2.5
#define ELE_ETA_BARREL  1.479
#define ELE_ETA_ENDCAP  2.5

// electrons ID bits and Pt & Eta cuts
#define ELELOOSE_BIT  1
#define ELEMEDIUM_BIT  2
#define ELETIGHT_BIT  3
#define ELE_PTCUT  10.
#define ELE_ETACUT  2.5


// Cut based ID
// Loose
//Barrel
#define FULL5x5_SIGMAIETAIETA_BARREL_LOOSE  0.0103
#define DETAIN_BARREL_LOOSE  0.0105
#define DPHIIN_BARREL_LOOSE  0.115
#define HOVERE_BARREL_LOOSE  0.104
#define RELISO_BARREL_LOOSE  0.0893
#define OOEMOOP_BARREL_LOOSE  0.102
#define D0_BARREL_LOOSE  0.0261
#define DZ_BARREL_LOOSE  0.41
#define EXPMISSINNERHITS_BARREL  2
#define PASSCONVERSIONVETO  1

//Endcap
#define FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE  0.0301
#define DETAIN_ENDCAP_LOOSE  0.00814
#define DPHIIN_ENDCAP_LOOSE  0.182
#define HOVERE_ENDCAP_LOOSE  0.0897
#define RELISO_ENDCAP_LOOSE  0.121
#define OOEMOOP_ENDCAP_LOOSE  0.126
#define D0_ENDCAP_LOOSE  0.118
#define DZ_ENDCAP_LOOSE  0.822
#define EXPMISSINNERHITS_ENDCAP 1

// Medium
//Barrel
#define FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM  0.0101
#define DETAIN_BARREL_MEDIUM  0.0103
#define DPHIIN_BARREL_MEDIUM  0.0336
#define HOVERE_BARREL_MEDIUM  0.0876
#define RELISO_BARREL_MEDIUM  0.0766
#define OOEMOOP_BARREL_MEDIUM  0.0174
#define D0_BARREL_MEDIUM  0.0118
#define DZ_BARREL_MEDIUM  0.373

//Endcap
#define FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM  0.0283
#define DETAIN_ENDCAP_MEDIUM  0.00733
#define DPHIIN_ENDCAP_MEDIUM  0.114
#define HOVERE_ENDCAP_MEDIUM  0.0678
#define RELISO_ENDCAP_MEDIUM  0.0678
#define OOEMOOP_ENDCAP_MEDIUM  0.0898
#define D0_ENDCAP_MEDIUM  0.0739
#define DZ_ENDCAP_MEDIUM  0.602

// Tight
//Barrel
#define FULL5x5_SIGMAIETAIETA_BARREL_TIGHT  0.0101
#define DETAIN_BARREL_TIGHT  0.00926
#define DPHIIN_BARREL_TIGHT  0.0336
#define HOVERE_BARREL_TIGHT  0.0597
#define RELISO_BARREL_TIGHT  0.0354
#define OOEMOOP_BARREL_TIGHT  0.012
#define D0_BARREL_TIGHT  0.0111
#define DZ_BARREL_TIGHT  0.0466

//Endcap
#define FULL5x5_SIGMAIETAIETA_ENDCAP_TIGHT  0.0279
#define DETAIN_ENDCAP_TIGHT  0.00724
#define DPHIIN_ENDCAP_TIGHT  0.0918
#define HOVERE_ENDCAP_TIGHT  0.0615
#define RELISO_ENDCAP_TIGHT  0.0646
#define OOEMOOP_ENDCAP_TIGHT  0.00999
#define D0_ENDCAP_TIGHT  0.0351
#define DZ_ENDCAP_TIGHT  0.417

// ## muons ##

#define MU_ETA_BARREL  1.2
#define MU_ETA_ENDCAP  2.5
#define GLOBALMUON_BIT  1
#define TRACKERMUON_BIT  2
#define PFMUON_BIT  5
#define MU_RELISO_LOOSE  0.2
#define MU_RELISO_TIGHT  0.12
#define MU_PTCUT  10.
#define MU_ETACUT  2.4

// ## general ##

#define N_TIGHTLEPTONS  3


// ### Z selection ###

#define MZ_MIN  60.
#define MZ_MAX  120.
#define PDG_ZMASS  91.1876
#define ZLEADINGLEPTON_PTCUT  20.


// ### W selection ###

#define WLEPTON_PTCUT  20.
#define WZ_DELTARCUT  0.1

// ### all combinations of 2 selected leptons inv. mass cut ###

#define MASS2LCUT  4.

// ### MET ###

#define METCUT  30.


// ### total inv. mass of 3 candidate leptons ###

#define MASS3LCUT  100.


// ### Jets ###

#define JET_PTCUT  30.
#define JET_ETACUT  2.5
#define ELEJET_DELTARMIN  0.1
#define MUJET_DELTARMIN  0.1
#define LEPTONJET_DELTARMIN  0.1
