#include "nusystematics/systproviders/CVspread_tool.hh"

#include "systematicstools/utility/FHiCLSystParamHeaderUtility.hh"

#include "RwFramework/GSyst.h"

 //GENIE-MC/Generator
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"

#include <cmath>

using namespace nusyst;
using namespace systtools;

// constructor passes up configuration object to base class for generic tool
// initialization and initialises our local copies of paramIds to unconfigured
//  flag values
CVspread::CVspread(fhicl::ParameterSet const &params)
    : IGENIESystProvider_tool(params) { //get values from associated fcl file

      //below are defined in .hh file
      pidx_Params.push_back(kParamUnhandled<size_t>);
      dial_infos = {"GENIE2NEUT"};
      CVs.push_back(0);
      Variations.emplace_back();
      ReWeightEngines.emplace_back();
} //close param


//Setting up the loading of histogram
void CVspread::LoadObject() {
    // Open the reference ROOT file
   RefRoot  = TFile::Open("reference3D.root", "READ");
    if (!RefRoot || RefRoot->IsZombie()) {
        std::cerr << "Error: Could not open -.root file which stores 'mapping object' reference." << std::endl;
        return;
    }

    // Retrieve the histogram
    ObjRef = RefRoot->Get<TH3D>("MapReference");
    if (!ObjRef) {
        std::cerr << "Error: Could not find 'mapping object' AKA histogram stored as a TObject." << std::endl;
    }
}


//Reads the fcl file - AKA, user deciding which parameter dials, and tweak values, and this reads it
SystMetaData CVspread::BuildSystMetaData(fhicl::ParameterSet const &ps,
                                            paramId_t firstId) {
  SystMetaData smd;

  // loop through the header struct
  for (auto const &prettyname : dial_infos) {
    SystParamHeader phdr; //short-hand for simplicity
    phdr.prettyName = prettyname; //keep things neat - fcl name for tweak

    // Set up parameter definition with a standard tool configuration form using helper function.
    phdr.systParamId = firstId++;

    // set any parameter-specific ParamHeader metadata here
    phdr.isSplineable = true;
    phdr.Variations = {0,1};

    // add it to the metadata list to pass back.
    smd.push_back(phdr);
  }

  // Put any options that you want to propagate to the ParamHeaders options
  tool_options.put("verbosity_level", ps.get<int>("verbosity_level", 0));

  return smd;
}

//Sets up the response from the fcl values; if set up correctly, reference histogram is loaded
bool CVspread::SetupResponseCalculator(
    fhicl::ParameterSet const &tool_options) {
  verbosity_level = tool_options.get<int>("verbosity_level", 0);

  // grab the pre-parsed param headers object
  SystMetaData const &md = GetSystMetaData();

  for (size_t i = 0; i < dial_infos.size(); ++i) {

    auto const &prettyname = dial_infos[i].prettyname;

    if (!HasParam(md, prettyname)) {
      if (verbosity_level > 1) {
        std::cout << "[INFO]: Don't have parameter " << prettyname
                  << " in SystMetaData. Skipping configuration." << std::endl;
      }
      continue;
    }
    pidx_Params[i] = GetParamIndex(md, prettyname);

    if (verbosity_level > 1) {
      std::cout << "[INFO]: Have parameter " << prettyname
                << " in SystMetaData with ParamId: " << pidx_Params[i]
                << ". Configuring." << std::endl;
    }

    auto phdr = md[pidx_Params[i]]; // param_md

    CVs[i] = phdr.centralParamValue; // the initial central value of parameter
    Variations[i] = phdr.paramVariations; // the variation of parameter


    //Temporary throw loop until fcl integer step resolved - clarification needed
    //currently set to return only two weights - one for tweak = 0 [GENIE], and = 1 [NEUT]
    if (Variations[i].size() > 2){
      std::cout << "Only one tweak value may be given at this time." << std::endl; 
      throw; 
    }

    // in the concrete version of this example we're going to configure a GENIE GReWeightNuXSecCCQE instance for each configured variation.
    for (auto v : Variations[i]) {  // loop through the variations of each param e.g (-3,-2...,2,3)std in fcl file
      
        if (v < 0 || v > 1) {
    //this will check to make sure param variations are within range
          std::cout << "Tweak out of range - dial range is between 0 and 1." << std::endl;
    throw;
  }
    }
  }
  LoadObject();
  // returning cleanly
  return true;
}

//Produces a weight value based on the fcl values, and returns the weight (as "resp")
event_unit_response_t CVspread::GetEventResponse(genie::EventRecord const &ev){

  event_unit_response_t resp;

  SystMetaData const &md = GetSystMetaData();


  //Below explained in GHepRecord.cxx
  //Returns the intial state probe particle i.e. the muon neutrino
  genie::GHepParticle *ISLep = ev.Probe(); //Probe() constructed in GHepRecord.h for GHepParticle.h
  
  //Returns the final state primary lepton
  genie::GHepParticle *FSLep = ev.FinalStatePrimaryLepton(); //FSPL() constructed in GHepRecord.h for GHepParticle.h

  genie::GHepParticle *TargetNucleon = ev.TargetNucleus(); //TargetNucleus() constructed in GHepRecord.h for GHepParticle.h

  /*  if (!FSLep || !ISLep) {
    throw incorrectly_generated()
        << "[ERROR]: Failed to find IS and FS lepton in event: "
        << ev.Summary()->AsString();
  }*/

  //P4() is a TLorentzVector constructed in GHepParticle.h; stores Enu_true, px, py and pz components for particle
  ISLep4 = *ISLep->P4(); //Points to initial state four-vector
  FSLepP4 = *FSLep->P4(); //Points to final state four-vector
  emTransfer = (ISLepP4 - FSLepP4); //Calculates the energy-momentum transfer

  //Initial state energy AKA true neutrino energy
  Enu_true = ISLepP4.E();
  //Enu_true = nu->fP.E(); / 1E3; //nuisance definition
  //Enu_true = ISLep->P4.E();

  //Calculates the Q2 (the negatove of energy-momentum transfer squared)
  Q2 = -emTransfer.Mag2(); //Mag2() just takes the square of the magntidue

  q0 = emTransfer.E();

  //Nuisance defines this as just the mass of the proton
  //rather than extracting the ID of the target nucleon, returning the mass accordingly

  //Will extract the nucleon mass from event so the process exists; but will use proton mass until Nuisnace is changed
  //m_nucleon = *TargetNucleon->Mass();

  //W_nuc_rest = sqrt((-Q2)+(2*m_nucleon*q0)+(m_nucleon*m_nucleon)); //W_nuc_rest calculated for target nucleon

  W_nuc_rest = sqrt((-Q2)+(2*m_proton*q0)+(m_proton*m_proton)); //W_nuc_rest matching nuisance



  //Determine bin from calculated W_nuc_rest (x), Q2 (y), and Enu_true (z)
  int binX = ObjRef->GetXaxis()->FindBin(W_nuc_rest);
  int binY = ObjRef->GetYaxis()->FindBin(Q2);
  int binZ = ObjRef->GetZaxis()->FindBin(Enu_true);

  int binREF = ObjRef->GetBin(binX, binY, binZ);
            
  // Get the bin content
  double binContent = ObjRef->GetBinContent(binREF);


//Determine bin weight based on tweaks

  // loop through and calculate weights
  for (size_t i = 0; i < dial_infos.size(); ++i) {
    // this parameter wasn't configured, nothing to do
    if (pidx_Params[i] == kParamUnhandled<size_t>) {
      continue;
    }

    auto param_id = md[pidx_Params[i]].systParamId;
    // initialize the response array with this paramId
    resp.push_back({param_id, {}});

    // loop through variations for this parameter
    for (auto tweak : Variations[i]) {

      // put the response weight for this variation of this parameter into the
      // response object
      double REFweight = 1 + tweak*(binContent-1);
      resp.back().responses.push_back(REFweight);
      if (verbosity_level > 3) {
        std::cout << "[DEBG]: For parameter " << md[pidx_Params[i]].prettyName
                  << " at variation[" << tweak << "] = " << Variations[i][tweak]
                  << " calculated weight: " << resp.back().responses.back()
                  << std::endl;
      }
    }
  }

  return resp; //returns the final calculated, tweak adjusted, weight(s)
}

//Nuisance Reference Info
//MCStudies/GenericFlux_Vectors.cxx useful reference

  //Returns the invariant hadronic mass
  //W = ev.Summary()->Kine().W(true);
  //Summary() constructed in GENIE~Framework/GHEP/GHepRecord.h
  //Kine() constructed in GENIE~Framework/Interaction/Interaction.h
  //W(bool) constructed in GENIE~Framework/Interaction/Kinematics.h

  //Invariant Mass assuming nucleon is at rest
  //Defined in nuisance, not nusystematics - must calculate internally
      // Get W_true with assumption of initial state nucleon at rest

//nuisance calculation
/*    
    
    float m_n = (float)PhysConst::mass_proton;
    // Q2 assuming nucleon at rest
    W_nuc_rest = sqrt(-Q2 + 2 * m_n * q0 + m_n * m_n);
    W = W_nuc_rest; // For want of a better thing to do
    // True Q2
    x = Q2 / (2 * m_n * q0);
    y = 1 - ELep / Enu_true;*/

  /*nuisance definitions
  m_n = extracted from PhysConst namespace defined in nuisance/src/Utils/PhysConst.h - only proton mass
  q0 = (nu->fP - lep->fP).E() / 1E3;
  Q2 = -1 * (nu->fP - lep->fP).Mag2() / 1E6;

  1Ex is a preset factor for nuisance for scaling the event histogram - not relevant for nusystematics
  */
