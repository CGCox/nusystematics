#pragma once

#include "nusystematics/interface/IGENIESystProvider_tool.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH3D.h"
#include "TLorentzVector.h"


#include <array>
#include <memory>
#include <string>

class CVspread : public nusyst::IGENIESystProvider_tool {

public:
  explicit CVspread(fhicl::ParameterSet const &);

  //'First' configuration step: tool configuration
  // - takes 'arbitrary' FHiCL configuration and
  //   produces SystMetaData object which can later be used to configure a
  //   specific set of parameter values to be calculated
  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &,
                                            systtools::paramId_t);

  //'Second' configuration step: parameter headers
  // - Reads the preconstructed SystMetaData produced by BuildSystMetaData
  //   to configure an instance of this class to calculate weights
  // - Recieves a copy of the tool_options instance constructed by 
  //   BuildSystMetaData as an argument
  bool SetupResponseCalculator(fhicl::ParameterSet const &);

  // Used to pass arbitrary FHiCL options from the tool configuration to the
  //   parameter headers.
  fhicl::ParameterSet GetExtraToolOptions() { return tool_options; }

  systtools::SystMetaData BuildSystMetaData(fhicl::ParameterSet const &, systtools::paramId_t);

  // Parameter-specific implementation goes in here
  systtools::event_unit_response_t GetEventResponse(genie::EventRecord const &);

  // Can add as much or as little stateful information here for use when
  // representing this instance as a string.
  std::string AsString() { return "CVspread"; }

  ~CVspread(){}

private:
  // arbitrary additional configuration from the tool configuration/parameter
  // headers can be storeds here
  fhicl::ParameterSet tool_options;

  // The ParamHeaders id of the free parameters provided by this systprovider
  std::vector<size_t> pidx_Params;

  // The configured variations to precalculate for the parameters
  std::vector<double> CVs;
  std::vector<std::vector<double>> Variations;
 
  void LoadObject();

  TFile *RefROOT; //TFile for the reference root file containing our 3D object
  TH3D *ObjRef; //Name of our 3D object (TH3D histogram)

  TLorentzVector ISLepP4, FSLepP4, emTransfer;

  struct DialInfo { //preface to a multi-generator option beyond NEUT
    std::string prettyname;
  };
  std::vector<DialInfo> dial_infos;

  // configurable verbosity as an example of some arbitrary systprovider
  // configuration
  int verbosity_level;

  float Enu_true, Q2, q0, m_nucleon, W_nuc_rest;

  const double m_proton = 0.93827203;   // Proton mass in GeV - set identically to nuisance

};