// -*- C++ -*-
//
// HJetsProcessInfo.h is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HJets_HJetsProcessInfo_H
#define HJets_HJetsProcessInfo_H

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace HJets {

using namespace ThePEG;

/**
 * Amplitude configuration to include
 */
struct AmplitudeInfo {

  /**
   * ij line quarks
   */
  pair<int,int> ijLine;

  /**
   * kl line quarks
   */
  pair<int,int> klLine;

  /**
   * Emissions on the ij line
   */
  pair<int,int> ijLineEmissions;

  /**
   * Emissions on the kl line
   */
  pair<int,int> klLineEmissions;

  /**
   * The sign factor for crossed fermions.
   */
  double fermionSign;

  /**
   * Wether we have a qqbar pair emitted
   */
  bool qqbarEmitted;

  /**
   * The exchanged boson
   */
  cPDPtr boson;

  /**
   * The boson mass
   */
  Energy bosonMass;

  /**
   * The boson width
   */
  Energy bosonWidth;

  /**
   * The left-handed coupling factor
   */
  complex<double> ijLeftCoupling;

  /**
   * The right-handed coupling factor
   */
  complex<double> ijRightCoupling;

  /**
   * The left-handed coupling factor
   */
  complex<double> klLeftCoupling;

  /**
   * The right-handed coupling factor
   */
  complex<double> klRightCoupling;

  /**
   * The Higgs coupling
   */
  complex<Energy> higgsCoupling;

  /**
   * The colour tensors to consider with respective weights
   */
  map<size_t,double> colourTensors;

  /**
   * Dump to ostream for diagnostics
   */
  void print(ostream&) const;

};

/**
 * Information on the generic H+2jets process to consider
 */
struct HJetsProcessInfo {

  /**
   * Get prototype configurations for the given physical process
   */
  static vector<AmplitudeInfo> configurations(const cPDVector& proc,
					      const StandardModelBase& sm,
					      bool complexMassScheme = true,
					      double kappaZ = 1.0,
					      double kappaW = 1.0);

  /**
   * Check if the given process corresponds to a valid H+>=2jets
   * configuration.
   */
  static bool canHandle(const PDVector& proc,
			const StandardModelBase& sm);

  /**
   * Generate all configurations to consider for evaluating the
   * amplitude.
   */
  static vector<AmplitudeInfo> getConfigurations(const cPDVector& proc,
						 const vector<int>& crossed2Proc,
						 const StandardModelBase& sm,
						 bool complexMassScheme = true,
						 double kappaZ = 1.0,
						 double kappaW = 1.0);

};

}

#endif // HJets_HJetsProcessInfo_H
