// -*- C++ -*-
//
// AmplitudeBase.h is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HJets_AmplitudeBase_H
#define HJets_AmplitudeBase_H
//
// This is the declaration of the AmplitudeBase class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxCurrents.h"
#include "HJetsProcessInfo.h"

namespace HJets {

using namespace ThePEG;
using namespace Herwig;

  //#define AMPVERBOSE

/**
 * \brief Tree level helicity amplitudes for 0 -> h + jets
 */
class AmplitudeBase: 
    public MatchboxAmplitude, public MatchboxCurrents {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  AmplitudeBase();

  /**
   * The destructor.
   */
  virtual ~AmplitudeBase();
  //@}

public:

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const;

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return 3; }

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const { return nPartons()-4; }

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {
    MatchboxCurrents::reset();
    MatchboxAmplitude::flushCaches();
  }

  /**
   * Return the number of partons attached to this amplitude.
   */
  virtual int nPartons() const = 0;

  /**
   * Return true, if one-loop contributions will be evaluated at amplitude level.
   */
  virtual bool oneLoopAmplitudes() const { return false; }

  /**
   * Return true, if one loop corrections are given in the convenctions
   * of everything expanded.
   */
  virtual bool isExpanded() const { return true; }

  //virtual bool isDR() const{return false;}

  /**
   * Return the virtual amplitude info
   */
  const map<int, pair<int, double> >& virtualInfo() const;

  /**
   * The virtuals info maps.
   */
  static map<cPDVector,map<int, pair<int, double> > >& virtualInfos();

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return mu2Factor*lastSHat(); }

protected:

  /**
   * Return the amplitude configurations to be considered.
   */
  vector<AmplitudeInfo>& amplitudeInfo();

  /**
   * Get the boson propagators and Higgs coupling factor
   */
  Complex bosonFactor(const AmplitudeInfo&) const;

  /**
   * Check if t-channel diagram is a neutral current diagram.
   */
  bool topologyOneIsNC() const;

  /**
   * Check if u-channel diagram is a neutral current diagram.
   */
  bool topologyTwoIsNC() const;

  /**
   * Check if t-channel diagram is a neutral current diagram.
   */
  bool topologyOneIsCC() const;

  /**
   * Check if u-channel diagram is a neutral current diagram.
   */
  bool topologyTwoIsCC() const;

  /**
   * Work out the couplings
   */
  void getCouplings(double*, double*, double*, double*, 
		    double*, double*, double*, double*, double*, double*, double*, double*, int*) const;


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /**
   * Return the factor to rescale the HZZ coupling
   */
  double kappaZ() const { return theKappaZ; }

  /**
   * Return the factor to rescale the HWW coupling
   */
  double kappaW() const { return theKappaW; }

protected:

  /**
   * Return the threshold for classifying points as unstable
   */
  double stableEpsilon() const { return theStableEpsilon; }

private:

  /**
   * The amplitude info map.
   */
  static map<cPDVector,vector<AmplitudeInfo> >& amplitudeInfos();

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AmplitudeBase & operator=(const AmplitudeBase &);

  /**
   * Wether or not the complex mass scheme should be used
   */
  bool complexMassScheme;

  /**
   * A rescaling factor for the t'Hooft mass
   */
  double mu2Factor;

  /**
   * A factor to rescale the HZZ coupling
   */
  double theKappaZ;

  /**
   * A factor to rescale the HWW coupling
   */
  double theKappaW;

  /**
   * The threshold for classifying points as unstable
   */
  double theStableEpsilon;

};

}

#endif /* HJets_AmplitudeBase_H */
