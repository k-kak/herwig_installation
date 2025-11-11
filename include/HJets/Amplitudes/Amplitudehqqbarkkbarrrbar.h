// -*- C++ -*-
//
// Amplitudehqqbarkkbarrrbar.h is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HJets_Amplitudehqqbarkkbarrrbar_H
#define HJets_Amplitudehqqbarkkbarrrbar_H
//
// This is the declaration of the Amplitudehqqbarkkbarrrbar class.
//

#include "AmplitudeBase.h"

namespace HJets {

using namespace ThePEG;
using namespace Herwig;

/**
 * \brief Tree level helicity amplitudes for 0 -> h q qbar k kbar g g
 */
class Amplitudehqqbarkkbarrrbar: 
    public AmplitudeBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Amplitudehqqbarkkbarrrbar();

  /**
   * The destructor.
   */
  virtual ~Amplitudehqqbarkkbarrrbar();
  //@}

public:

  /**
   * Return the number of partons attached to this amplitude.
   */
  virtual int nPartons() const { return 6; }

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector& proc) const;

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t, const vector<int>&, Complex&);

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Amplitudehqqbarkkbarrrbar & operator=(const Amplitudehqqbarkkbarrrbar &);

};

}

#endif /* HJets_Amplitudehqqbarkkbarrrbar_H */
