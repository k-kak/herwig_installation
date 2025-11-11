// -*- C++ -*-
//
// Amplitudehqqbarkkbargg.cc is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Amplitudehqqbarkkbargg class.
//

#include "Amplitudehqqbarkkbargg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace HJets;

Amplitudehqqbarkkbargg::Amplitudehqqbarkkbargg() {}

Amplitudehqqbarkkbargg::~Amplitudehqqbarkkbargg() {}

IBPtr Amplitudehqqbarkkbargg::clone() const {
  return new_ptr(*this);
}

IBPtr Amplitudehqqbarkkbargg::fullclone() const {
  return new_ptr(*this);
}

bool Amplitudehqqbarkkbargg::canHandle(const PDVector& proc) const {
  bool gotAGluon = false;
  for ( PDVector::const_iterator p = proc.begin(); p != proc.end(); ++p )
    if ( (**p).id() == ParticleID::g ) {
      gotAGluon = true;
      break;
    }
  if ( !gotAGluon )
    return false;
  return AmplitudeBase::canHandle(proc);
}

void Amplitudehqqbarkkbargg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  // inform the amplitude cache how many legs we are handling
  nPoints(7);

  // everything is expressed in dimensionless quantities, referring to
  // sqrt(shat) of the collission as the mass unit
  amplitudeScale(sqrt(lastSHat()));

  // the higgs momentum
  Lorentz5Momentum ph = amplitudeMomentum(0);
  momentum(0,ph,true,ph.mass());

  // the parton momenta
  for ( size_t k = 1; k < 7; ++k )
    momentum(k,amplitudeMomentum(k),true);

  MatchboxAmplitude::prepareAmplitudes(me);

}


Complex Amplitudehqqbarkkbargg::evaluate(size_t colourTensor, 
					 const vector<int>& helicities,
					 Complex& largeN) {

  // get configurations we need to consider
  const vector<AmplitudeInfo>& confs = amplitudeInfo();

  // sum of all contributions
  Complex result = 0;
  largeN = 0.;

  for ( vector<AmplitudeInfo>::const_iterator c = confs.begin();
	c != confs.end(); ++c ) {

    map<size_t,double>::const_iterator cx = c->colourTensors.find(colourTensor);
    if ( cx == c->colourTensors.end() )
      continue;

    Complex diagram = 0.;

    Complex ijL = c->ijLeftCoupling;
    Complex ijR = c->ijRightCoupling;
    Complex klL = c->klLeftCoupling;
    Complex klR = c->klRightCoupling;

    int i = c->ijLine.first;
    int j = c->ijLine.second;
    int k = c->klLine.first;
    int l = c->klLine.second;

    if ( c->ijLineEmissions.first > 0 && c->ijLineEmissions.second > 0 &&
	 c->klLineEmissions.first < 0 && c->klLineEmissions.second < 0 ) {
      int g1 = c->ijLineEmissions.first;
      int g2 = c->ijLineEmissions.second;
      if ( ijL != 0. && klL != 0. ) {
	diagram +=
	  ijL*klL*qqbarggLeftCurrent(i,helicities[i],
				     j,helicities[j],
				     g1,helicities[g1],
				     g2,helicities[g2]).
	  dot(qqbarLeftCurrent(k,helicities[k],
			       l,helicities[l]));
      }
      if ( ijL != 0. && klR != 0. ) {
	diagram +=
	  ijL*klR*qqbarggLeftCurrent(i,helicities[i],
				     j,helicities[j],
				     g1,helicities[g1],
				     g2,helicities[g2]).
	  dot(qqbarRightCurrent(k,helicities[k],
				l,helicities[l]));
      }
      if ( ijR != 0. && klL != 0. ) {
	diagram +=
	  ijR*klL*qqbarggRightCurrent(i,helicities[i],
				      j,helicities[j],
				      g1,helicities[g1],
				      g2,helicities[g2]).
	  dot(qqbarLeftCurrent(k,helicities[k],
			       l,helicities[l]));
      }
      if ( ijR != 0. && klR != 0. ) {
	diagram +=
	  ijR*klR*qqbarggRightCurrent(i,helicities[i],
				      j,helicities[j],
				      g1,helicities[g1],
				      g2,helicities[g2]).
	  dot(qqbarRightCurrent(k,helicities[k],
				l,helicities[l]));
      }
    } else if ( c->ijLineEmissions.first < 0 && c->ijLineEmissions.second < 0 &&
                c->klLineEmissions.first > 0 && c->klLineEmissions.second > 0 ) {
      int g1 = c->klLineEmissions.first;
      int g2 = c->klLineEmissions.second;
      if ( ijL != 0. && klL != 0. ) {
	diagram +=
	  ijL*klL*qqbarLeftCurrent(i,helicities[i],
				   j,helicities[j]).
	  dot(qqbarggLeftCurrent(k,helicities[k],
				 l,helicities[l],
				 g1,helicities[g1],
				 g2,helicities[g2]));
      }
      if ( ijL != 0. && klR != 0. ) {
	diagram +=
	  ijL*klR*qqbarLeftCurrent(i,helicities[i],
				   j,helicities[j]).
	  dot(qqbarggRightCurrent(k,helicities[k],
				  l,helicities[l],
				  g1,helicities[g1],
				  g2,helicities[g2]));
      }
      if ( ijR != 0. && klL != 0. ) {
	diagram +=
	  ijR*klL*qqbarRightCurrent(i,helicities[i],
				    j,helicities[j]).
	  dot(qqbarggLeftCurrent(k,helicities[k],
				 l,helicities[l],
				 g1,helicities[g1],
				 g2,helicities[g2]));
      }
      if ( ijR != 0. && klR != 0. ) {
	diagram +=
	  ijR*klR*qqbarRightCurrent(i,helicities[i],
				    j,helicities[j]).
	  dot(qqbarggRightCurrent(k,helicities[k],
				  l,helicities[l],
				  g1,helicities[g1],
				  g2,helicities[g2]));
      }
    } else if ( c->ijLineEmissions.first > 0 && c->ijLineEmissions.second < 0 &&
		c->klLineEmissions.first > 0 && c->klLineEmissions.second < 0 ) {
      int g1 = c->ijLineEmissions.first;
      int g2 = c->klLineEmissions.first;
      if ( ijL != 0. && klL != 0. ) {
	diagram +=
	  ijL*klL*qqbargLeftCurrent(i,helicities[i],
				    j,helicities[j],
				    g1,helicities[g1]).
	  dot(qqbargLeftCurrent(k,helicities[k],
				l,helicities[l],
				g2,helicities[g2]));
      }
      if ( ijL != 0. && klR != 0. ) {
	diagram +=
	  ijL*klR*qqbargLeftCurrent(i,helicities[i],
				    j,helicities[j],
				    g1,helicities[g1]).
	  dot(qqbargRightCurrent(k,helicities[k],
				 l,helicities[l],
				 g2,helicities[g2]));
      }
      if ( ijR != 0. && klL != 0. ) {
	diagram +=
	  ijR*klL*qqbargRightCurrent(i,helicities[i],
				     j,helicities[j],
				     g1,helicities[g1]).
	  dot(qqbargLeftCurrent(k,helicities[k],
				l,helicities[l],
				g2,helicities[g2]));
      }
      if ( ijR != 0. && klR != 0. ) {
	diagram +=
	  ijR*klR*qqbargRightCurrent(i,helicities[i],
				     j,helicities[j],
				     g1,helicities[g1]).
	  dot(qqbargRightCurrent(k,helicities[k],
				 l,helicities[l],
				 g2,helicities[g2]));
      }
    } else assert(false);

    diagram *= bosonFactor(*c)*cx->second  * c->fermionSign;

    result += diagram;

    if ( cx->second > 0. )
      largeN += diagram;

  }

  largeN *= (4.*Constants::pi*SM().alphaS());
  return (4.*Constants::pi*SM().alphaS())*result;
  
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Amplitudehqqbarkkbargg::persistentOutput(PersistentOStream &) const {}

void Amplitudehqqbarkkbargg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Amplitudehqqbarkkbargg,AmplitudeBase>
  describeHJetsAmplitudehqqbarkkbargg("HJets::Amplitudehqqbarkkbargg", "HwMatchboxBuiltin.so HJets.so");

void Amplitudehqqbarkkbargg::Init() {

  static ClassDocumentation<Amplitudehqqbarkkbargg> documentation
    ("Helicity amplitudes for 0 -> h q qbar k kbar g g");

}

