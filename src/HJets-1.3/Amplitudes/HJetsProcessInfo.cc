// -*- C++ -*-
//
// HJetsProcessInfo.cc is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "HJetsProcessInfo.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace HJets;

void AmplitudeInfo::print(ostream& os) const {

  if ( ijLineEmissions.first >= 0 &&
       ijLineEmissions.second >= 0 ) {
    os << "     " 
       << ijLineEmissions.first 
       << "       " 
       << ijLineEmissions.second
       << "     \n"
       << "      \\     /      \n"
       << "       \\   /       \n"
       << "        \\ /        \n";
  }


  if ( ijLineEmissions.first >= 0 &&
       ijLineEmissions.second < 0 ) {
    os << "         " 
       << ijLineEmissions.first 
       << "         \n"
       << "         |         \n"
       << "         |         \n"
       << "         |         \n";
  }

  if ( ijLineEmissions.first >= 0 ||
       ijLineEmissions.second >= 0 ) {
    os << ijLine.first 
       << "---<----0----<---"
       << ijLine.second << "\n";
  } else {
    os << ijLine.first 
       << "---<---------<---"
       << ijLine.second << "\n";
  }

  os << "         \\ " << boson->PDGName() << "\n"
     << "         /         \n"
     << "         \\........0\n"
     << "         /         \n"
     << "         \\         \n";

  if ( klLineEmissions.first >= 0 ||
       klLineEmissions.second >= 0 ) {
    os << klLine.first 
       << "---<----0----<---"
       << klLine.second << "\n";
  } else {
    os << klLine.first 
       << "---<---------<---"
       << klLine.second << "\n";
  }

  if ( klLineEmissions.first >= 0 &&
       klLineEmissions.second >= 0 ) {
    os << "        / \\        \n"
       << "       /   \\       \n"
       << "      /     \\      \n"
       << "     " 
       << klLineEmissions.first 
       << "       " 
       << klLineEmissions.second
       << "     \n";
  }

  if ( klLineEmissions.first >= 0 &&
       klLineEmissions.second < 0 ) {
    os << "         |         \n"
       << "         |         \n"
       << "         |         \n"
       << "         " 
       << klLineEmissions.first 
       << "         \n";
  }

  os << "x " << fermionSign << "\n\n";

}

bool isZ(long a, long b) {
  return a + b == 0 && abs(a) < 6;
}

bool isZPair(long a, long b,
	     long c, long d) {
  return isZ(a,b) && isZ(c,d);
}

bool isWplus(long a, long b) {
  if ( abs(a)%2 != 0 )
    swap(a,b);
  if ( a < 0 || b > 0 )
    return false;
  return
    a == ParticleID::u && b == ParticleID::dbar ||
    a == ParticleID::c && b == ParticleID::sbar;
}

bool isWminus(long a, long b) {
  return isWplus(-a,-b);
}

bool isWPair(long a, long b,
	     long c, long d) {
  return 
    (isWplus(a,b) && isWminus(c,d)) ||
    (isWminus(a,b) && isWplus(c,d));
}

vector<AmplitudeInfo> HJetsProcessInfo::configurations(const cPDVector& proc,
						       const StandardModelBase& sm,
						       bool complexMassScheme,
						       double kappaZ,
						       double kappaW) {

  /*
  cerr << "check if can handle ";
  for ( cPDVector::const_iterator p = proc.begin(); p != proc.end(); ++p )
    cerr << (**p).PDGName() << " ";
  cerr << "\n" << flush;
  */

  vector<AmplitudeInfo> res;

  cPDVector xproc = proc;

  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();

  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();

  cPDVector::iterator h = xproc.begin();
  for ( ; h != xproc.end(); ++h )
    if ( (**h).id() == ParticleID::h0 )
      break;

  if ( h == xproc.end() )
    return res;

  map<vector<pair<int,cPDPtr> >, vector<pair<int,cPDPtr> > > coreAndEmissions;

  vector<pair<int,cPDPtr> > proto;
  for ( size_t k = 0; k < xproc.size(); ++k ) {
    if ( xproc[k]->id() == ParticleID::h0 )
      continue;
    proto.push_back(pair<int,cPDPtr>(k,xproc[k]));
  }

  if ( proto.size() == 4 ) {
    coreAndEmissions.insert(make_pair(proto,vector<pair<int,cPDPtr> >()));
  }

  if ( proto.size() > 4 ) {

    bool gotGlue = false;

    for ( vector<pair<int,cPDPtr> >::const_iterator p = proto.begin();
	  p != proto.end(); ++p )
      if ( p->second->id() == ParticleID::g ) {
	gotGlue = true;
	break;
      }

    if ( gotGlue ) {

      vector<pair<int,cPDPtr> > core;
      vector<pair<int,cPDPtr> > emissions;
      for ( vector<pair<int,cPDPtr> >::const_iterator p = proto.begin();
	    p != proto.end(); ++p )
	if ( p->second->id() == ParticleID::g ) {
	  emissions.push_back(*p);
	} else {
	  core.push_back(*p);
	}

      if ( core.size() != 4 || emissions.size() > 2 )
	return res;

      coreAndEmissions.insert(make_pair(core,emissions));

    } else {

      if ( proto.size() != 6 )
	return res;

      for ( vector<pair<int,cPDPtr> >::const_iterator i = proto.begin();
	    i != proto.end(); ++i ) {
	vector<pair<int,cPDPtr> >::const_iterator j = i; ++j;
	for ( ; j != proto.end(); ++j ) {
	  if ( abs(i->second->id()) < 6 && abs(j->second->id()) < 6 &&
	       i->second->id() + j->second->id() == 0 ) {
	    vector<pair<int,cPDPtr> > core;
	    vector<pair<int,cPDPtr> > emissions;
	    for ( vector<pair<int,cPDPtr> >::const_iterator k = proto.begin();
		  k != proto.end(); ++k )
	      if ( k != i && k != j )
		core.push_back(*k);
	    emissions.push_back(*i);
	    emissions.push_back(*j);
	    coreAndEmissions.insert(make_pair(core,emissions));
	  }
	}
      }

    }

  }

  Energy MZ = sm.getParticleData(ParticleID::Z0)->hardProcessMass();
  Energy GZ = sm.getParticleData(ParticleID::Z0)->hardProcessWidth();
  Energy MW = sm.getParticleData(ParticleID::Wplus)->hardProcessMass();
  Energy GW = sm.getParticleData(ParticleID::Wplus)->hardProcessWidth();

  
  complex<double> alphaQED = sm.alphaEMMZ(); 
  complex<double> c2tw = complex<double>(1.,0.) - sm.sin2ThetaW();
  complex<double> s2tw = sm.sin2ThetaW();

  if ( complexMassScheme ) {
    c2tw = complex<Energy2>(sqr(MW),-MW*GW)/complex<Energy2>(sqr(MZ),-MZ*GZ);
    s2tw = complex<double>(1.,0.) - c2tw;
    double re1 = sm.fermiConstant()*sqr(MW);
    double im1 = sm.fermiConstant()*MW*GW;
    alphaQED = sqrt(2.)*complex<double>(re1,-im1)*s2tw/Constants::pi; 
  }

  complex<double> stw = sqrt(s2tw);
  complex<double> ctw = sqrt(c2tw);
  
  complex<double> gem = sqrt(4.*Constants::pi*alphaQED);
 

  map<vector<pair<int,cPDPtr> >, vector<pair<int,cPDPtr> > > coresSorted;
  for ( map<vector<pair<int,cPDPtr> >, vector<pair<int,cPDPtr> > >::const_iterator c =
	  coreAndEmissions.begin(); c != coreAndEmissions.end(); ++c ) {


    bool valid = true;
    vector<pair<int,cPDPtr> > coreSorted(4,pair<int,cPDPtr>(-1,cPDPtr()));
    for ( size_t k = 0; k < 4; ++k ) {
      if ( c->first[k].second->id() > 0 )
	if ( coreSorted[0].first < 0 ) {
	  coreSorted[0].first = c->first[k].first;
	  coreSorted[0].second = c->first[k].second;
	} else if ( coreSorted[2].first < 0 ) {
	  coreSorted[2].first = c->first[k].first;
	  coreSorted[2].second = c->first[k].second;
	} else {
	  valid = false;
	}
      if ( c->first[k].second->id() < 0 )
	if ( coreSorted[1].first < 0 ) {
	  coreSorted[1].first = c->first[k].first;
	  coreSorted[1].second = c->first[k].second;
	} else if ( coreSorted[3].first < 0 ) {
	  coreSorted[3].first = c->first[k].first;
	  coreSorted[3].second = c->first[k].second;
	} else {
	  valid = false;
	}
    }

    if ( !valid )
      continue;

    if ( coreSorted[0].first > coreSorted[2].first ) {
      swap(coreSorted[0],coreSorted[2]);
      swap(coreSorted[1],coreSorted[3]);
    }

    vector<pair<int,cPDPtr> > emissionsSorted = c->second;
    if ( emissionsSorted.size() == 2 ) {
      if ( emissionsSorted[0].second->id() < emissionsSorted[1].second->id() )
	swap(emissionsSorted[0],emissionsSorted[1]);
      if ( emissionsSorted[0].second->id() == emissionsSorted[1].second->id() ) {
	assert(emissionsSorted[0].second->id() == ParticleID::g);
	if ( emissionsSorted[0].first > emissionsSorted[1].first )
	  swap(emissionsSorted[0],emissionsSorted[1]);
      }
    }
    coresSorted.insert(make_pair(coreSorted,emissionsSorted));
  }
  coreAndEmissions = coresSorted;

  /*
  cerr << "cores sorted are:\n";
  for ( map<vector<pair<int,cPDPtr> >, vector<pair<int,cPDPtr> > >::const_iterator c =
	  coreAndEmissions.begin(); c != coreAndEmissions.end(); ++c ) {
    for ( vector<pair<int,cPDPtr> >::const_iterator p = c->first.begin();
	  p != c->first.end(); ++p ) {
      cerr << p->first << " " << p->second->PDGName() << " ";
    }
    cerr << " -- ";
    for ( vector<pair<int,cPDPtr> >::const_iterator p = c->second.begin();
	  p != c->second.end(); ++p ) {
      cerr << p->first << " " << p->second->PDGName() << " ";
    }
    cerr << "\n";
  }
  */

  for ( map<vector<pair<int,cPDPtr> >, vector<pair<int,cPDPtr> > >::const_iterator c =
	  coreAndEmissions.begin(); c != coreAndEmissions.end(); ++c ) {

    vector<AmplitudeInfo> cores;
    AmplitudeInfo core;
    core.ijLineEmissions.first = -1;
    core.ijLineEmissions.second = -1;
    core.klLineEmissions.first = -1;
    core.klLineEmissions.second = -1;
    core.qqbarEmitted = false;
    core.fermionSign = 1.;
    if ( c->second.size() == 2 )
      if ( abs(c->second[0].second->id()) < 6 )
	core.qqbarEmitted = true;

    if ( isZPair(c->first[0].second->id(),c->first[1].second->id(),
		 c->first[2].second->id(),c->first[3].second->id()) ||
	 isZPair(c->first[0].second->id(),c->first[3].second->id(),
		 c->first[1].second->id(),c->first[2].second->id()) ) {

      double Q; double I3;
      core.boson = sm.getParticleData(ParticleID::Z0);
      core.bosonMass = MZ;
      core.bosonWidth = GZ;
      core.higgsCoupling = kappaZ*gem*MZ/(stw*ctw);
      if( complexMassScheme){
	complex<Energy2> MZ2c = complex<Energy2>(sqr(MZ),-MZ*GZ);
	complex<Energy> MZc = GeV*sqrt(MZ2c/GeV2);
	core.higgsCoupling = kappaZ*gem*MZc/(stw*ctw);
      }
      
      Q = c->first[0].second->iCharge()/3.;
      I3 = c->first[0].second->id() % 2 == 0 ? 0.5 : -0.5;
      core.ijRightCoupling = -gem*Q*stw/ctw;
      core.ijLeftCoupling = gem*(I3-s2tw*Q)/(stw*ctw);
      Q = c->first[2].second->iCharge()/3.;
      I3 = c->first[2].second->id() % 2 == 0 ? 0.5 : -0.5;
      core.klRightCoupling = -gem*Q*stw/ctw;
      core.klLeftCoupling = gem*(I3-s2tw*Q)/(stw*ctw);

      if ( isZPair(c->first[0].second->id(),c->first[1].second->id(),
		   c->first[2].second->id(),c->first[3].second->id()) ) {
	core.ijLine.first = c->first[0].first;
	core.ijLine.second = c->first[1].first;
	core.klLine.first = c->first[2].first;
	core.klLine.second = c->first[3].first;
	cores.push_back(core);
      }

      if ( isZPair(c->first[0].second->id(),c->first[3].second->id(),
		   c->first[1].second->id(),c->first[2].second->id()) ) {
	core.ijLine.first = c->first[0].first;
	core.ijLine.second = c->first[3].first;
	core.klLine.first = c->first[2].first;
	core.klLine.second = c->first[1].first;
	cores.push_back(core);
      }

    }

    if ( isWPair(c->first[0].second->id(),c->first[1].second->id(),
		 c->first[2].second->id(),c->first[3].second->id()) ||
	 isWPair(c->first[0].second->id(),c->first[3].second->id(),
		 c->first[1].second->id(),c->first[2].second->id()) ) {

      core.bosonMass = MW;
      core.bosonWidth = GW;
      core.higgsCoupling = kappaW*gem*MW/stw;
      
      if( complexMassScheme){
	complex<Energy2> MW2c = complex<Energy2>(sqr(MW),-MW*GW);
	complex<Energy> MWc = GeV*sqrt(MW2c/GeV2);
	core.higgsCoupling = kappaW*gem*MWc/stw;
      }

      core.ijLeftCoupling = gem/(sqrt(2.)*stw);
      core.ijRightCoupling = 0.;
      core.klLeftCoupling = gem/(sqrt(2.)*stw);
      core.klRightCoupling = 0.;

      if ( isWPair(c->first[0].second->id(),c->first[1].second->id(),
		   c->first[2].second->id(),c->first[3].second->id()) ) {
	if ( isWplus(c->first[0].second->id(),c->first[1].second->id()) ) {
	  core.boson = sm.getParticleData(ParticleID::Wplus);
	} else {
	  core.boson = sm.getParticleData(ParticleID::Wminus);
	}
	core.ijLine.first = c->first[0].first;
	core.ijLine.second = c->first[1].first;
	core.klLine.first = c->first[2].first;
	core.klLine.second = c->first[3].first;
	cores.push_back(core);
      }

      if ( isWPair(c->first[0].second->id(),c->first[3].second->id(),
		   c->first[1].second->id(),c->first[2].second->id()) ) {
	if ( isWplus(c->first[0].second->id(),c->first[3].second->id()) ) {
	  core.boson = sm.getParticleData(ParticleID::Wplus);
	} else {
	  core.boson = sm.getParticleData(ParticleID::Wminus);
	}
	core.ijLine.first = c->first[0].first;
	core.ijLine.second = c->first[3].first;
	core.klLine.first = c->first[2].first;
	core.klLine.second = c->first[1].first;
	cores.push_back(core);
      }

    }

    for ( vector<AmplitudeInfo>::const_iterator conf = cores.begin();
	  conf != cores.end(); ++conf ) {

      AmplitudeInfo full = *conf;

      if ( c->second.size() == 0 ) {
	res.push_back(full);
      }
      
      if ( c->second.size() == 1 ) 
	if ( c->second[0].second->id() == ParticleID::g ) {

	  full.ijLineEmissions.first = c->second[0].first;
	  res.push_back(full);
	  full.ijLineEmissions.first = -1;
	  full.klLineEmissions.first = c->second[0].first;
	  res.push_back(full);
	  full.klLineEmissions.first = -1;

	}

      if ( c->second.size() == 2 ) 
	if ( c->second[0].second->id() == ParticleID::g &&
	     c->second[1].second->id() == ParticleID::g ) {

	  full.ijLineEmissions.first = c->second[0].first;
	  full.ijLineEmissions.second = c->second[1].first;
	  res.push_back(full);
	  swap(full.ijLineEmissions.first,full.ijLineEmissions.second);
	  res.push_back(full);
	  full.ijLineEmissions.first = -1;
	  full.ijLineEmissions.second = -1;

	  full.klLineEmissions.first = c->second[0].first;
	  full.klLineEmissions.second = c->second[1].first;
	  res.push_back(full);
	  swap(full.klLineEmissions.first,full.klLineEmissions.second);
	  res.push_back(full);
	  full.klLineEmissions.first = -1;
	  full.klLineEmissions.second = -1;

	  full.ijLineEmissions.first = c->second[0].first;
	  full.klLineEmissions.first = c->second[1].first;
	  res.push_back(full);
	  swap(full.ijLineEmissions.first,full.klLineEmissions.first);
	  res.push_back(full);
	  full.ijLineEmissions.first = -1;
	  full.klLineEmissions.first = -1;

	}

      if ( c->second.size() == 2 ) 
	if ( abs(c->second[0].second->id()) < 6 &&
	     c->second[0].second->id() + c->second[1].second->id() == 0 ) {

	  full.ijLineEmissions.first = c->second[0].first;
	  full.ijLineEmissions.second = c->second[1].first;
	  res.push_back(full);
	  full.ijLineEmissions.first = -1;
	  full.ijLineEmissions.second = -1;

	  full.klLineEmissions.first = c->second[0].first;
	  full.klLineEmissions.second = c->second[1].first;
	  res.push_back(full);
	  full.klLineEmissions.first = -1;
	  full.klLineEmissions.second = -1;

	}

    }

  }

  return res;

}


bool HJetsProcessInfo::canHandle(const PDVector& proc,
				 const StandardModelBase& sm) {

  cPDVector cproc;
  copy(proc.begin(), proc.end(), back_inserter(cproc));

  vector<AmplitudeInfo> confs = configurations(cproc,sm);

  return !confs.empty();

}

vector<AmplitudeInfo> HJetsProcessInfo::getConfigurations(const cPDVector& proc,
							  const vector<int>& crossed2proc,
							  const StandardModelBase& sm,
							  bool complexMassScheme,
							  double kappaZ,
							  double kappaW) {

  vector<AmplitudeInfo> procConfigs = configurations(proc,sm,complexMassScheme,kappaZ,kappaW);

  /*
  cerr << "configurations before crossing\n\n";

  for ( vector<AmplitudeInfo>::const_iterator c = procConfigs.begin();
	c != procConfigs.end(); ++c ) {
    c->print(cerr);
    cerr << "\n";
  }
  */

  vector<AmplitudeInfo> ampConfigs;

  map<int,int> proc2crossed;
  for ( size_t crossedx = 0; crossedx < crossed2proc.size(); ++crossedx )
    proc2crossed[crossed2proc[crossedx]] = crossedx;
  proc2crossed[-1] = -1;

  for ( vector<AmplitudeInfo>::const_iterator c = procConfigs.begin();
	c != procConfigs.end(); ++c ) {
    AmplitudeInfo crossed = *c;
    crossed.ijLine.first = proc2crossed[crossed.ijLine.first];
    crossed.ijLine.second = proc2crossed[crossed.ijLine.second];
    crossed.klLine.first = proc2crossed[crossed.klLine.first];
    crossed.klLine.second = proc2crossed[crossed.klLine.second];
    crossed.ijLineEmissions.first = proc2crossed[crossed.ijLineEmissions.first];
    crossed.ijLineEmissions.second = proc2crossed[crossed.ijLineEmissions.second];
    crossed.klLineEmissions.first = proc2crossed[crossed.klLineEmissions.first];
    crossed.klLineEmissions.second = proc2crossed[crossed.klLineEmissions.second];

    // work out the fermion signs

    set<pair<int,int> > ordering;

    ordering.insert(crossed.ijLine);
    ordering.insert(crossed.klLine);

    if ( crossed.qqbarEmitted ) {
      if ( crossed.ijLineEmissions.first >= 0 &&
	   crossed.ijLineEmissions.second >= 0 )
	ordering.insert(crossed.ijLineEmissions);
      else if ( crossed.klLineEmissions.first >= 0 &&
		crossed.klLineEmissions.second >= 0 )
	ordering.insert(crossed.klLineEmissions);
      else assert(false);
    }

    set<pair<int,int> >::iterator first = ordering.begin();
    set<pair<int,int> >::iterator second = first; ++second;
    set<pair<int,int> >::iterator third = ordering.end();
    if ( crossed.qqbarEmitted ) {
      third = second; ++third;
    }

    if ( !crossed.qqbarEmitted ) {

      if ( *first == pair<int,int>(1,2) && 
	   *second == pair<int,int>(3,4) ) {
	crossed.fermionSign = 1.;
      } else if ( *first == pair<int,int>(1,4) && 
		  *second == pair<int,int>(3,2) ) {
	crossed.fermionSign = -1.;
      } else assert(false);

    } else {

      if ( *first == pair<int,int>(1,2) && 
	   *second == pair<int,int>(3,4) &&
	   *third == pair<int,int>(5,6) ) {
	crossed.fermionSign = 1.;
      } else if ( *first == pair<int,int>(1,4) && 
		  *second == pair<int,int>(3,2) &&
		  *third == pair<int,int>(5,6) ) {
	crossed.fermionSign = -1.;
      } else if ( *first == pair<int,int>(1,6) && 
		  *second == pair<int,int>(3,4) &&
		  *third == pair<int,int>(5,2) ) {
	crossed.fermionSign = -1.;
      } else if ( *first == pair<int,int>(1,2) && 
		  *second == pair<int,int>(3,6) &&
		  *third == pair<int,int>(5,4) ) {
	crossed.fermionSign = -1.;
      } else if ( *first == pair<int,int>(1,4) && 
		  *second == pair<int,int>(3,6) &&
		  *third == pair<int,int>(5,2) ) {
	crossed.fermionSign = 1.;
      } else if ( *first == pair<int,int>(1,6) && 
		  *second == pair<int,int>(3,2) &&
		  *third == pair<int,int>(5,4) ) {
	crossed.fermionSign = 1.;
      } else assert(false);
      
    }

    ampConfigs.push_back(crossed);

  }

  return ampConfigs;

}
