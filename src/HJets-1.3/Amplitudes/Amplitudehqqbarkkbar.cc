// -*- C++ -*-
//
// Amplitudehqqbarkkbar.cc is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Amplitudehqqbarkkbar class.
//

#include "Amplitudehqqbarkkbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

extern "C"
{
  void qqhqqvbf_(double PBAR1[4], double PBAR2[4], double PBAR3[4],
		 double PBAR4[4], double PBAR5[4],
		 int SIGN[5], double ReCL1[2], double ReCR1[2],
		 double ReCL3[2], double ReCR3[2],
		 double ImCL1[2], double ImCR1[2],
		 double ImCL3[2], double ImCR3[2],
		 double Mass[2], double Width[2], int IsVaW[2],
		 int& divMax,
		 double ReGHVV[2], double ImGHVV[2], double& scale,
		 double& born, double virt[3]);
}

inline void F77_qqhqqvbf(double PBAR1[4], double PBAR2[4], double PBAR3[4],
		 double PBAR4[4], double PBAR5[4],
		 int SIGN[5], double ReCL1[2], double ReCR1[2],
		 double ReCL3[2], double ReCR3[2],
		 double ImCL1[2], double ImCR1[2],
		 double ImCL3[2], double ImCR3[2],
		 double Mass[2], double Width[2], int IsVaW[2], int& divMax,
		 double ReGHVV[2], double ImGHVV[2], double& scale,
		 double& born, double virt[3]){
  qqhqqvbf_(PBAR1,PBAR2,PBAR3,PBAR4,PBAR5,SIGN,ReCL1,ReCR1,ReCL3,ReCR3,
          ImCL1,ImCR1,ImCL3,ImCR3,Mass,Width,IsVaW,divMax,ReGHVV,ImGHVV,scale,born,virt);
  return;
}



using namespace HJets;

Amplitudehqqbarkkbar::Amplitudehqqbarkkbar() {}

Amplitudehqqbarkkbar::~Amplitudehqqbarkkbar() {}

IBPtr Amplitudehqqbarkkbar::clone() const {
  return new_ptr(*this);
}

IBPtr Amplitudehqqbarkkbar::fullclone() const {
  return new_ptr(*this);
}

void Amplitudehqqbarkkbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  // inform the amplitude cache how many legs we are handling
  nPoints(5);

  // everything is expressed in dimensionless quantities, referring to
  // sqrt(shat) of the collission as the mass unit
  amplitudeScale(sqrt(lastSHat()));

  // the higgs momentum
  Lorentz5Momentum ph = amplitudeMomentum(0);
  momentum(0,ph,true,ph.mass());

  // the parton momenta
  for ( size_t k = 1; k < 5; ++k )
    momentum(k,amplitudeMomentum(k),true);

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex Amplitudehqqbarkkbar::evaluate(size_t colourTensor, 
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

    if ( c->ijLineEmissions.first < 0 && c->ijLineEmissions.second < 0 &&
	 c->klLineEmissions.first < 0 && c->klLineEmissions.second < 0 ) {
      if ( ijL != 0. && klL != 0. ) {
	diagram +=
	  ijL*klL*qqbarLeftCurrent(i,helicities[i],
				   j,helicities[j]).
	  dot(qqbarLeftCurrent(k,helicities[k],
			       l,helicities[l]));
      }
      if ( ijL != 0. && klR != 0. ) {
	diagram +=
	  ijL*klR*qqbarLeftCurrent(i,helicities[i],
				   j,helicities[j]).
	  dot(qqbarRightCurrent(k,helicities[k],
				l,helicities[l]));
      }
      if ( ijR != 0. && klL != 0. ) {
	diagram +=
	  ijR*klL*qqbarRightCurrent(i,helicities[i],
				    j,helicities[j]).
	  dot(qqbarLeftCurrent(k,helicities[k],
			       l,helicities[l]));
      }
      if ( ijR != 0. && klR != 0. ) {
	diagram +=
	  ijR*klR*qqbarRightCurrent(i,helicities[i],
				    j,helicities[j]).
	  dot(qqbarRightCurrent(k,helicities[k],
				l,helicities[l]));
	
      }
    } else assert(false);

    diagram *= bosonFactor(*c)*cx->second * c->fermionSign;

    result += diagram;

    if ( cx->second > 0. )
      largeN += diagram;

  }

  return result;

}


inline void convert(const Lorentz5Momentum& p, double* res) {
  res[0] = p.t()/GeV;
  res[1] = p.x()/GeV;
  res[2] = p.y()/GeV;
  res[3] = p.z()/GeV;
}

double Amplitudehqqbarkkbar::oneLoopInterference() const {

  if ( !calculateOneLoopInterference() )
    return lastOneLoopInterference();
  
  const map<int, pair<int, double> >& transmap = virtualInfo();

  double pbar1[4];
  double pbar2[4];
  double pbar3[4];
  double pbar4[4];
  double pbar5[4];

  convert(meMomenta()[transmap.find(1)->second.first],pbar1);
  convert(meMomenta()[transmap.find(2)->second.first],pbar2);
  convert(meMomenta()[transmap.find(3)->second.first],pbar3);
  convert(meMomenta()[transmap.find(4)->second.first],pbar4);
  convert(meMomenta()[transmap.find(5)->second.first],pbar5);



  int sign[5];
  for ( size_t k = 0; k < 5; ++k )
    sign[k] = transmap.find(k+1)->second.second;

  double Recl1[2];
  double Recr1[2];
  double Recl3[2];
  double Recr3[2];

  double Reghvv[2];
  double Imcl1[2];
  double Imcr1[2];
  double Imcl3[2];
  double Imcr3[2];

  double Imghvv[2];

  double mass[2];
  double width[2];

  int isVaW[2];

  getCouplings(Recl1,Recr1,Recl3,Recr3,
               Imcl1,Imcr1,Imcl3,Imcr3,Reghvv,Imghvv,mass,width,isVaW);

  //int NF = nLight();
  int divMax = 2;

  double born;
  double virt[3];
  double scale = mu2()/GeV2;

  F77_qqhqqvbf (pbar1,pbar2,pbar3,pbar4,pbar5,sign,Recl1,Recr1,Recl3,Recr3,
		Imcl1,Imcr1,Imcl3,Imcr3,
		mass,width,isVaW,divMax,Reghvv,Imghvv,scale,born,virt);

  double res = (SM().alphaS()/(2.*Constants::pi))*virt[0]*(lastSHat()/GeV2);

 
 
#ifdef AMPVERBOSE
  const cPDVector& proc = mePartonData();
  const cPDVector& proc1 = amplitudePartonData();
  map<int, pair<int, double> > XMap = virtualInfos().find(proc)->second;

  generator()->log() << "c2="<< virt[2]/born << "\n";
  generator()->log() << "c1="<< virt[1]/born << "\n";
  generator()->log() << "c0="<< virt[0]/born << "\n";
  born*=lastSHat()/GeV2;

  generator()->log() << "--------------------- \n";
  generator()->log() << "rb=" << born/me2() << "  Base(TF) Diagram:"
		     << XMap.find(1)->second.first << (XMap.find(1)->second.second > 0 ? "+" : "-") << proc[XMap.find(1)->second.first]->PDGName() << " "
		     << XMap.find(3)->second.first << (XMap.find(3)->second.second > 0 ? "+" : "-") << proc[XMap.find(3)->second.first]->PDGName() << " -> "
		     << XMap.find(2)->second.first << (XMap.find(2)->second.second > 0 ? "+" : "-") << proc[XMap.find(2)->second.first]->PDGName() << " "
		     << XMap.find(4)->second.first << (XMap.find(4)->second.second > 0 ? "+" : "-") << proc[XMap.find(4)->second.first]->PDGName() << " "
		     << XMap.find(5)->second.first << (XMap.find(5)->second.second > 0 ? "+" : "-") << proc[XMap.find(5)->second.first]->PDGName() << "\n";

  generator()->log() << "rb=" << born/me2() << "  Base(SP) Diagram:"
		     << crossingMap()[0] << " " << proc1[0]->PDGName() << " "
		     << crossingMap()[1] << " " << proc1[1]->PDGName() << " "
		     << crossingMap()[2] << " " << proc1[2]->PDGName() << " "
		     << crossingMap()[3] << " " << proc1[3]->PDGName() << " "
		     << crossingMap()[4] << " " << proc1[4]->PDGName() << "\n";
  generator()->log() << "OneisNC= " << topologyOneIsNC()
		     << " OneisCC= " << topologyOneIsCC()
		     << " TwoisNC= " << topologyTwoIsNC()
		     << " TwoisCC= " << topologyTwoIsCC()
		     << "\n";
 

  generator()->log() << "sign[0]="<< sign[0]<< "\n";
  generator()->log() << "sign[1]="<< sign[1]<< "\n";
  generator()->log() << "sign[2]="<< sign[2]<< "\n";
  generator()->log() << "sign[3]="<< sign[3]<< "\n";
  generator()->log() << "sign[4]="<< sign[4]<< "\n";

  generator()->log() << "Recl1[0]="<< Recl1[0] << "  Recl3[0]=" << Recl3[0] << "\n";
  generator()->log() << "Recr1[0]="<< Recr1[0] << "  Recr3[0]=" << Recr3[0] << "\n";
  generator()->log() << "Imcl1[0]="<< Imcl1[0] << "  Imcl3[0]=" << Imcl3[0] << "\n";
  generator()->log() << "Imcr1[0]="<< Imcr1[0] << "  Imcr3[0]=" << Imcr3[0] << "\n";

  generator()->log() << "Mass[0]="<< mass[0] << "\n";

  generator()->log() << "Recl1[1]="<< Recl1[1] << "  Recl3[0]=" << Recl3[1] << "\n";
  generator()->log() << "Recr1[1]="<< Recr1[1] << "  Recr3[0]=" << Recr3[1] << "\n";
  generator()->log() << "Imcl1[1]="<< Imcl1[1] << "  Imcl3[0]=" << Imcl3[1] << "\n";
  generator()->log() << "Imcr1[1]="<< Imcr1[1] << "  Imcr3[0]=" << Imcr3[1] << "\n";

  generator()->log() << "Mass[1]="<< mass[1] << "\n";


  generator()->log() << "born=" << born << "\n";
 
  generator()->log() << "me2=" << me2() << "\n";
  
#endif 

  /*
  double mm = me2();
  born *= lastSHat()/GeV2;
  double delta = abs((mm-born)/(mm+born));
  if ( mm < 0.0 || born < 0.0 ) {
    cerr << "in ";
    for ( cPDVector::const_iterator p = mePartonData().begin();
	  p != mePartonData().end(); ++p )
      cerr << (**p).PDGName() << " ";
    cerr << setprecision(20);
    cerr << "\nfatal : mm = " << mm << " born = " << born << "\n" << flush;
  }
  if ( delta > 1.e-12 ) {
    cerr << "in ";
    for ( cPDVector::const_iterator p = mePartonData().begin();
	  p != mePartonData().end(); ++p )
      cerr << (**p).PDGName() << " ";
    cerr << setprecision(20);
    cerr << "\nproblem : mm = " << mm << " born = " << born 
	 << " -> delta = " << delta << "\n" << flush;
  }
  */

  lastOneLoopInterference(res);

  return lastOneLoopInterference();
  
}

double Amplitudehqqbarkkbar::oneLoopDoublePole() const {

  
  const map<int, pair<int, double> >& transmap = virtualInfo();

  double pbar1[4];
  double pbar2[4];
  double pbar3[4];
  double pbar4[4];
  double pbar5[4];

  convert(meMomenta()[transmap.find(1)->second.first],pbar1);
  convert(meMomenta()[transmap.find(2)->second.first],pbar2);
  convert(meMomenta()[transmap.find(3)->second.first],pbar3);
  convert(meMomenta()[transmap.find(4)->second.first],pbar4);
  convert(meMomenta()[transmap.find(5)->second.first],pbar5);

  

  int sign[5];
  for ( size_t k = 0; k < 5; ++k )
    sign[k] = transmap.find(k+1)->second.second;
  

  double Recl1[2];
  double Recr1[2];
  double Recl3[2];
  double Recr3[2];

  double Reghvv[2];
  double Imcl1[2];
  double Imcr1[2];
  double Imcl3[2];
  double Imcr3[2];

  double Imghvv[2];

  double mass[2];
  double width[2];

  int IsVaW[2];

  getCouplings(Recl1,Recr1,Recl3,Recr3,
               Imcl1,Imcr1,Imcl3,Imcr3,Reghvv,Imghvv,mass,width,IsVaW);

  //int NF = nLight();
  int divMax = 2;

  double born;
  double virt[3];
  double scale = mu2()/GeV2;

  F77_qqhqqvbf (pbar1,pbar2,pbar3,pbar4,pbar5,sign,Recl1,Recr1,Recl3,Recr3,
		Imcl1,Imcr1,Imcl3,Imcr3,mass,width,IsVaW,divMax,Reghvv,Imghvv,scale,born,virt);

  

  double res = (SM().alphaS()/(2.*Constants::pi))*virt[2]*(lastSHat()/GeV2);
  
#ifdef AMPVERBOSE
  generator()->log() << "double pole="<< res/4./3./3. << "\n";
#endif
  
  return res;

  
}

double Amplitudehqqbarkkbar::oneLoopSinglePole() const {

  
  const map<int, pair<int, double> >& transmap = virtualInfo();

  double pbar1[4];
  double pbar2[4];
  double pbar3[4];
  double pbar4[4];
  double pbar5[4];

  convert(meMomenta()[transmap.find(1)->second.first],pbar1);
  convert(meMomenta()[transmap.find(2)->second.first],pbar2);
  convert(meMomenta()[transmap.find(3)->second.first],pbar3);
  convert(meMomenta()[transmap.find(4)->second.first],pbar4);
  convert(meMomenta()[transmap.find(5)->second.first],pbar5);

  

  int sign[5];
  for ( size_t k = 0; k < 5; ++k )
    sign[k] = transmap.find(k+1)->second.second;


  double Recl1[2];
  double Recr1[2];
  double Recl3[2];
  double Recr3[2];

  double Reghvv[2];
  double Imcl1[2];
  double Imcr1[2];
  double Imcl3[2];
  double Imcr3[2];

  double Imghvv[2];

  double mass[2];
  double width[2];

  int IsVaW[2];

  getCouplings(Recl1,Recr1,Recl3,Recr3,
               Imcl1,Imcr1,Imcl3,Imcr3,Reghvv,Imghvv,mass,width,IsVaW);

  //int NF = nLight();
  int divMax = 2;

  double born;
  double virt[3];
  double scale = mu2()/GeV2;

  F77_qqhqqvbf (pbar1,pbar2,pbar3,pbar4,pbar5,sign,Recl1,Recr1,Recl3,Recr3,
		Imcl1,Imcr1,Imcl3,Imcr3,
		mass,width,IsVaW,divMax,Reghvv,Imghvv,scale,born,virt);

  

  double res = (SM().alphaS()/(2.*Constants::pi))*virt[1]*(lastSHat()/GeV2);
  
  
#ifdef AMPVERBOSE
  generator()->log() << "single pole="<< res/4./3./3. << "\n";
#endif

  return res;

  
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Amplitudehqqbarkkbar::persistentOutput(PersistentOStream &) const {}

void Amplitudehqqbarkkbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Amplitudehqqbarkkbar,AmplitudeBase>
  describeHJetsAmplitudehqqbarkkbar("HJets::Amplitudehqqbarkkbar", "HwMatchboxBuiltin.so HJets.so");

void Amplitudehqqbarkkbar::Init() {

  static ClassDocumentation<Amplitudehqqbarkkbar> documentation
    ("Helicity amplitudes for 0 -> h q qbar k kbar");

}

