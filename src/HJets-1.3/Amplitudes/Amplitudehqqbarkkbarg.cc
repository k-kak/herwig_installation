// -*- C++ -*-
//
// Amplitudehqqbarkkbarg.cc is a part of HJets
// Copyright (C) 2011-2012 
// Ken Arnold, Francisco Campanario, Terrance Figy, Simon Platzer and Malin Sjodahl
//
// HJets is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Amplitudehqqbarkkbarg class.
//

#include "Amplitudehqqbarkkbarg.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

extern "C"
{
  void qqhqqgvbf_(double PBAR1[4], double PBAR2[4], double PBAR3[4],
		  double PBAR4[4], double PBAR5[4], double PBAR6[4],
		  int SIGN[6], double ReCL1[2], double ReCR1[2],
                  double ReCL3[2], double ReCR3[2],
		  double ImCL1[2], double ImCR1[2],
                  double ImCL3[2], double ImCR3[2],
		  double Mass[2], double Width[2], int IsVaW[2],
		  int& divMax, int gaugegood[3], int gaugebad[3], double& gaugethres,
                  double& gaugeacc,
		  double ReGHVV[2], double ImGHVV[2], double& scale, double& NF,
		  double& born, double virt[3]);
}

inline void F77_qqhqqgvbf(double PBAR1[4], double PBAR2[4], double PBAR3[4],
		  double PBAR4[4], double PBAR5[4], double PBAR6[4],
		  int SIGN[6], double ReCL1[2], double ReCR1[2],
                  double ReCL3[2], double ReCR3[2],
		  double ImCL1[2], double ImCR1[2],
                  double ImCL3[2], double ImCR3[2],
		  double Mass[2], double Width[2], int IsVaW[2],
		  int& divMax, int gaugegood[3], int gaugebad[3], double& gaugethres,
		  double& gaugeacc,
		  double ReGHVV[2], double ImGHVV[2], double& scale, double& NF,
		  double& born, double virt[3]){
  qqhqqgvbf_(PBAR1,PBAR2,PBAR3,PBAR4,PBAR5,PBAR6,SIGN,ReCL1,ReCR1,ReCL3,ReCR3,
	     ImCL1,ImCR1,ImCL3,ImCR3,Mass,Width,IsVaW,divMax,
             gaugegood,gaugebad,gaugethres,gaugeacc,ReGHVV,ImGHVV,scale,NF,born,virt);
  return;
}

using namespace HJets;

struct GaugeStatistics {

  unsigned long gaugePassed[3];

  unsigned long gaugeFailed[3];

  GaugeStatistics();

  ~GaugeStatistics();

};

GaugeStatistics::GaugeStatistics() {
  for ( unsigned int i = 0; i < 3; ++i ) {
    gaugePassed[i] = 0;
    gaugeFailed[i] = 0;
  }
}

GaugeStatistics::~GaugeStatistics() {
  if ( (double)gaugePassed[0] > 0 ||
       (double)gaugePassed[1] > 0 ||
       (double)gaugePassed[2] > 0 ) {
    cerr << "\ngauge check summary\n"
	 << "--------------------------------------------------------------------------------\n";
    if ( gaugePassed[0] > 0 )
      cerr << 100*((double)gaugeFailed[0]/(double)gaugePassed[0]) << "% points failing box gauge checks\n";
    if ( gaugePassed[1] > 0 )
      cerr << 100*((double)gaugeFailed[1]/(double)gaugePassed[1]) << "% points failing abelian hexagon gauge checks\n";
    if ( gaugePassed[2] > 0 )
      cerr << 100*((double)gaugeFailed[2]/(double)gaugePassed[2]) << "% points failing non-abelian hexagon gauge checks\n";
    cerr  << "--------------------------------------------------------------------------------\n"
	  << flush;
  } else if ( (double)gaugeFailed[0] > 0 ||
	      (double)gaugeFailed[1] > 0 ||
	      (double)gaugeFailed[2] > 0 ) {
    cerr << "\ngauge check fatal: no points passed, some failed!\n" << flush;
  }
}

static GaugeStatistics theGaugeStatistics;

Amplitudehqqbarkkbarg::Amplitudehqqbarkkbarg()
  : theGaugeThreshold(0.1), theAccuracyCheck(false) {}




Amplitudehqqbarkkbarg::~Amplitudehqqbarkkbarg() {}

IBPtr Amplitudehqqbarkkbarg::clone() const {
  return new_ptr(*this);
}

IBPtr Amplitudehqqbarkkbarg::fullclone() const {
  return new_ptr(*this);
}

void Amplitudehqqbarkkbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  // inform the amplitude cache how many legs we are handling
  nPoints(6);

  // everything is expressed in dimensionless quantities, referring to
  // sqrt(shat) of the collission as the mass unit
  amplitudeScale(sqrt(lastSHat()));

  // the higgs momentum
  Lorentz5Momentum ph = amplitudeMomentum(0);
  momentum(0,ph,true,ph.mass());

  // the parton momenta
  for ( size_t k = 1; k < 6; ++k )
    momentum(k,amplitudeMomentum(k),true);

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex Amplitudehqqbarkkbarg::evaluate(size_t colourTensor, 
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

    if ( c->ijLineEmissions.first > 0 && c->ijLineEmissions.second < 0 &&
	 c->klLineEmissions.first < 0 && c->klLineEmissions.second < 0 ) {
      int g = c->ijLineEmissions.first;
      if ( ijL != 0. && klL != 0. ) {
	diagram +=
	  ijL*klL*qqbargLeftCurrent(i,helicities[i],
				    j,helicities[j],
				    g,helicities[g]).
	  dot(qqbarLeftCurrent(k,helicities[k],
			       l,helicities[l]));
      }
      if ( ijL != 0. && klR != 0. ) {
	diagram +=
	  ijL*klR*qqbargLeftCurrent(i,helicities[i],
				    j,helicities[j],
				    g,helicities[g]).
	  dot(qqbarRightCurrent(k,helicities[k],
				l,helicities[l]));
      }
      if ( ijR != 0. && klL != 0. ) {
	diagram +=
	  ijR*klL*qqbargRightCurrent(i,helicities[i],
				     j,helicities[j],
				     g,helicities[g]).
	  dot(qqbarLeftCurrent(k,helicities[k],
			       l,helicities[l]));
      }
      if ( ijR != 0. && klR != 0. ) {
	diagram +=
	  ijR*klR*qqbargRightCurrent(i,helicities[i],
				     j,helicities[j],
				     g,helicities[g]).
	  dot(qqbarRightCurrent(k,helicities[k],
				l,helicities[l]));
      }
    } else if ( c->ijLineEmissions.first < 0 && c->ijLineEmissions.second < 0 &&
		c->klLineEmissions.first > 0 && c->klLineEmissions.second < 0 ) {
      int g = c->klLineEmissions.first;
      if ( ijL != 0. && klL != 0. ) {
	diagram +=
	  ijL*klL*qqbarLeftCurrent(i,helicities[i],
				   j,helicities[j]).
	  dot(qqbargLeftCurrent(k,helicities[k],
				l,helicities[l],
				g,helicities[g]));
      }
      if ( ijL != 0. && klR != 0. ) {
	diagram +=
	  ijL*klR*qqbarLeftCurrent(i,helicities[i],
				   j,helicities[j]).
	  dot(qqbargRightCurrent(k,helicities[k],
				 l,helicities[l],
				 g,helicities[g]));
      }
      if ( ijR != 0. && klL != 0. ) {
	diagram +=
	  ijR*klL*qqbarRightCurrent(i,helicities[i],
				    j,helicities[j]).
	  dot(qqbargLeftCurrent(k,helicities[k],
				l,helicities[l],
				g,helicities[g]));
      }
      if ( ijR != 0. && klR != 0. ) {
	diagram +=
	  ijR*klR*qqbarRightCurrent(i,helicities[i],
				    j,helicities[j]).
	  dot(qqbargRightCurrent(k,helicities[k],
				 l,helicities[l],
				 g,helicities[g]));
      }
    } else assert(false);

    diagram *= bosonFactor(*c)*cx->second  * c->fermionSign;

    result += diagram;

    if ( cx->second > 0. )
      largeN += diagram;

  }

  largeN *= sqrt(4.*Constants::pi*SM().alphaS());
  return sqrt(4.*Constants::pi*SM().alphaS())*result;

}

inline void convert(const Lorentz5Momentum& p, double* res) {
  res[0] = p.t()/GeV;
  res[1] = p.x()/GeV;
  res[2] = p.y()/GeV;
  res[3] = p.z()/GeV;
}

double Amplitudehqqbarkkbarg::oneLoopInterference() const {

  if ( !calculateOneLoopInterference() )
    return lastOneLoopInterference();

  const map<int, pair<int, double> >& transmap = virtualInfo();

  double pbar1[4];
  double pbar2[4];
  double pbar3[4];
  double pbar4[4];
  double pbar5[4];
  double pbar6[4];

  convert(meMomenta()[transmap.find(1)->second.first],pbar1);
  convert(meMomenta()[transmap.find(3)->second.first],pbar3);
  convert(meMomenta()[transmap.find(2)->second.first],pbar2);
  convert(meMomenta()[transmap.find(4)->second.first],pbar4);
  convert(meMomenta()[transmap.find(5)->second.first],pbar5);
  convert(meMomenta()[transmap.find(6)->second.first],pbar6);

  int sign[6];
  for ( size_t k = 0; k < 6; ++k )
    sign[k] = (int) transmap.find(k+1)->second.second;

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
               Imcl1,Imcr1,Imcl3,Imcr3,Reghvv,Imghvv,
               mass,width,IsVaW);

  double NF = nLight();
  int gaugegood[3];
  int gaugebad[3];
  int divMax = 2;

  double born;
  double virt[3];
  double scale = mu2()/GeV2;

  double gt = theGaugeThreshold;

  double gacc(0.);

#ifdef AMPVERBOSE
  generator()->log()<<"scale="<<scale << "\n";
#endif

  F77_qqhqqgvbf (pbar1,pbar2,pbar3,pbar4,pbar5,pbar6,sign,Recl1,Recr1,Recl3,Recr3,
	     Imcl1,Imcr1,Imcl3,Imcr3,mass,width,IsVaW,divMax,gaugegood,gaugebad,gt,gacc,
             Reghvv,Imghvv,scale,NF,born,virt);

  for (int i = 0; i < 3; i++) {
    theGaugeStatistics.gaugePassed[i] += gaugegood[i];
    theGaugeStatistics.gaugeFailed[i] += gaugebad[i];
  }

  

  double res = 2.*sqr(SM().alphaS())*virt[0]*pow(lastSHat()/GeV2,2);

  if (gaugebad[0]+gaugebad[1]+gaugebad[2]) {
    throw Veto();
  }

  if( theAccuracyCheck ){
    generator()->log() << "checkacc "<< log10 (abs(gacc))<<"\n";
  }

#ifdef AMPVERBOSE
  const cPDVector& proc = mePartonData();
  map<int, pair<int, double> > XMap = virtualInfos().find(proc)->second;
  generator()->log() << setprecision(20);
  generator()->log() << "--------------------- \n";
  generator()->log() << "c2="<< virt[2]/born << "\n";
  generator()->log() << "c1="<< virt[1]/born << "\n";
  generator()->log() << "c0="<< virt[0]/born << "\n";
  double gs2(SM().alphaS()*4.*Constants::pi);
  born*=pow(lastSHat()/GeV2,2)*gs2;

  generator()->log() << "rb=" << born/me2() << "  Base Diagram:"
       << (XMap.find(1)->second.second > 0 ? "+" : "-") << proc[XMap.find(1)->second.first]->PDGName() << " "
       << (XMap.find(3)->second.second > 0 ? "+" : "-") << proc[XMap.find(3)->second.first]->PDGName() << " -> "
       << (XMap.find(2)->second.second > 0 ? "+" : "-") << proc[XMap.find(2)->second.first]->PDGName() << " "
       << (XMap.find(4)->second.second > 0 ? "+" : "-") << proc[XMap.find(4)->second.first]->PDGName() << " "
       << (XMap.find(5)->second.second > 0 ? "+" : "-") << proc[XMap.find(5)->second.first]->PDGName() << " "
       << (XMap.find(6)->second.second > 0 ? "+" : "-") << proc[XMap.find(6)->second.first]->PDGName() << "\n";
  generator()->log() << "OneisNC= " << topologyOneIsNC()
       << " OneisCC= " << topologyOneIsCC()
       << " TwoisNC= " << topologyTwoIsNC()
       << " TwoisCC= " << topologyTwoIsCC()
       << "\n";
 

  generator()->log() << "sign[0]="<< sign[0]<< "  xmap(1)="<< XMap.find(1)->second.second << "\n";
  generator()->log() << "sign[2]="<< sign[2]<< "  xmap(3)="<< XMap.find(3)->second.second << "\n";
  generator()->log() << "sign[1]="<< sign[1]<< "\n";
  generator()->log() << "sign[3]="<< sign[3]<< "\n";
  generator()->log() << "sign[4]="<< sign[4]<< "\n";
  generator()->log() << "sign[5]="<< sign[5]<< "\n";
  

  generator()->log() << "Recl1[0]="<< Recl1[0] << "  Recl3[0]=" << Recl3[0] << "\n";
  generator()->log() << "Recr1[0]="<< Recr1[0] << "  Recr3[0]=" << Recr3[0] << "\n";
  generator()->log() << "Imcl1[0]="<< Imcl1[0] << "  Imcl3[0]=" << Imcl3[0] << "\n";
  generator()->log() << "Imcr1[0]="<< Imcr1[0] << "  Imcr3[0]=" << Imcr3[0] << "\n";
  generator()->log() << "Reghvv[0]="<< Reghvv[0]<< " Imghvv[0]=" << Imghvv[0] << "\n";
  generator()->log() << "Mass[0]="<< mass[0] << "\n";

  generator()->log() << "Recl1[1]="<< Recl1[1] << "  Recl3[0]=" << Recl3[1] << "\n";
  generator()->log() << "Recr1[1]="<< Recr1[1] << "  Recr3[0]=" << Recr3[1] << "\n";
  generator()->log() << "Imcl1[1]="<< Imcl1[1] << "  Imcl3[0]=" << Imcl3[1] << "\n";
  generator()->log() << "Imcr1[1]="<< Imcr1[1] << "  Imcr3[0]=" << Imcr3[1] << "\n";

  generator()->log() << "Reghvv[1]="<< Reghvv[1] << " Imghvv[1]="<<Imghvv[1] << "\n";

  generator()->log() << "Mass[1]="<< mass[1] << "\n";


  generator()->log() << "born=" << born << "\n";
  generator()->log() << "me2=" << me2() << "\n";
  generator()->log() << "--------------------- \n";
#endif

  lastOneLoopInterference(res);

  return lastOneLoopInterference();

}


double Amplitudehqqbarkkbarg::oneLoopDoublePole() const {

  const map<int, pair<int, double> >& transmap = virtualInfo();

  double pbar1[4];
  double pbar2[4];
  double pbar3[4];
  double pbar4[4];
  double pbar5[4];
  double pbar6[4];

  convert(meMomenta()[transmap.find(1)->second.first],pbar1);
  convert(meMomenta()[transmap.find(3)->second.first],pbar3);
  convert(meMomenta()[transmap.find(2)->second.first],pbar2);
  convert(meMomenta()[transmap.find(4)->second.first],pbar4);
  convert(meMomenta()[transmap.find(5)->second.first],pbar5);
  convert(meMomenta()[transmap.find(6)->second.first],pbar6);

  int sign[6];
  for ( size_t k = 0; k < 6; ++k )
    sign[k] = (int) transmap.find(k+1)->second.second;
  
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
		Imcl1,Imcr1,Imcl3,Imcr3,Reghvv,Imghvv,
		mass,width,IsVaW);
  
   
   double NF = nLight();
   int gaugegood[3];
   int gaugebad[3];
   int divMax = 2;
   
   double born;
   double virt[3];
   double scale = mu2()/GeV2;

   double gt = theGaugeThreshold;
   double gacc(0.);
   
#ifdef AMPVERBOSE
  generator()->log()<<"scale="<<scale << "\n";
#endif

  F77_qqhqqgvbf (pbar1,pbar2,pbar3,pbar4,pbar5,pbar6,sign,Recl1,Recr1,Recl3,Recr3,
	     Imcl1,Imcr1,Imcl3,Imcr3,mass,width,IsVaW,divMax,gaugegood,gaugebad,
             gt,gacc,Reghvv,Imghvv,
	     scale,NF,born,virt);

 
  double res = 2.*sqr(SM().alphaS())*virt[2]*pow(lastSHat()/GeV2,2);

  if (gaugebad[0]+gaugebad[1]+gaugebad[2]) {
    throw Veto();
  }

  /*generator()->log() << "c2="<< virt[2]/born << "\n";
  generator()->log() << "c1="<< virt[1]/born << "\n";
  generator()->log() << "c0="<< virt[0]/born << "\n";
  */
  return res;

}

double Amplitudehqqbarkkbarg::oneLoopSinglePole() const {

  const map<int, pair<int, double> >& transmap = virtualInfo();

  double pbar1[4];
  double pbar2[4];
  double pbar3[4];
  double pbar4[4];
  double pbar5[4];
  double pbar6[4];

  convert(meMomenta()[transmap.find(1)->second.first],pbar1);
  convert(meMomenta()[transmap.find(3)->second.first],pbar3);
  convert(meMomenta()[transmap.find(2)->second.first],pbar2);
  convert(meMomenta()[transmap.find(4)->second.first],pbar4);
  convert(meMomenta()[transmap.find(5)->second.first],pbar5);
  convert(meMomenta()[transmap.find(6)->second.first],pbar6);

  int sign[6];
  for ( size_t k = 0; k < 6; ++k )
    sign[k] = (int) transmap.find(k+1)->second.second;

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
               Imcl1,Imcr1,Imcl3,Imcr3,Reghvv,Imghvv,
               mass,width,IsVaW);
  
   double NF = nLight();
   int gaugegood[3];
   int gaugebad[3];
   int divMax = 2;
   
   double born;
   double virt[3];
   double scale = mu2()/GeV2;

   double gt = theGaugeThreshold;
   double gacc(0.);
   
#ifdef AMPVERBOSE
  generator()->log()<<"scale="<<scale << "\n";
#endif

  F77_qqhqqgvbf (pbar1,pbar2,pbar3,pbar4,pbar5,pbar6,sign,Recl1,Recr1,Recl3,Recr3,
	     Imcl1,Imcr1,Imcl3,Imcr3,
	     mass,width,IsVaW,divMax,gaugegood,gaugebad,gt,gacc,Reghvv,Imghvv,
	     scale,NF,born,virt);

  double res = 2.*sqr(SM().alphaS())*virt[1]*pow(lastSHat()/GeV2,2);
  
  if (gaugebad[0]+gaugebad[1]+gaugebad[2]) {
    throw Veto();
  }

  /*
  generator()->log() << "NLight= "<<NF<< "\n";
  generator()->log() << "c2="<< virt[2]/born << "\n";
  generator()->log() << "c1="<< virt[1]/born << "\n";
  generator()->log() << "c0="<< virt[0]/born << "\n";

  const cPDVector& proc = mePartonData();
  double CA1 = 3.;
  double CF1 = 4./3.;
  double msq =  me2();
  double two(2.);
  double cTaTb = colourCorrelatedME2(pair<int,int>(0,1))*
    (proc[0]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/abs((meMomenta()[0] + meMomenta()[1]).m2()));
  double cT1T2 = colourCorrelatedME2(pair<int,int>(2,3))*
    (proc[2]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/abs((meMomenta()[2] + meMomenta()[3]).m2()));
  double cT1T3 = colourCorrelatedME2(pair<int,int>(2,4))*
   (proc[2]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[2] + meMomenta()[4]).m2())));
  double cT2T3 = colourCorrelatedME2(pair<int,int>(3,4))*
  (proc[3]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[3] + meMomenta()[4]).m2())));
  double cTaT1 = colourCorrelatedME2(pair<int,int>(0,2))*
  (proc[0]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[0] + meMomenta()[2]).m2())));
  double cTaT2 = colourCorrelatedME2(pair<int,int>(0,3))*
  (proc[0]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[0] + meMomenta()[3]).m2())));
  double cTaT3 = colourCorrelatedME2(pair<int,int>(0,4))*
						     (proc[0]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[0] + meMomenta()[4]).m2())));
  double cTbT1 = colourCorrelatedME2(pair<int,int>(1,2))*
						     (proc[1]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[1] + meMomenta()[2]).m2())));
  double cTbT2 = colourCorrelatedME2(pair<int,int>(1,3))*
						     (proc[1]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[1] + meMomenta()[3]).m2())));
  double cTbT3 = colourCorrelatedME2(pair<int,int>(1,4))*
  (proc[1]->id() == ParticleID::g ? CA1 : CF1)*log(scale*GeV2/(abs((meMomenta()[1] + meMomenta()[4]).m2())));
 
  generator()->log() <<"CA1"<<CA1<< "\n";
  
  double gammaQuark = (3./2.)*CF1;
  double gammaGluon = (11./6.)*CA1 - (1./3.)*NF;

  double cfac(4.*gammaQuark+gammaGluon);
    
  double res1 =-cfac*msq+2.*(cT1T2 + cT1T3 + cTaT1 + cTbT1 + cT2T3 + cTaT2+ cTbT2 + cTaT3 + cTbT3 + cTaTb);


  double gs2(SM().alphaS()*4.*Constants::pi);
  born*=pow(lastSHat()/GeV2,2)*gs2;
  generator()->log()<<"res1/born="<<res1/me2()<<"\n";
  */

  /*
  double mm = me2();
  born *= pow(lastSHat()/GeV2,2);
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

  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void Amplitudehqqbarkkbarg::persistentOutput(PersistentOStream & os) const {
  os << theGaugeThreshold;
}

void Amplitudehqqbarkkbarg::persistentInput(PersistentIStream & is, int) {
  is >> theGaugeThreshold;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<Amplitudehqqbarkkbarg,AmplitudeBase>
  describeHJetsAmplitudehqqbarkkbarg("HJets::Amplitudehqqbarkkbarg", "HwMatchboxBuiltin.so HJets.so");

void Amplitudehqqbarkkbarg::Init() {

  static ClassDocumentation<Amplitudehqqbarkkbarg> documentation
    ("Helicity amplitudes for 0 -> h q qbar k kbar g");


  static Parameter<Amplitudehqqbarkkbarg,double> interfaceGaugeThreshold
    ("GaugeThreshold",
     "The threshold for gauge check failure.",
     &Amplitudehqqbarkkbarg::theGaugeThreshold, 0.1, 0.0, 0,
     false, false, Interface::lowerlim);

  static Switch<Amplitudehqqbarkkbarg,bool> interfaceAccuracyCheck
    ("AccuracyCheck",
     "Turn on check of accuracy of box,pentagon, and hexagons.",
     &Amplitudehqqbarkkbarg::theAccuracyCheck, false, false, false);
  static SwitchOption interfaceAccuracyCheckOn
    (interfaceAccuracyCheck,
     "On",
     "On",
     true);
  static SwitchOption interfaceAccuracyCheckOff
    (interfaceAccuracyCheck,
     "Off",
     "Off",
     false);


}

