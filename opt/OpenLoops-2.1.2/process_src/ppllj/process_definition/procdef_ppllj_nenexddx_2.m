
FeynArtsProcess = {F[1, {1}], -F[1, {1}]} -> {-F[4, {1}], F[4, {1}]};

SortExternal = True;

OpenLoopsModel = "SM";

CreateTopologiesOptions = {
  ExcludeTopologies -> {Snails, WFCorrectionCTs, TadpoleCTs, Loops[5]},
  Adjacencies -> {3, 4}
};

InsertFieldsOptions = {
  Model -> {"SMQCD", "SMQCDR2"},
  GenericModel -> "Lorentz",
  InsertionLevel -> {Particles},
  Restrictions -> {ExcludeParticles -> {S[2 | 3], V[4]}, NoQuarkMixing}
};

UnitaryGauge = True;

ColourCorrelations = Automatic;

OTFColourCorrelations = Automatic;

SpinCorrelatedHardFactor = Automatic;

SubProcessName = Automatic;

SelectCoupling = MemberQ[{2}, Exponent[#1, eQED]] & ;

SelectInterference = {
  eQED -> {4}
};

SelectTreeDiagrams = True & ;

SelectLoopDiagrams = True & ;

SelectCTDiagrams = True & ;

ReplaceOSw = False;

SetParameters = {
  CKMORDER -> 0,
  nc -> 3,
  nf -> 6,
  nfl -> 3,
  MU -> 0,
  MD -> 0,
  MS -> 0,
  LeadingColour -> 0,
  POLSEL -> 1
};

ChannelMap = {
  {"nenexbbx", "MB=0"},
  {"nlnlxddx"},
  {"nlnlxbbx"},
  {"nenexssx"}
};

Approximation = "";

QED = 0;

ForceLoops = Automatic;

ForceLoopsInclude = Automatic;

NonZeroHels = Null;

OnTheFlyMode = Automatic;

noQCD = False;

noEW = False;

SetExternalSubtrees = {};
