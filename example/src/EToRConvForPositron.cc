#include "EToRConvForPositron.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

EToRConvForPositron::EToRConvForPositron() 
  : EnergyToRangeConverter(),
    Mass(0.0),
    Z(-1.),  
    taul(0.0),
    ionpot(0.0),
    ionpotlog(-1.0e-10),
    bremfactor(0.1)
{    
  theParticle =  G4ParticleTable::GetParticleTable()->FindParticle("e+");
  if (theParticle ==0) {

      G4cout << " EToRConvForPositron::EToRConvForPositron() ";
      G4cout << " Positron is not defined !!" << G4endl;

  } else {
    Mass = theParticle->GetPDGMass();
  }
}

EToRConvForPositron::~EToRConvForPositron()
{ 
}




// **********************************************************************
// ************************* ComputeLoss ********************************
// **********************************************************************
G4double EToRConvForPositron::ComputeLoss(G4double AtomicNumber,
					    G4double KineticEnergy) 
{
  const  G4double cbr1=0.02, cbr2=-5.7e-5, cbr3=1., cbr4=0.072;
  const  G4double Tlow=10.*keV, Thigh=1.*GeV;

  //  calculate dE/dx for electrons
  if( std::fabs(AtomicNumber-Z)>0.1 ) {
    Z = AtomicNumber;
    taul = Tlow/Mass;
    ionpot = 1.6e-5*MeV*std::exp(0.9*std::log(Z))/Mass;
    ionpotlog = std::log(ionpot);
  } 

  G4double tau = KineticEnergy/Mass;
  G4double dEdx;

 if(tau<taul)
  {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double     f = 2.*std::log(taul)
                     -(6.*taul+1.5*tsq-taul*(1.-tsq/3.)/t2-tsq*(0.5-tsq/12.)/
                       (t2*t2))/(t1*t1);
    dEdx = (std::log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;
    G4double clow = dEdx*std::sqrt(taul);
    dEdx = clow/std::sqrt(KineticEnergy/Mass);

  } else {
   G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 2.*std::log(tau)
                 - (6.*tau+1.5*tsq-tau*(1.-tsq/3.)/t2-tsq*(0.5-tsq/12.)/
                     (t2*t2))/(t1*t1);
    dEdx = (std::log(2.*tau+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;

    // loss from bremsstrahlung follows
    G4double cbrem = (cbr1+cbr2*Z)
                       *(cbr3+cbr4*std::log(KineticEnergy/Thigh));
    cbrem = Z*(Z+1.)*cbrem*tau/beta2;
    cbrem *= bremfactor ;
    dEdx += twopi_mc2_rcl2*cbrem;
  }
  return dEdx;
}


