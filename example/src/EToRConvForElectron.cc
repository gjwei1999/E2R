#include "EToRConvForElectron.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4PhysicsLogVector.hh"

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

EToRConvForElectron::EToRConvForElectron() 
  : EnergyToRangeConverter(),
    Mass(0.0),
    Z(-1.),  
    taul(0.0),
    ionpot(0.0),
    ionpotlog(-1.0e-10),
    bremfactor(0.1)
{    
  theParticle =  G4ParticleTable::GetParticleTable()->FindParticle("e-");
  if (theParticle ==0) {

      G4cout << " EToRConvForElectron::EToRConvForElectron() ";
      G4cout << " Electron is not defined !!" << G4endl;

  } else {
    Mass = theParticle->GetPDGMass();
  }
}

EToRConvForElectron::~EToRConvForElectron()
{ 
}


// **********************************************************************
// ************************* ComputeLoss ********************************
// **********************************************************************
G4double EToRConvForElectron::ComputeLoss(G4double AtomicNumber,
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

  if(tau<taul) {
    G4double t1 = taul+1.;
    G4double t2 = taul+2.;
    G4double tsq = taul*taul;
    G4double beta2 = taul*t2/(t1*t1);
    G4double f = 1.-beta2+std::log(tsq/2.)
                  +(0.5+0.25*tsq+(1.+2.*taul)*std::log(0.5))/(t1*t1);
    dEdx = (std::log(2.*taul+4.)-2.*ionpotlog+f)/beta2;
    dEdx = twopi_mc2_rcl2*Z*dEdx;
    G4double clow = dEdx*std::sqrt(taul);
    dEdx = clow/std::sqrt(KineticEnergy/Mass);

  } else {
    G4double t1 = tau+1.;
    G4double t2 = tau+2.;
    G4double tsq = tau*tau;
    G4double beta2 = tau*t2/(t1*t1);
    G4double f = 1.-beta2+std::log(tsq/2.)
                   +(0.5+0.25*tsq+(1.+2.*tau)*std::log(0.5))/(t1*t1);
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
