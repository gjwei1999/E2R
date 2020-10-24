#ifndef EToRConvForPositron_h
#define EToRConvForPositron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "EnergyToRangeConverter.hh"


class EToRConvForPositron : public EnergyToRangeConverter
{
  public:
  //  constructor
  EToRConvForPositron();

  public:
  //  destructor
  virtual ~EToRConvForPositron();

  protected:
    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy
                                ) ;

  protected:
    G4double Mass;
    G4double Z;  
    G4double taul;
    G4double ionpot;
    G4double ionpotlog;
    G4double bremfactor;
};


#endif









