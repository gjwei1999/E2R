#ifndef EToRConvForProton_h
#define EToRConvForProton_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "EnergyToRangeConverter.hh"


class EToRConvForProton : public EnergyToRangeConverter
{
  public: 
  //  constructor
  EToRConvForProton();

  public:
  //  destructor
  virtual ~EToRConvForProton();

  virtual G4double Convert(G4double rangeCut, const G4Material* material);

  // reset Loss Table and Range Vectors
  virtual void Reset();

  protected:
    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy
                                ) ;

  protected:
    G4double Mass;
    G4double Z;  
    G4double tau0;
    G4double taul;
    G4double taum;
    G4double ionpot;
    G4double ca;
    G4double cba;
    G4double cc;  
};


#endif









