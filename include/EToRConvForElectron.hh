//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
#ifndef EToRConvForElectron_h
#define EToRConvForElectron_h 1

#include "globals.hh"
#include "G4ios.hh"
#include <vector>

#include "EnergyToRangeConverter.hh"


class EToRConvForElectron : public EnergyToRangeConverter
{
  public: 
  //  constructor
  EToRConvForElectron();

  public:
  //  destructor
  virtual ~EToRConvForElectron();

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









