#ifndef EnergyToRangeConverter_h
#define EnergyToRangeConverter_h 1

#include "globals.hh"
#include <cmath>
#include "G4ios.hh"
#include <vector>

#include "G4ParticleDefinition.hh"

#include "G4PhysicsTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
class G4PhysicsLogVector;

class EnergyToRangeConverter
{
  public: // with description
  //  constructor
  EnergyToRangeConverter();

  //  copy constructor
//  EnergyToRangeConverter(const EnergyToRangeConverter &right);

//  EnergyToRangeConverter & operator=(const EnergyToRangeConverter &right);

  public:
  //  destructor
  virtual ~EnergyToRangeConverter();

  // equal opperators
//  G4int operator==(const EnergyToRangeConverter &right) const;
//  G4int operator!=(const EnergyToRangeConverter &right) const;

  public: // with description 
  // calculate energy cut from given range cut for the material
  virtual G4double Convert(G4double rangeCut, const G4Material* material);
      
      
  //  set energy range for all particle type
  static void SetEnergyRange(G4double lowedge, G4double highedge);

  //  get energy range for all particle type
  static G4double GetLowEdgeEnergy();
  static G4double GetHighEdgeEnergy();

  //  get/set max cut energy for all particle type
  static G4double GetMaxEnergyCut();
  static void SetMaxEnergyCut(G4double value);
  
  // return pointer to the particle type which this converter takes care
  const G4ParticleDefinition* GetParticleType() const;

  // return the Loss Table
  const  G4PhysicsTable* GetLossTable() const;   
   //-------------- Loss Table ------------------------------------------
   // theLossTable is a collection of loss vectors for all elements.
   // Each loss vector has energy loss values (cross section values
   // for neutral particles) which are calculated by
   // ComputeLoss(G4double AtomicNumber,G4double KineticEnergy).
   // ComputeLoss method is pure virtual and should be provided for each 
   // particle type

  // reset Loss Table and Range Vectors
  virtual void Reset();
    
 protected:

    static G4double               LowestEnergy, HighestEnergy;
    static G4double               MaxEnergyCut; 
//    G4double                      fMaxEnergyCut;
   
    const G4ParticleDefinition*   theParticle;
    typedef G4PhysicsTable        G4LossTable;
    G4LossTable*                  theLossTable;
    G4int                         NumberOfElements;
  
    typedef G4PhysicsLogVector    G4LossVector;
    const G4int                   TotBin;

  protected:// with description  
    virtual void BuildLossTable();

    virtual G4double ComputeLoss(G4double AtomicNumber,
                                 G4double KineticEnergy
			       ) =0;

  //-------------- Range Table ------------------------------------------
  protected:
    typedef G4PhysicsLogVector G4RangeVector;

    virtual void BuildRangeVector(const G4Material* aMaterial,
                                   G4RangeVector* rangeVector);

//    std::vector< G4RangeVector* > fRangeVectorStore;   
      
  protected:
    G4double ConvertEnergyToCut(
                                        G4RangeVector* theRangeVector,
                                        G4double       kineticenergy 
                                      ) const;

};

inline 
 const G4ParticleDefinition* EnergyToRangeConverter::GetParticleType() const
{
   return theParticle;
}
#endif








