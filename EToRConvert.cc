#include "EToRConvForElectron.hh"
#include "EToRConvForGamma.hh"
#include "EToRConvForPositron.hh"
#include "EToRConvForProton.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"


G4double EToRConvert(G4double KineticEnergy, const G4Material* material, const G4ParticleDefinition* ptcl);


G4double EToRConvert(G4double KineticEnergy, const G4Material* material, const G4ParticleDefinition* ptcl){

	G4double range;

	if(ptcl->GetParticleName()=="e-"){
		EnergyToRangeConverter* converter = new EToRConvForElectron();
        	range = converter->Convert(KineticEnergy, material);
		G4cout<<"The electron in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
		return range;
	}
	else 
	if(ptcl->GetParticleName()=="gamma"){
		EnergyToRangeConverter* converter = new EToRConvForGamma();
		range = converter->Convert(KineticEnergy, material);
		G4cout<<"The gamma in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
		return range;
	}
	else if(ptcl->GetParticleName()=="e+"){
		EnergyToRangeConverter* converter = new EToRConvForPositron();
		range = converter->Convert(KineticEnergy, material);
		G4cout<<"The positron in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
		return range;
	}
	else if(ptcl->GetParticleName()=="proton"){
		EnergyToRangeConverter* converter = new EToRConvForProton();
		range = converter->Convert(KineticEnergy, material);
		G4cout<<"The proton in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
		return range;
	}
	else{
		G4cout << " The particle can only be e-, e+, gamma, or proton." <<G4endl;
		return 0;
	}
	
}

