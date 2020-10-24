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
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


#include "EToRConvert.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;









G4cout <<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;

	gStyle->SetPalette(53);
	

	TString path = "/home/jiaweiguo/work/test/build/"; //plot path should be modified
	TString filename = path + "test.root";

	TFile * fs = new TFile(filename.Data(), "RECREATE");

	TH1D * h0 = new TH1D("electron_Pb", "", 1000.0, 0.0, 100.0);
	TH1D * h1 = new TH1D("gamma_Pb", "", 1000.0, 0.0, 100.0);
	TH1D * h2 = new TH1D("positron_Pb", "", 1000.0, 0.0, 100.0);
	TH1D * h3 = new TH1D("proton_Pb", "", 1000.0, 0.0, 100.0);

	h0->SetDirectory(fs);
	h1->SetDirectory(fs);
	h2->SetDirectory(fs);
	h3->SetDirectory(fs);
	
	auto nistManager = G4NistManager::Instance();
        G4String materialName="G4_Pb";			//define the material
        auto material = nistManager->FindOrBuildMaterial(materialName);
 
	G4double range;

	//define the particles
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  	G4String particleName0 = "e-" ;
  	G4ParticleDefinition* ptcl0 = particleTable->FindParticle(particleName0);

	G4String particleName1 = "gamma" ;
        G4ParticleDefinition* ptcl1 = particleTable->FindParticle(particleName1);

	G4String particleName2 = "e+" ;
        G4ParticleDefinition* ptcl2 = particleTable->FindParticle(particleName2);

	G4String particleName3 = "proton" ;
        G4ParticleDefinition* ptcl3 = particleTable->FindParticle(particleName3);

	G4double egy;

	int i=0;
	for(i=0;i<1000;i++){
		egy = 0.1*i;
		range = EToRConvert(egy,material,ptcl0);	//the function
		G4cout<<"The electron in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl; 
		h0->Fill(egy,range);
	}

	for(i=0;i<1000;i++){
                egy = 0.1*i;
                range = EToRConvert(egy,material,ptcl1);	//the function
                G4cout<<"The gamma in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
                h1->Fill(egy,range);
        }

	for(i=0;i<1000;i++){
                egy = 0.1*i;
                range = EToRConvert(egy,material,ptcl2);	//the function
                G4cout<<"The positron in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
                h2->Fill(egy,range);
        }


	for(i=0;i<1000;i++){
                egy = 0.1*i;
                range = EToRConvert(egy,material,ptcl3);	//the function
                G4cout<<"The proton in material "<< material->GetName() <<" can move "<<range<<" mm."<<G4endl;
                h3->Fill(egy,range);
        }

	for(i=0;i<1000;i++){
		h0->SetBinError(i,0.0);
		h1->SetBinError(i,0.0);
		h2->SetBinError(i,0.0);
		h3->SetBinError(i,0.0);
	}

	h0->GetXaxis()->SetTitle("E(MeV)");
	h0->GetYaxis()->SetTitle("Range(mm)");
	h0->SetTitle("electron_Pb");

	h1->GetXaxis()->SetTitle("E(MeV)");
        h1->GetYaxis()->SetTitle("Range(mm)");
        h1->SetTitle("gamma_Pb");

	h2->GetXaxis()->SetTitle("E(MeV)");
        h2->GetYaxis()->SetTitle("Range(mm)");
        h2->SetTitle("positron_Pb");

	h3->GetXaxis()->SetTitle("E(MeV)");
        h3->GetYaxis()->SetTitle("Range(mm)");
        h3->SetTitle("proton_Pb");
	
	
	fs->Write();


G4cout <<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<G4endl;















}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

