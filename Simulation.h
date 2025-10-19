// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <omp.h>
#include <algorithm>
#include <set>
#include <memory>
#include "SimulationBox.h"
#include "Grid.h"
#include "ParticleInfo.h"
#include "MiscStructs.h"
#include "Filament.h"
#include "SpringBondInfo.h"

using namespace std;

class Simulation
{
	private:
		const double pi = 3.14159265359;
        const double boltzmannConstant = 1.38e-5; // [pN-um/K]
        
        vector<double> unfolding_rate_constants = {2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5, 2.5e-5};
        vector<double> unfolding_alpha = {2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114, 2.5/4.114};
        
        
        double temperature;
        double dt;
        double zeta;
        double zetaNucleus;
        int numSteps;
        int numThreads;
        int snapshotStep;
        int updateGridStep;
        double kappa;
        double kappa_Nuc;
        double persistenceLength;
        double L_zero;
        double eta;
        double thermal_force;
        unsigned seed;
        double totalSimulationTime;
        double snapshotTime;
        double currSimTime;
        double K_Area_nuc;
        double K_Area_CM;
        double AreaNucleusZero;
        double AreaCMZero;
        double F_Mag_Front;
        double Ratio;
        double F_Mag_back;
        double areaForceNucleus;
        double areaForceCortex;
        double OMEGA;
        double d_seal;
        int numNucBeads;
        int numCMBeads;
        int numLeadingEdgeBeads;
        double EquilibriumAreaUpper = 0.0;
        double EquilibriumAreaLower = 0.0;
        std::vector<int> CM_List;
        std::vector<int> Nuc_List;
        std::vector<int> upperCompartment;
        std::vector<int> lowerCompartment;
        int CBM_1 = -1;
        int CBN_1 = -1;
        int CBN_2 = -1;
        int CBM_2 = -1;
        double timeSRS = 0.0;
        double timeLag = 123.36953; 

        double omega;
        double K_Stretch;
        double K_Stretch_nuc;
        double kCenter;
        double kExcl_CM;
        double kExcl_Obs;

        double d_Excl;
        double thresholdDistCB;

        bool thermalForcesOn;

        bool excludedVolOn;
        bool CB_On;
		bool blebbing;
		bool changeDirection;
        double beadRadius;
        
        std::vector <int> removalList;

        std::vector <std::vector <int>> excludedVolumeList;

        std::string outputFileName;

        std::string positionFileName;
        std::string bondFileName;
        std::string angleFileName;

        SimulationBox simBox;
        Grid particleGrid;

        ParticleInfo pinfo;

        std::vector <Filament> filaments;
        std::vector <std::unique_ptr<BondInfo>> bondTypes;

        TaggedVector f_TaggedVector;

        std::uniform_real_distribution<> dis;
        std::normal_distribution<> gaussianDis;

        std::vector <std::mt19937> generators;

        std::vector<std::vector <std::vector <Coordinate>>> tempForces;

        // 0 - bond, 1 - angle, 2 - excluded volume cell beads, 3 - area force, 4 - nuclear centering force, 5 - protrusive force, 6 - retractive force, 7 - area force plasma membrane 1, 8 - area force pm 2, 9 - excluded volume obs beads
        int num_force_types = 10;
        std::vector<std::vector <Coordinate>> tempStresses;
        
        Coordinate r_nuc_com = Coordinate (0.0, 0.0, 0.0);

	public:
		Simulation();

		void readParameterFile(std::string file_name);
        void run();
        void writeXYZFile(std::vector<std::string>& xyzArray);
        void writeBondFile(std::vector<std::string>& bondArray);
        void writeEndingFile(std::ofstream& outputEND);
        double calcForces();
        void putAllParticlesInGrid();
        void moveParticles();
        void readAngleFile();
        void readBondFile();
        void readXYZ();
        void initDerivedValues();
		void polaritySwitch();
        void shiftBeadTypes(int shiftAmount);					  

        double bondForceMag();

        double AreaConservationDerivativeXb(Coordinate a, Coordinate b, Coordinate c);
        double AreaConservationDerivativeYb(Coordinate a, Coordinate b, Coordinate c);

        double MeasureAreaNucleus();
        double MeasureAreaCM();
        double MeasureAreaTriangle(Coordinate a, Coordinate b);
        std::vector<int> compartmentalizeUpperSys(int a, int b, int c, int d);
        std::vector<int> compartmentalizeLowerSys(int a, int b, int c, int d);
        double MeasureArea(std::vector<int> perimeter_list);
        double conserveAreaOfCompartments(std::vector<int> perimeter_list, double EquilibriumArea, double areaConservationK, int forceType);
        double conserveAreaOfCompartmentsRing(std::vector<int> perimeter_list, std::vector<int> perimeter_list_two, double EquilibriumArea, double areaConservationK, int forceType);
        double calcDeterminant(Coordinate firstRow, Coordinate secondRow, Coordinate thirdRow);
        Coordinate closestDistanceBetweenLines(Coordinate a0, Coordinate a1, Coordinate b0, Coordinate b1, Coordinate &aClose, Coordinate &bClose);
        
        double HeavisideTheta(double val);
        double DiracDelta(double val);
};


#endif
