// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#include "Simulation.h"

Simulation::Simulation()
{
    //Default parameter values:
    temperature = 300.0;
    eta = 21281.24398;
    dt = 3.08E-3;
    totalSimulationTime = 10000.0;
    snapshotTime = 30.84;
    persistenceLength = 369.5652174; 
    numThreads = 1;
    thermalForcesOn = false;
    K_Stretch = 1000000.0;
    d_Excl = 0.2608696; 
    excludedVolOn = true;
    K_Area_CM = 6000.0; 
    K_Area_nuc = 300000.0; 
    AreaNucleusZero = 0.0;
    F_Mag_Front = 50.0;
    Ratio = 0.1176;
    d_seal = 1.3043478; 
    numNucBeads = 100;
    numCMBeads = 200; 
    numLeadingEdgeBeads = 41; 
    for (int i = numNucBeads; i < (numNucBeads + numCMBeads); i++) 
    {
        CM_List.push_back(i);

    }
        std::cout << endl; //Defined CM_List - conserve area in it.
    for (int i = 0; i < numNucBeads; i++)
    {
        Nuc_List.push_back(i);

    }
    
    beadRadius = 0.0434783; 
    
    
    currSimTime = 0.0;
    
    outputFileName = "";
    
    initDerivedValues();

    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    dis = std::uniform_real_distribution<> (0.0, 1.0);
    gaussianDis = std::normal_distribution<> (0.0, 1.0);
	for (int i = 0; i < omp_get_max_threads(); i++)
	{
		generators.emplace_back(std::mt19937(gen));
	}

}

void Simulation::initDerivedValues()
{
    zeta = 1.0*(6*pi*eta*beadRadius);
    zetaNucleus = 0.75*(6*pi*eta*beadRadius);
    
    snapshotStep = (int)round(snapshotTime / dt);
    numSteps = (int)round(totalSimulationTime / dt);
    kappa = temperature * boltzmannConstant * persistenceLength;
    kappa_Nuc = 1.0*(temperature * boltzmannConstant * persistenceLength);
	
    thermal_force = sqrt(2*boltzmannConstant*temperature*zeta/dt);
    F_Mag_back = F_Mag_Front/Ratio;
}

void Simulation::readParameterFile(std::string fileName)
{

    ifstream inputfile (fileName);
    double timeBetweenDeltaSteps = 0.0;
    
    double boxLengths[] = {600.0, 600.0, 5.0};
    
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        std::string line;
        std::string delimiter = " ";
              
		while(std::getline(inputfile, line)) 
		{	

            size_t position = line.find(delimiter);
            std::string category = line.substr(0, position);
            
            line.erase(0, position + delimiter.length());
            
            size_t positionEndLine;
            std::string value;
            if(((positionEndLine = line.find("\r")) != std::string::npos))
            {
                
                value = line.substr(0, line.size()-1);
            }
            else
            {
                value = line.substr(0, line.size());
            }

            //Parameter values in input.txt:
            if(category.compare("dt") == 0)
            {
                dt = std::stod(value);
            }
            else if(category.compare("totalSimulationTime") == 0)
            {
                totalSimulationTime = std::stod(value);
            }
            else if(category.compare("snapshotTime") == 0)
            {
                snapshotTime = std::stod(value);
            }
            else if(category.compare("persistenceLength") == 0)
            {
                persistenceLength = std::stod(value);
            }
            else if(category.compare("temperature") == 0)
            {
                temperature = std::stod(value);
            }
            else if(category.compare("eta") == 0)
            {
                eta = std::stod(value);
            }
            else if(category.compare("positions") == 0)
            {
                positionFileName = value;
            }
            else if(category.compare("bonds") == 0)
            {
                bondFileName = value;
            }
            else if(category.compare("angles") == 0)
            {
                angleFileName = value;
            }
            else if(category.compare("numThreads") == 0)
            {
                numThreads = std::stoi(value);
            }
            else if(category.compare("thermalForcesOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> thermalForcesOn;
            }
            else if(category.compare("K_Stretch") == 0)
            {
                K_Stretch = std::stod(value);
            }
            else if(category.compare("K_Stretch_nuc") == 0)
            {
                K_Stretch_nuc = std::stod(value);
            }
            else if(category.compare("kCenter") == 0)
            {
                kCenter = std::stod(value);
            }
            else if(category.compare("kExcl_CM") == 0)
            {
                kExcl_CM = std::stod(value);
            }
            else if(category.compare("kExcl_Obs") == 0)
            {
                kExcl_Obs = std::stod(value);
            }
            else if(category.compare("d_Excl") == 0)
            {
                d_Excl = std::stod(value);
            }
            else if(category.compare("thresholdDistCB") == 0)
            {
                thresholdDistCB = std::stod(value);
            }
            else if(category.compare("outputName") == 0)
            {
                outputFileName = value;
            }
            else if(category.compare("excludedVolOn") == 0)
            {
                istringstream(value) >> std::boolalpha >> excludedVolOn;
            }
            else if(category.compare("CB_On") == 0)
            {
                istringstream(value) >> std::boolalpha >> CB_On;
            }
			else if(category.compare("blebbing") == 0)
			{
				istringstream(value) >> std::boolalpha >> blebbing;
			}
			else if(category.compare("changeDirection") == 0)
			{
				istringstream(value) >> std::boolalpha >> changeDirection;
			}
            else if(category.compare("beadRadius") == 0)
            {
                beadRadius = std::stod(value);
            }
            else if(category.compare("L_zero") == 0)
            {
                L_zero = std::stod(value);
            }
            else if(category.compare("K_Area_nuc")==0)
            {
                K_Area_nuc = std::stod(value);
            }
            else if(category.compare("K_Area_CM")==0)
            {
                K_Area_CM = std::stod(value);
            }
            else if(category.compare("F_Mag_Front")==0)
            {
                F_Mag_Front = std::stod(value);
            }
            else if(category.compare("Ratio")==0)
            {
                Ratio = std::stod(value);
            }
            else if(category.compare("OMEGA")==0)
            {
                OMEGA = std::stod(value);
            }
            else if(category.compare("boxDimensions") == 0)
            {
            
                std::string tmp_str = value;
                std::string delimiter = " ";
                
                int count = 0;
                size_t pos = 0;
                std::string token;
                while ((pos = tmp_str.find(delimiter)) != std::string::npos)
                {
                    token = tmp_str.substr(0, pos);
                    tmp_str.erase(0, pos + delimiter.length());
                    
                    if(count >= 0 && count <= 3)
                    {
                        boxLengths[count] = std::stod(token);
                    }
                    count++;
                }
                if(count == 2)
                {
                    boxLengths[count] = std::stod(tmp_str);
                }
            }
            else
            {
                throw std::runtime_error("Unknown value in input file associated with: " + category);
            }
            
            linecount++;
			
		}
        
		inputfile.close();
    }
    else
    {

        throw std::runtime_error("main could not open: " + fileName);
    }
    
    initDerivedValues();
    
    Coordinate boxSize = {boxLengths[0], boxLengths[1], boxLengths[2]};
    bool periodicX = false; 
    bool periodicY = false;
    bool periodicZ = false;
    
    simBox = SimulationBox(boxSize, periodicX, periodicY, periodicZ);
    particleGrid = Grid(boxSize, periodicX, periodicY, periodicZ);
    particleGrid.setBinSize(L_zero*3.0);
    
    if(!positionFileName.empty())
    {
        readXYZ();
    }
    
    if(bondFileName.empty() && !positionFileName.empty())
    {
   
        bondFileName = positionFileName.substr(0, positionFileName.length() - 4) + ".bnd";
        
        cout << "Assuming bond file is named: " << bondFileName << endl;
    }
    readBondFile();
    
    if(angleFileName.empty() && !positionFileName.empty())
    {
        angleFileName = positionFileName.substr(0, positionFileName.length() - 4) + ".ang";
        
        cout << "Assuming angle file is named: " << angleFileName << endl;
    }
    readAngleFile();
    double updateGridTime = 1.0; 
    updateGridStep = (int)round(updateGridTime / dt);

    bondTypes.push_back(std::unique_ptr<BondInfo> (new SpringBondInfo(K_Stretch_nuc, L_zero)));
    
    bondTypes.push_back(std::unique_ptr<BondInfo> (new SpringBondInfo(K_Stretch, L_zero)));

    bondTypes.push_back(std::unique_ptr<BondInfo> (new SpringBondInfo(kCenter, 0.0))); 

    AreaNucleusZero = MeasureAreaNucleus();
    AreaCMZero = MeasureAreaCM();
    std::cout << AreaCMZero << " " << AreaNucleusZero << endl;
    int numExtraBeads = 5;
    int numExtraBeadsTwo = 15;
    int numMinusNuc = 0; 
    int plusPseudopod = 0;

    for (int p = 0; p < pinfo.getNumParticles(); p++)
    {
        if(p >= numCMBeads - numExtraBeads - numExtraBeadsTwo - numMinusNuc + plusPseudopod && p < numCMBeads - numExtraBeads - numMinusNuc + plusPseudopod){
            pinfo.setBeadToFilamentByIndex(p, 2);
        }
        else if(p >= numCMBeads - numExtraBeads - numMinusNuc + plusPseudopod && p <= numCMBeads + numExtraBeads - numMinusNuc + plusPseudopod){
            pinfo.setBeadToFilamentByIndex(p, 3); 
        }
        else if(p > numCMBeads + numExtraBeads - numMinusNuc + plusPseudopod && p <= numCMBeads + numExtraBeads + numExtraBeadsTwo - numMinusNuc + plusPseudopod){
            pinfo.setBeadToFilamentByIndex(p, 4);
        }

        else if(p < numNucBeads+numCMBeads && p > numNucBeads - 1){
            pinfo.setBeadToFilamentByIndex(p, 5);
        }
    }

}

void Simulation::readXYZ()
{
    ifstream inputfile (positionFileName);
    
    int numParticles = 0;
    
    if (inputfile.is_open())
	{
		int linecount = 0;
        
        int currFilamentType = -1;
        int newFilamentTag = -1;
        
        std::string line;
        
		while(std::getline(inputfile, line))
		{
            std::istringstream iss(line);
            
            if(linecount == 0)
            {
                iss >> numParticles;
            }
            else if(linecount == 1)
            {
                std::string tempString;
                iss >> tempString;
                
                int pos = tempString.find("=");
                
                std::string timeString = tempString.substr(pos+1,tempString.length());
                
                currSimTime = std::stod(timeString);

            }
            else if(linecount > 1)
            {
                int tag, type;
                double x,y,z;
                int viscosityType;
                
                iss >> type >> tag >> x >> y >> z >> viscosityType;
                
                if(type != currFilamentType)
				{
					Filament newFil;
					filaments.push_back(newFil);
                    newFilamentTag = f_TaggedVector.add();
					
                    currFilamentType = type;
				}

				Coordinate newCoordinate = {x, y, z};
				
				Particle newParticle(x, y, z);
                newParticle.setType(type); 
                
                int numStressMeasurements = 2;
                
				int newTag = pinfo.addParticle(newParticle, numStressMeasurements);
                if(newTag != tag)
                {
                    pinfo.setTagAtPos(pinfo.getNumParticles()-1, tag);
                    newTag = tag;
                }
                
                pinfo.setViscosityTypeByIndex(pinfo.getNumParticles()-1, viscosityType);
                
				int priorTag = filaments[newFilamentTag].addParticleBack(newTag);
            }
            
            linecount++;
        }
        
        if(linecount-2 != numParticles)
        {
            std::cout << "Error. " << numParticles << " particles desired in simulation but " << linecount-2 << " particles in " << this->positionFileName;
        }
        
        inputfile.close();
    }
    else
    {
        cout << "Could not open " + positionFileName << endl;
    }
}

void Simulation::readBondFile()
{
    ifstream inputfile (bondFileName);
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        int numBonds = 0;
        int countBonds = 0;
        
        int numParticles = pinfo.getNumParticles();
        
        std::string line;
        
        int maxCurrTag = pinfo.getCurrMaxTag();
        
		while(std::getline(inputfile, line))
		{
            std::istringstream iss(line);
            
            if(linecount == 0)
            {
                iss >> numBonds;
            }
            else if(linecount > 1)
            {
                int particleI, particleJ, bondType;
                double bond_creationTime, bond_destructionTime;
                
                iss >> bondType >> particleI >> particleJ >> bond_creationTime >> bond_destructionTime;
                
                if(particleI < 0 || particleI > maxCurrTag || particleJ < 0 || particleJ > maxCurrTag)
                {
                    throw std::runtime_error(std::to_string(particleI) + " or " + std::to_string(particleJ) + " tags are larger than current max tag " + std::to_string(maxCurrTag));
                }

                pinfo.addBond(Bond {particleI, particleJ, bondType, bond_creationTime, bond_destructionTime});
                countBonds++;
            }
            
            linecount++;
        }

        if(countBonds != numBonds)
        {
            std::cout << "Error. " << numBonds << " bonds desired in simulation but " << countBonds << " bonds in " << this->bondFileName;
        }
        
        inputfile.close();
    }
}

void Simulation::readAngleFile()
{
    ifstream inputfile (angleFileName);
    if (inputfile.is_open())
	{
		int linecount = 0;
		
        int numAngles = 0;
        int countAngles = 0;
        
        std::string line;
        
		while(std::getline(inputfile, line))
		{
            std::istringstream iss(line);
            
            if(linecount == 0)
            {
                iss >> numAngles;
            }
            else if(linecount > 1)
            {
                int particleI, particleJ, particleK, angleType;
                iss >> angleType >> particleI >> particleJ >> particleK;

                pinfo.addAngle(Angle {particleI, particleJ, particleK, angleType});
                countAngles++;

            }
            
            linecount++;
        }
        
        if(countAngles != numAngles)
        {
            throw std::runtime_error("Number of angles desired in simulation is " + std::to_string(numAngles) + " but " + std::to_string(countAngles) + " angles in " + angleFileName);
        }
        
        inputfile.close();
    }
}

void Simulation::moveParticles()
{
    int numParticles = pinfo.getNumParticles();
    
    #pragma omp parallel for
    for(int i = 0; i < numParticles; i++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(i);
        Coordinate currForce = pinfo.getForceByIndex(i);
        
        int filamentTag = pinfo.getBeadToFilamentByIndex(i);
        
        Coordinate nextPosition = {0.0, 0.0, 0.0};

        if(filamentTag != 1)
        {
            double vx, vy, vz;
            if(filamentTag == 0)
            {
                vx = currForce.x / (zetaNucleus); 
                vy = currForce.y / (zetaNucleus);
                vz = currForce.z / (zetaNucleus);
            }
            else
            {
                vx = currForce.x / zeta;
                vy = currForce.y / zeta;
                vz = currForce.z / zeta;
            }
                    
            nextPosition.x = currPosition.x + vx * dt;
            nextPosition.y = currPosition.y + vy * dt;
            nextPosition.z = currPosition.z + vz * dt;
            
            if(std::isnan(nextPosition.x))
            {
                throw std::runtime_error("Particle " + std::to_string(i) + " has a NaN position in the x direction: " + std::to_string(nextPosition.x));
            }
            
            nextPosition = simBox.periodicWrap(nextPosition);
            
            pinfo.setPosByIndex(i, nextPosition);
        }
	}
}

void Simulation::putAllParticlesInGrid()
{
    particleGrid.clearGrid();
    
    int numParticles = pinfo.getNumParticles();
    for(int i = numParticles-1; i > -1; i--)
    {
        int tempTag = pinfo.getTagAtIndex(i);

        if(tempTag >= 0)
        {
            int newBin = particleGrid.putInGrid(tempTag, pinfo.getPosByIndex(i));
        }
    }

    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        pinfo.setNeighborsByIndex(p, particleGrid.getNeighboringTags(pinfo.getTagAtIndex(p)));
    }

    std::vector <int> possibleBondNeighborTypes = {0,1}; 
    if(excludedVolOn)
    {
        #pragma omp parallel for
        for(int b = 0; b < pinfo.getNumBonds(); b++)
        {
            pinfo.resetBondNeighborTagsByIndex(b, possibleBondNeighborTypes);
            
            Bond bondB = pinfo.getBondByIndex(b);
            Coordinate b_i = pinfo.getPos(bondB.i);
            Coordinate b_j = pinfo.getPos(bondB.j);
                
            std::vector <int> newNeighboringTags = pinfo.getBondNeighborTagsByIndex(b);

            pinfo.setBondNeighborTagsByIndex(b, newNeighboringTags);
        }
    }
}

//Compartmentalizing the cell to simulate intrinsic ring and nuclear piston:
std::vector<int> Simulation::compartmentalizeUpperSys(int a, int b, int c, int d)
{

    int StartingBead = a; 

    int p = StartingBead;
    int direction = 1; 
    std::vector<int>perimeter_list; 


    while(perimeter_list.size() == 0 || p != StartingBead){ 

        if(p == a){ 

            perimeter_list.push_back(p);
            p = b;
            direction = -1;
        }
        else if(p == c){

            perimeter_list.push_back(p);
            p = d;
            direction = 1;
        }

        {
            perimeter_list.push_back(p);
        }
        
        p = p + direction;

        if(p >= (numNucBeads + numCMBeads)){ 
            p = numNucBeads;
        }
        else if(p == -1){
            p = numNucBeads - 1;
        }

    }

    return perimeter_list;
}

std::vector<int> Simulation::compartmentalizeLowerSys(int a, int b, int c, int d)
{
    
    int StartingBead = d;

    int p = StartingBead;
    int direction = 1;
    std::vector<int>perimeter_list;


    while(perimeter_list.size() == 0 || p != StartingBead){
        if(p == d){

            perimeter_list.push_back(p);
            p = c; 
            direction = -1;
        }
        else if(p == b){

            perimeter_list.push_back(p);
            p = a;
            direction = 1;
        }

        {
            perimeter_list.push_back(p);
        }
        
        p = p + direction;

        if(p >= (numNucBeads + numCMBeads)){
            p = numNucBeads;
        }
        else if(p == -1){
            p = numNucBeads - 1;
        }

    }

    return perimeter_list;
}

double Simulation::MeasureArea(std::vector<int> perimeter_list) //Takes perimeter_list made in compartmentalizeUpper(Lower)System and measures the compartment's area.
{
    double Area = 0.0;

    #pragma omp parallel for reduction (+:Area)
    for(int i = 0; i<perimeter_list.size(); i++)
    {

        int p = perimeter_list[i];
        Coordinate currPosition = pinfo.getPos(p);

        int a, b;
        if(i == perimeter_list.size()-1)
        {
            a = perimeter_list[i];
            b = perimeter_list[0];

        }
        else
        {
            a = perimeter_list[i];
            b = perimeter_list[i+1];

        }

                Area += MeasureAreaTriangle(pinfo.getPos(a),pinfo.getPos(b));
    }
    return Area;

}
double Simulation::conserveAreaOfCompartments(std::vector<int> perimeter_list, double EquilibriumArea, double K_Area, int forceType) //Conserves area of perimenter_list with given initial area
{
   double currArea = MeasureArea(perimeter_list);
   double avgAreaForce = 0.0;
   int countArea = 0;

   double AreaForcePrefactor = K_Area*(currArea - EquilibriumArea) / EquilibriumArea;
   #pragma omp parallel for reduction (+:avgAreaForce,countArea)
   for(int i = 0; i<perimeter_list.size(); i++)
   {
       int a, b, c;
       if(i == perimeter_list.size() - 2)
       {
            a = perimeter_list[i];
            b = perimeter_list[i+1];
            c = perimeter_list[0];
       }
       else if(i == perimeter_list.size()-1)
       {
            a = perimeter_list[i];
            b = perimeter_list[0];
            c = perimeter_list[1];
       }
       else
       {
            a = perimeter_list[i];
            b = perimeter_list[i+1];
            c = perimeter_list[i+2];
       }
     

        if(MeasureAreaTriangle(pinfo.getPos(a), pinfo.getPos(b)) != 0.0 && MeasureAreaTriangle(pinfo.getPos(b), pinfo.getPos(c)) != 0.0)
        {

            auto forceB = Coordinate {-AreaForcePrefactor*AreaConservationDerivativeXb(pinfo.getPos(a),pinfo.getPos(b),pinfo.getPos(c)), -AreaForcePrefactor*AreaConservationDerivativeYb(pinfo.getPos(a),pinfo.getPos(b),pinfo.getPos(c)), 0.0};
            avgAreaForce += forceB.getMagnitude();
            countArea++;

            tempForces[omp_get_thread_num()][forceType][pinfo.getIndexOfTag(b)] = tempForces[omp_get_thread_num()][forceType][pinfo.getIndexOfTag(b)] + forceB;
        }
   }

   return avgAreaForce / (double)countArea;

}

double Simulation::conserveAreaOfCompartmentsRing(std::vector<int> perimeter_list, std::vector<int> perimeter_list_two, double EquilibriumArea, double K_Area, int forceType)
{
   double currArea = MeasureArea(perimeter_list);
   double currAreaTwo = MeasureArea(perimeter_list_two); 
   double avgAreaForce = 0.0;
   int countArea = 0; 

   double AreaForcePrefactor1 = K_Area*(currArea - currAreaTwo - EquilibriumArea) / EquilibriumArea; 
   #pragma omp parallel for reduction (+:avgAreaForce,countArea) 
   for(int i = 0; i<perimeter_list.size(); i++) 
   {
       int a, b, c;
       if(i == perimeter_list.size() - 2)
       {
            a = perimeter_list[i];
            b = perimeter_list[i+1];
            c = perimeter_list[0];
       }
       else if(i == perimeter_list.size()-1)
       {
            a = perimeter_list[i];
            b = perimeter_list[0];
            c = perimeter_list[1];
       }
       else
       {
            a = perimeter_list[i];
            b = perimeter_list[i+1];
            c = perimeter_list[i+2];
       }

        if(MeasureAreaTriangle(pinfo.getPos(a), pinfo.getPos(b)) != 0.0 && MeasureAreaTriangle(pinfo.getPos(b), pinfo.getPos(c)) != 0.0)
        {

            auto forceB = Coordinate {-AreaForcePrefactor1*AreaConservationDerivativeXb(pinfo.getPos(a),pinfo.getPos(b),pinfo.getPos(c)), -AreaForcePrefactor1*AreaConservationDerivativeYb(pinfo.getPos(a),pinfo.getPos(b),pinfo.getPos(c)), 0.0};
            avgAreaForce += forceB.getMagnitude();

			countArea++;

            tempForces[omp_get_thread_num()][forceType][pinfo.getIndexOfTag(b)] = tempForces[omp_get_thread_num()][forceType][pinfo.getIndexOfTag(b)] + forceB;
        }
        else 
        {
            std::cout << "zero area? " << b << std::endl;
        }
   }

   double AreaForcePrefactor2 = -K_Area*(currArea - currAreaTwo - EquilibriumArea) / EquilibriumArea;
   #pragma omp parallel for reduction (+:avgAreaForce,countArea)
   for(int i = 0; i<perimeter_list_two.size(); i++)
   {
       int a, b, c;
       if(i == perimeter_list_two.size() - 2)
       {
            a = perimeter_list_two[i];
            b = perimeter_list_two[i+1];
            c = perimeter_list_two[0];
       }
       else if(i == perimeter_list_two.size()-1)
       {
            a = perimeter_list_two[i];
            b = perimeter_list_two[0];
            c = perimeter_list_two[1];
       }
       else
       {
            a = perimeter_list_two[i];
            b = perimeter_list_two[i+1];
            c = perimeter_list_two[i+2];
       }

        if(MeasureAreaTriangle(pinfo.getPos(a), pinfo.getPos(b)) != 0.0 && MeasureAreaTriangle(pinfo.getPos(b), pinfo.getPos(c)) != 0.0)
        {

            auto forceB = Coordinate {-AreaForcePrefactor2*AreaConservationDerivativeXb(pinfo.getPos(a),pinfo.getPos(b),pinfo.getPos(c)), -AreaForcePrefactor2*AreaConservationDerivativeYb(pinfo.getPos(a),pinfo.getPos(b),pinfo.getPos(c)), 0.0};
            avgAreaForce += forceB.getMagnitude(); 

			countArea++;

            tempForces[omp_get_thread_num()][forceType][pinfo.getIndexOfTag(b)] = tempForces[omp_get_thread_num()][forceType][pinfo.getIndexOfTag(b)] + forceB;
        }
        else
        {
            std::cout << "zero area? " << b << std::endl;
        }
   }
   return avgAreaForce / (double)countArea;
}

double Simulation::calcDeterminant(Coordinate firstRow, Coordinate secondRow, Coordinate thirdRow)
{
    double x = secondRow.y*thirdRow.z - thirdRow.y*secondRow.z;
    double y = secondRow.x*thirdRow.z - thirdRow.x*secondRow.z;
    double z = secondRow.x*thirdRow.y - thirdRow.x*secondRow.y;

    return firstRow.x*x - firstRow.y*y + firstRow.z*z;
}

Coordinate Simulation::closestDistanceBetweenLines(Coordinate a0, Coordinate a1, Coordinate b0, Coordinate b1, Coordinate &aClose, Coordinate &bClose)
{
    
    Coordinate d_kl;

    Coordinate A = a1 - a0;
    Coordinate B = b1 - b0;
    double magA = A.getMagnitude();
    double magB = B.getMagnitude();

    Coordinate _A = A / magA;
    Coordinate _B = B / magB;

    Coordinate cross = _A.crossProduct(_B);
    double denom = cross.getMagnitude() * cross.getMagnitude();

    if (denom == 0)
    {
        double d0 = _A * (b0-a0);

        double d1 = _A * (b1-a0);

        if (d0 <= 0 && 0 >= d1)
        {
            if (fabs(d0) < fabs(d1))
            {

                aClose = a0;
                bClose = b0;

                d_kl = (a0 - b0);

                return d_kl;
            }

            aClose = a0;
            bClose = b1;

            d_kl = (a0 - b1);

            return d_kl;
        }

        else if (d0 >= magA && magA <= d1)
        {
            if (fabs(d0) < fabs(d1))
            {

                aClose = a1;
                bClose = b0;

                d_kl = (a1 - b0);

                return d_kl;
            }

            aClose = a1;
            bClose = b1;

            d_kl = (a1 - b1);

            return d_kl;
        }

        aClose = a0 + A*0.5;
        bClose = b0 + B*0.5;

        d_kl = (((d0 * _A) + a0) - b0);
        return d_kl;
    }

    Coordinate t = (b0 - a0);
    double detA = calcDeterminant(t, _B, cross);
    double detB = calcDeterminant(t, _A, cross);

    double t0 = detA / denom;
    double t1 = detB / denom;

    Coordinate pA = a0 + (_A * t0); 
    Coordinate pB = b0 + (_B * t1); 

    if (t0 < 0)
        pA = a0;
    else if (t0 > magA)
        pA = a1;

    if (t1 < 0)
        pB = b0;
    else if (t1 > magB)
        pB = b1;

    double dot;
    if (t0 < 0 || t0 > magA)
    {
        dot = _B * (pA - b0);
        if (dot < 0)
            dot = 0;
        else if (dot > magB)
            dot = magB;
        pB = b0 + (_B * dot);
    }
    if (t1 < 0 || t1 > magB)
    {
        dot = _A*(pB - a0);
        if (dot < 0)
            dot = 0;
        else if (dot > magA)
            dot = magA;
        pA = a0 + (_A * dot);
    }

    aClose = pA;
    bClose = pB;

    d_kl = (pA - pB);
    return d_kl;
}

// Function that simulates cell turning and allows for small changes in direction of migration by shifting bead types as a repolarization mechanism:
void Simulation::shiftBeadTypes(int shiftAmount)
{
    std::vector <int> beadTypeList;

    for(int p = numNucBeads; p < numNucBeads + numCMBeads; p++)
    {
        beadTypeList.push_back(pinfo.getBeadToFilament(p));
    }

    int count = 0;
    int currPos = numNucBeads + shiftAmount;

    while(count < numCMBeads)
    {
        while(currPos <= numNucBeads)
        {
            currPos = currPos + numCMBeads;
        }

        while(currPos >= numNucBeads + numCMBeads)
        {
            currPos = currPos - numCMBeads;
        }

        pinfo.setBeadToFilamentByIndex(currPos, beadTypeList[count]);

        count++;
        currPos++;
    }
}				
																																   
double Simulation::calcForces()
{
    std::mt19937& threadEngine = generators[omp_get_thread_num()];
    int numParticles = pinfo.getNumParticles();
    {
        for(int thread = 0; thread < omp_get_max_threads(); thread++)
        {
            int oldSize = tempForces[thread][0].size();
            
            for(int k = 0; k < num_force_types; k++)
            {
                tempForces[thread][k].resize(pinfo.getNumParticles());
                
                for(int i = 0; i < pinfo.getNumParticles(); i++)
                {
                    tempForces[thread][k][i] = Coordinate {0.0, 0.0, 0.0};
                }
            }
            
            
            tempStresses[thread].resize(pinfo.getNumParticles());
            
            for(int i = oldSize; i < pinfo.getNumParticles(); i++)
            {

                tempStresses[thread][i] = Coordinate {0.0, 0.0, 0.0};
            }
        }
    }
    
    double currArea = MeasureAreaNucleus();
    double currAreaCM = MeasureAreaCM();
    #pragma omp parallel for
    for(int p = 0; p < numParticles; p++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(p);

        int currTag = pinfo.getTagAtIndex(p);
        std::vector <int> iNeighbors = pinfo.getNeighbors(currTag);

        std::vector <struct HalfBond> tempBondList = pinfo.getBondTagsByIndex(p);
		double phase = M_PI;	
        // Stretching forces from bonds:						   
        for(int b = 0; b < tempBondList.size(); b++)
        {
            Bond firstBond = pinfo.getBond(tempBondList[b].bondTag);

            if(currTag < tempBondList[b].j)
            {

                    Coordinate ri = pinfo.getPos(currTag);
                    Coordinate rj = pinfo.getPos(tempBondList[b].j);
                    Coordinate rij = simBox.calcDisplacement(ri, rj);


                    if(firstBond.bondTypeIndex >= 0)
                    {
                        double rijMag = rij.getMagnitude();

                        double bondDestructionTime = firstBond.destructionTime;
						
						Coordinate rijUnit = rij / rijMag;

						Coordinate bondForceVector = bondTypes[firstBond.bondTypeIndex]->calcForce(rijMag)*rijUnit;

						int currPosJ = pinfo.getIndexOfTag(tempBondList[b].j);
                        
                        // Bleb expansion and retraction is simulated as oscillations in the spring constant of bonds in the flanking regions of the leading edge:
						if(blebbing){
							
							double bondForceMag = bondForceVector.getMagnitude(); 
							Coordinate bondForceUnit = bondForceVector.getUnitCoord();

							if (pinfo.getBeadToFilamentByIndex(p) >= 2 && pinfo.getBeadToFilamentByIndex(p) <= 4)
							{
								double K_Stretch_Osc = K_Stretch;
								double K_Stretch_min = 10000.0; 
								double A = (K_Stretch_min + K_Stretch)/(K_Stretch - K_Stretch_min); 
	
								if(pinfo.getBeadToFilamentByIndex(p) == 2)
								{
									K_Stretch_Osc = (K_Stretch/(A+1.0))*(A+cos(OMEGA*currSimTime + phase)); 
								}
								else if(pinfo.getBeadToFilamentByIndex(p) == 4)
								{
									K_Stretch_Osc = (K_Stretch/(A+1.0))*(A+cos(OMEGA*currSimTime + (phase+M_PI)));
								}

								bondForceMag = (K_Stretch_Osc/K_Stretch)*bondForceMag;

							}
								
							
							tempForces[omp_get_thread_num()][0][p] = tempForces[omp_get_thread_num()][0][p] + bondForceMag*bondForceUnit;
							tempForces[omp_get_thread_num()][0][currPosJ] = tempForces[omp_get_thread_num()][0][currPosJ] - bondForceMag*bondForceUnit;
						}
                        
						else{
							tempForces[omp_get_thread_num()][0][p] = tempForces[omp_get_thread_num()][0][p] + bondForceVector;
                            tempForces[omp_get_thread_num()][0][currPosJ] = tempForces[omp_get_thread_num()][0][currPosJ] - bondForceVector;
						}
						 
                    }
            }
        }

        std::vector <struct Angle> tempAngleList = pinfo.getAngleTagsByIndex(p); 
        //Bending forces from angles:
        for(int la = 0; la < tempAngleList.size(); la++)
        {
            if(currTag == tempAngleList[la].i)
            {
                if(currTag < numNucBeads)
                {

                    Angle ang = pinfo.getAngle(tempAngleList[la].angleTypeIndex);

                    Coordinate a = pinfo.getPos(ang.i);
                    Coordinate b = pinfo.getPos(ang.j);
                    Coordinate c = pinfo.getPos(ang.k);

                    Coordinate ab = simBox.calcDisplacement(a, b);
                    double abMag = ab.getMagnitude();
                    Coordinate abUnit = ab / abMag;

                    Coordinate cb = simBox.calcDisplacement(c, b);
                    double cbMag = cb.getMagnitude();
                    Coordinate cbUnit = cb / cbMag;

                    double averageMag = (abMag + cbMag) * 0.5;
                    double invAverageMag = 1.0 / averageMag;

                    double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;

                    Coordinate forceA = -kappa_Nuc * invAverageMag * (cb / cbMag / abMag + dotProd / cbMag * (-ab / (abMag*abMag*abMag)));
                    Coordinate forceC = -kappa_Nuc * invAverageMag * (ab / abMag / cbMag + dotProd / abMag * (-cb / (cbMag*cbMag*cbMag)));

                    Coordinate forceB = forceA + forceC;

                    int aIndex = pinfo.getIndexOfTag(ang.i);
                    int bIndex = pinfo.getIndexOfTag(ang.j);
                    int cIndex = pinfo.getIndexOfTag(ang.k);

                    tempForces[omp_get_thread_num()][1][aIndex] = tempForces[omp_get_thread_num()][1][aIndex] + forceA;
                    tempForces[omp_get_thread_num()][1][bIndex] = tempForces[omp_get_thread_num()][1][bIndex] - forceB;
                    tempForces[omp_get_thread_num()][1][cIndex] = tempForces[omp_get_thread_num()][1][cIndex] + forceC;
                }
                if(currTag >= numNucBeads)
                {
                    Angle ang = pinfo.getAngle(tempAngleList[la].angleTypeIndex);

                    Coordinate a = pinfo.getPos(ang.i);
                    Coordinate b = pinfo.getPos(ang.j);
                    Coordinate c = pinfo.getPos(ang.k);

                    Coordinate ab = simBox.calcDisplacement(a, b);
                    double abMag = ab.getMagnitude();
                    Coordinate abUnit = ab / abMag;

                    Coordinate cb = simBox.calcDisplacement(c, b);
                    double cbMag = cb.getMagnitude();
                    Coordinate cbUnit = cb / cbMag;

                    double averageMag = (abMag + cbMag) * 0.5;
                    double invAverageMag = 1.0 / averageMag;

                    double dotProd = ab.x*cb.x + ab.y*cb.y + ab.z*cb.z;

                    Coordinate forceA = -kappa * invAverageMag * (cb / cbMag / abMag + dotProd / cbMag * (-ab / (abMag*abMag*abMag)));
                    Coordinate forceC = -kappa * invAverageMag * (ab / abMag / cbMag + dotProd / abMag * (-cb / (cbMag*cbMag*cbMag)));

                    Coordinate forceB = forceA + forceC;

                    int aIndex = pinfo.getIndexOfTag(ang.i);
                    int bIndex = pinfo.getIndexOfTag(ang.j);
                    int cIndex = pinfo.getIndexOfTag(ang.k);

                    tempForces[omp_get_thread_num()][1][aIndex] = tempForces[omp_get_thread_num()][1][aIndex] + forceA;
                    tempForces[omp_get_thread_num()][1][bIndex] = tempForces[omp_get_thread_num()][1][bIndex] - forceB;
                    tempForces[omp_get_thread_num()][1][cIndex] = tempForces[omp_get_thread_num()][1][cIndex] + forceC;

                }
            }
        }
        
	}

    // Nucleus-centering force:
	if(numNucBeads > 0) 																			 
	{
		Coordinate pos_nuc_sum = Coordinate {0.0,0.0,0.0};
		for(int p = 0; p < numNucBeads; p++)
		{
			pos_nuc_sum = pos_nuc_sum + pinfo.getPosByIndex(p);
		}

		r_nuc_com = pos_nuc_sum/(1.0*numNucBeads);

		Coordinate pos_CM_sum = Coordinate {0.0,0.0,0.0};
		for(int p = numNucBeads; p < (numNucBeads + numCMBeads); p++) 
		{
			pos_CM_sum = pos_CM_sum + pinfo.getPosByIndex(p);
		}

		Coordinate r_CM_com = pos_CM_sum/(1.0*numCMBeads);

		Coordinate r_disp = r_CM_com - r_nuc_com;
		Coordinate bond_com_Force = bondTypes[2] -> calcForce(r_disp.getMagnitude())*r_disp.getUnitCoord();

		double number_of_nuclear_membrane_beads = 1.0*numNucBeads;

		for(int p = 0; p < numNucBeads; p++)
		{

			tempForces[omp_get_thread_num()][4][p] = tempForces[omp_get_thread_num()][4][p] - bond_com_Force/number_of_nuclear_membrane_beads;
		}
		for(int p = numNucBeads; p < (numNucBeads + numCMBeads); p++) 
		{

			tempForces[omp_get_thread_num()][4][p] = tempForces[omp_get_thread_num()][4][p] + bond_com_Force/(1.0*numCMBeads); 
		}
	}
    
    if(excludedVolOn)
    {
 
        int currNumBonds = pinfo.getNumBonds(); 
        
        std::vector <int> countPerThread (omp_get_max_threads());

        std::vector<std::tuple<int, int, double>> NeighborPairList;
        std::vector<std::tuple<int, int, double>> obstaclePairList;
        #pragma omp parallel for schedule(static, 1)
        for(int b = 0; b < currNumBonds; b++)
        {            
            std::vector <int> currNeighborBondTags = pinfo.getBondNeighborTagsByIndex(b);
            Bond firstBond = pinfo.getBondByIndex(b);
            
            int currPosI, currPosJ;
            Coordinate ri, rj;
            
            try
            {
                currPosI = pinfo.getIndexOfTag(firstBond.i);
                currPosJ = pinfo.getIndexOfTag(firstBond.j);
                
                ri = pinfo.getPosByIndex(currPosI);
                rj = pinfo.getPosByIndex(currPosJ);
            }
            catch (const std::runtime_error& error)
            {
                #pragma omp critical
                {
                    cout << b << " " << firstBond.i << " " << firstBond.j << " " << currNumBonds << " " << currPosI << " " << currPosJ << endl;
                    exit(0);
                }
            }

            Coordinate R_k = simBox.calcDisplacement(ri, rj);
            double R_kSquared = R_k*R_k;
            
            Coordinate P_k = 0.5 * R_k + rj;
            P_k = simBox.periodicWrap(P_k);

            for(int neigh = 0; neigh < currNeighborBondTags.size(); neigh++)
            {
                if(currNeighborBondTags[neigh] < 0)
                {
                    cout<<currNeighborBondTags[neigh]<< endl;
                    }
                if(currNeighborBondTags[neigh] > pinfo.getNumParticles())
                {
                    cout<<currNeighborBondTags[neigh]<< endl;
                    }

                if(currNeighborBondTags[neigh] < pinfo.getBondTagAtPos(b))
                {

                    Coordinate tempStress = {0.0, 0.0, 0.0};

                    Bond secondBond = pinfo.getBond(currNeighborBondTags[neigh]);

                    int currPosI2, currPosJ2;
                    Coordinate ri2, rj2;

                    try
                    {
                        currPosI2 = pinfo.getIndexOfTag(secondBond.i);
                        currPosJ2 = pinfo.getIndexOfTag(secondBond.j);

                        ri2 = pinfo.getPosByIndex(currPosI2);
                        rj2 = pinfo.getPosByIndex(currPosJ2);
                    }
                    catch(const runtime_error& error)
                    {
                        #pragma omp critical
                        {
                            cout << "second" << " " << secondBond.i << " " << secondBond.j << " " << currPosI2 << " " << currPosJ2 << endl;
                            exit(0);
                        }
                    }

                    Coordinate R_l = simBox.calcDisplacement(ri2, rj2);

                    Coordinate P_l = simBox.periodicWrap(0.5 * R_l + rj2);

                    double R_lSquared = R_l*R_l;

                    double R_klSquared = R_k*R_l;

                    Coordinate closestK, closestL;
                    Coordinate d_kl = closestDistanceBetweenLines(ri, rj, ri2, rj2, closestK, closestL); 
                    double minDist = d_kl.getMagnitude();

                    if(minDist < thresholdDistCB) // Checking for neighboring filament beads that are within d_seal distance of each other to create a neighbor list for sealing:
                    {
                        if(firstBond.i >= numNucBeads && secondBond.i < numNucBeads)
                        {
                            double pairDistance = simBox.calcDistance(pinfo.getPos(firstBond.i),pinfo.getPos(secondBond.i));
                            if(pairDistance < d_seal)
                            {
                                # pragma omp critical
                                {
                                NeighborPairList.push_back(std::make_tuple(firstBond.i, secondBond.i, pairDistance));
                                }
                            }
                            pairDistance = simBox.calcDistance(pinfo.getPos(firstBond.i),pinfo.getPos(secondBond.j));
                            if(pairDistance < d_seal)
                            {
                                # pragma omp critical
                                {
                                NeighborPairList.push_back(std::make_tuple(firstBond.i, secondBond.j, pairDistance));
                                }
                            }
                            pairDistance = simBox.calcDistance(pinfo.getPos(firstBond.j),pinfo.getPos(secondBond.i));
                            if(pairDistance < d_seal)
                            {
                                # pragma omp critical
                                {
                                NeighborPairList.push_back(std::make_tuple(firstBond.j, secondBond.i, pairDistance));
                                }
                            }
                            pairDistance = simBox.calcDistance(pinfo.getPos(firstBond.j),pinfo.getPos(secondBond.j));
                            if(pairDistance < d_seal)
                            {
                                # pragma omp critical
                                {
                                NeighborPairList.push_back(std::make_tuple(firstBond.j, secondBond.j, pairDistance));
                                }
                            }
                        }
                    }

                    // Excluded volume forces between beads of cell and nucleus to prevent overlap:
                    if(minDist < d_Excl) 
                    {

                        double forceMag = 0.0;

						forceMag = kExcl_CM * (d_Excl - minDist);

                        Coordinate Fr = forceMag * d_kl.getUnitCoord();

                        double length_k = R_k.getMagnitude();

                        double firstSegmentLength = simBox.calcDistance(closestK, rj);
                        double secondSegmentLength = length_k - firstSegmentLength;

                        Coordinate Falpha = firstSegmentLength / length_k * (Fr);
                        Coordinate Fbeta = secondSegmentLength / length_k * (Fr);

                        double length_l = R_l.getMagnitude();

                        firstSegmentLength = simBox.calcDistance(closestL, rj2);
                        secondSegmentLength = length_l - firstSegmentLength;

                        Coordinate Falpha2 = firstSegmentLength / length_l * (-Fr);
                        Coordinate Fbeta2 = secondSegmentLength / length_l * (-Fr);

                        Falpha.z = 0.0;
                        Fbeta.z = 0.0;
                        Falpha2.z = 0.0;
                        Fbeta2.z = 0.0;
                        tempForces[omp_get_thread_num()][2][currPosI] = tempForces[omp_get_thread_num()][2][currPosI] + Falpha;
                        tempForces[omp_get_thread_num()][2][currPosJ] = tempForces[omp_get_thread_num()][2][currPosJ] + Fbeta;
                        tempForces[omp_get_thread_num()][2][currPosI2] = tempForces[omp_get_thread_num()][2][currPosI2] + Falpha2;
                        tempForces[omp_get_thread_num()][2][currPosJ2] = tempForces[omp_get_thread_num()][2][currPosJ2] + Fbeta2;

                    }
                }
            }
            if(pinfo.getBeadToFilament(firstBond.i) != 0)
            {
                std::vector <int> currObstacleNeighbors_i = pinfo.getNeighbors(firstBond.i); 
                std::vector <int> currObstacleNeighbors_j = pinfo.getNeighbors(firstBond.j); 

                std::set <int> currObstacleSet;
                for(int i = 0; i < currObstacleNeighbors_i.size(); i++)
                {
                    currObstacleSet.insert(currObstacleNeighbors_i[i]);
                }
                for(int i = 0; i < currObstacleNeighbors_j.size(); i++)
                {
                    currObstacleSet.insert(currObstacleNeighbors_j[i]);
                }

                for(auto it = currObstacleSet.begin(); it != currObstacleSet.end(); it++) 
                {
                    int p_tag = *it;
                    int p = pinfo.getIndexOfTag(p_tag);

                    int tempType = pinfo.getBeadToFilamentByIndex(p);

                    if(tempType == 1)
                    {

                        Coordinate static_pos = pinfo.getPosByIndex(p);

                        Coordinate static_pos_A = Coordinate {static_pos.x, static_pos.y, -1.0};
                        Coordinate static_pos_B = Coordinate {static_pos.x, static_pos.y, 1.0};

                        Coordinate R_l = simBox.calcDisplacement(static_pos_A, static_pos_B);
                        Coordinate P_l = simBox.periodicWrap(0.5 * R_l + static_pos_B);

                        double R_lSquared = R_l*R_l;

                        double R_klSquared = R_k*R_l;

                        Coordinate closestK, closestL;
                        Coordinate d_kl = closestDistanceBetweenLines(ri, rj, static_pos_A, static_pos_B, closestK, closestL);
                        double minDist = d_kl.getMagnitude();

                        if(minDist < thresholdDistCB)
                        {
                            double pairDistance = simBox.calcDistance(ri, static_pos);

                            if(pairDistance < d_seal)
                            {
                                # pragma omp critical
                                {
                                obstaclePairList.push_back(std::make_tuple(firstBond.i, p, pairDistance));
                                }
                            }
                            pairDistance = simBox.calcDistance(rj, static_pos);
                            if(pairDistance < d_seal)
                            {
                                # pragma omp critical
                                {
                                obstaclePairList.push_back(std::make_tuple(firstBond.j, p, pairDistance));
                                }
                            }
                        }

                        // Excluded volume forces between beads of cell and obstacles:
                        if(minDist < d_Excl) 
                        {
                            double forceMag = 0.0;

							forceMag = kExcl_Obs * (d_Excl - minDist);

                            Coordinate Fr = forceMag * d_kl.getUnitCoord();

                            double length_k = R_k.getMagnitude();

                            double firstSegmentLength = simBox.calcDistance(closestK, rj);
                            double secondSegmentLength = length_k - firstSegmentLength;

                            Coordinate Falpha = firstSegmentLength / length_k * (Fr);
                            Coordinate Fbeta = secondSegmentLength / length_k * (Fr);

                            double length_l = R_l.getMagnitude();

                            firstSegmentLength = simBox.calcDistance(closestL, static_pos_B);
                            secondSegmentLength = length_l - firstSegmentLength;

                            Coordinate Falpha2 = firstSegmentLength / length_l * (-Fr);
                            Coordinate Fbeta2 = secondSegmentLength / length_l * (-Fr);

                            Falpha.z = 0.0;
                            Fbeta.z = 0.0;


                            tempForces[omp_get_thread_num()][9][currPosI] = tempForces[omp_get_thread_num()][9][currPosI] + Falpha;
                            tempForces[omp_get_thread_num()][9][currPosJ] = tempForces[omp_get_thread_num()][9][currPosJ] + Fbeta;
                        }
                    }

                }

            }
        }

        areaForceNucleus = conserveAreaOfCompartments(Nuc_List, AreaNucleusZero, K_Area_nuc, 3); //area conservation force due to nucleus

        // Mechanism of Calcium ion release from nuclear membrane and perinuclear ER when nucleus is highly compressed.
        // Calcium ion release triggers an increase in the cortical tension of the cell.
        if(areaForceNucleus > 3.0*(F_Mag_Front/Ratio)){ 
            double a_CalciumIon  = 1.5;
            F_Mag_back = a_CalciumIon*(F_Mag_Front/Ratio);
        }
        else{
            F_Mag_back = F_Mag_Front/Ratio;
        }

        if(CB_On)
        {
            int tempCBM_1 = -1; //points to membrane bead making the first seal, -1 if there is no seal
            int tempCBN_1 = -1; //points to nucleus bead making the first seal, -1 otherwise
            int tempCBN_2 = -1;
            int tempCBM_2 = -1;
            
            //sorts NPL from smallest to largest
            std::sort(std::begin(NeighborPairList), std::end(NeighborPairList), [](tuple<int, int, double> const &t1, tuple<int, int, double> const &t2)
                        {
                            return get<2>(t1) < get<2>(t2); 
                        }
                    );
            std::sort(std::begin(obstaclePairList), std::end(obstaclePairList), [](tuple<int, int, double> const &t1, tuple<int, int, double> const &t2)
                        {
                            return get<2>(t1) < get<2>(t2);
                        }
                    );
            for(int o = 0; o < obstaclePairList.size(); o++) //iterate thru OPS and see if 1st & 2nd elements are also in NPL.
            {
                int PM_bead_obstacle = std::get<0>(obstaclePairList[o]);
                int obstacleBeadIndex = std::get<1>(obstaclePairList[o]);


                for(int i = 0; i < NeighborPairList.size(); i++) //Go through all nearest neighbors if there is at least 1 seal
                {
                    int PM_bead = std::get<0>(NeighborPairList[i]);
                    int Nuclear_bead = std::get<1>(NeighborPairList[i]);

                    if(PM_bead == PM_bead_obstacle)
                    {

                        if(tempCBM_1 == -1)
                        {   
                            // When nearest neighbors form, assign first seal
                            //sort NPL, assign seal for smallest pair distance
                            //https://stackoverflow.com/questions/23030267/custom-sorting-a-vector-of-tuples
                            tempCBM_1 = std::get<0>(NeighborPairList[i]);
                            
                            tempCBN_1 = std::get<1>(NeighborPairList[i]);

                        }

                        else if(tempCBM_2 == -1) //If 1st seal is found search for a 2nd seal.
                        {
                            int tempX = std::get<0>(NeighborPairList[i]);
                            int tempY = tempCBM_1;

                            int circDist = std::min(fabs(tempX-tempY), (numCMBeads - 1) - (fabs(tempX-tempY)-1));

                            Coordinate sealDisplacement1 = simBox.calcDisplacement(pinfo.getPos(tempCBM_1), pinfo.getPos(tempCBN_1));
                            Coordinate sealDisplacementHat1 = sealDisplacement1.getUnitCoord(); 

                            Coordinate sealDisplacement2 = simBox.calcDisplacement(pinfo.getPos(std::get<0>(NeighborPairList[i])), pinfo.getPos(std::get<1>(NeighborPairList[i])));
                            Coordinate sealDisplacementHat2 = sealDisplacement2.getUnitCoord();

                            double dotProd = sealDisplacement1*sealDisplacement2;

                            if(circDist > (int)(0.3*numNucBeads) && dotProd < 0.0)
                            {
                                tempCBN_2 = std::get<1>(NeighborPairList[i]);

                                tempCBM_2 = std::get<0>(NeighborPairList[i]);

                                o = obstaclePairList.size();
                                break;
                            }
                        }
                    }

                }
            }
            // Once a seal is formed, it is maintained for a certain time period (timeLag) to avoid rapid fluctuations in the sealing region.
            // If the seal breaks, then the sealing region is reset.
            if(CBN_2 == -1 || currSimTime - timeSRS > timeLag) 
            {
                CBM_1 = tempCBM_1;
                CBN_1 = tempCBN_1;
                CBN_2 = tempCBN_2;
                CBM_2 = tempCBM_2;
                timeSRS = currSimTime;
            }

            if(tempCBN_2 == -1) 
            {
                CBM_1 = -1;
                CBN_1 = -1;
                CBN_2 = -1;
                CBM_2 = -1;
            }

            if(CBN_2 == -1)  //If there are no seals, keep conserving the cortex annulus as one compartment.
            {

                areaForceCortex = conserveAreaOfCompartmentsRing(CM_List, Nuc_List, AreaCMZero - AreaNucleusZero, K_Area_CM, 7); 
                EquilibriumAreaUpper = 0.0; 
                EquilibriumAreaLower = 0.0;

            }
            else // If there are seals, compartmentalize cell cortex into upper and lower compartments and conserve their areas separately.
            {

                upperCompartment = compartmentalizeUpperSys(CBM_1, CBN_1, CBN_2, CBM_2);
                lowerCompartment = compartmentalizeLowerSys(CBM_1, CBN_1, CBN_2, CBM_2);


                for(int i = 0; i < upperCompartment.size(); i++)
                {
                    if(upperCompartment[i] == (int)(0.5*numNucBeads)) 
                    {
                        std::vector<int> tmpCompartment = upperCompartment; 
                        upperCompartment = lowerCompartment;
                        lowerCompartment = tmpCompartment; 
                        break;
                    }
                }

                if(EquilibriumAreaUpper == 0.0)
                {
                    EquilibriumAreaUpper = MeasureArea(upperCompartment);
                    EquilibriumAreaLower = MeasureArea(lowerCompartment);
                    double TempCompartmentAreaUpper = (AreaCMZero - AreaNucleusZero)*(EquilibriumAreaUpper/(EquilibriumAreaUpper + EquilibriumAreaLower)); 
                    double TempCompartmentAreaLower = (AreaCMZero - AreaNucleusZero)*(EquilibriumAreaLower/(EquilibriumAreaUpper + EquilibriumAreaLower)); 
                    EquilibriumAreaUpper = TempCompartmentAreaUpper;
                    EquilibriumAreaLower = TempCompartmentAreaLower;

                }

                conserveAreaOfCompartments(upperCompartment, EquilibriumAreaUpper, K_Area_CM, 7);

                conserveAreaOfCompartments(lowerCompartment, EquilibriumAreaLower, K_Area_CM, 8);

            }
        }
        else
        {
            areaForceCortex = conserveAreaOfCompartmentsRing(CM_List, Nuc_List, AreaCMZero - AreaNucleusZero, K_Area_CM, 7);
        }
    }

    double simBoxYMax = simBox.getBoxLength(1)*0.5;

    #pragma omp parallel for
    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
  
        std::mt19937& threadEngine = generators[omp_get_thread_num()];
        
        Coordinate currPos = pinfo.getPosByIndex(p);
            
        if(temperature > 0.0 && thermalForcesOn == true)
        {

            Coordinate c = {thermal_force*gaussianDis(threadEngine), thermal_force*gaussianDis(threadEngine), 0.0};
            pinfo.setForceByIndex(p, c);
        }
        else
        {
            pinfo.setForceByIndex(p, Coordinate {0.0, 0.0, 0.0});
        }

        double phase = M_PI;

        // Contractile force on cell cortex beads:
        if (pinfo.getBeadToFilamentByIndex(p) == 5){

            int p_1 = p - 1;
            int p_plus1 = p + 1;
            if(p_1 == (numNucBeads - 1)){

            p_1 = numNucBeads + numCMBeads - 1;

            }

            if(p_plus1 == (numNucBeads + numCMBeads)){

            p_plus1 = numNucBeads;

            }
            Coordinate P_minus1 = pinfo.getPosByIndex(p_1);
            Coordinate P_plus1 = pinfo.getPosByIndex(p_plus1);
            Coordinate r =  simBox.calcDisplacement(P_minus1, P_plus1);

            Coordinate F_unit = Coordinate (-r.y, r.x, 0.0); 
            F_unit = F_unit.getUnitCoord();

            tempForces[omp_get_thread_num()][6][p] = F_unit*F_Mag_back + tempForces[omp_get_thread_num()][6][p];
            
        }

        pinfo.setStressByIndex(p, Coordinate {0.0, 0.0, 0.0}, 0);
        pinfo.setStressByIndex(p, Coordinate {0.0, 0.0, 0.0}, 1);

        Coordinate totalForce = {0.0, 0.0, 0.0};
        Coordinate totalStress = {0.0, 0.0, 0.0};

        bool ExcVolIsActing = false; 
        Coordinate ExcVolForce = Coordinate {0.0, 0.0, 0.0}; 
        for(int k = 0; k < num_force_types; k++) 
        {
            Coordinate tempSum = {0.0, 0.0, 0.0}; 
            for(int thread = 0; thread < omp_get_max_threads(); thread++)
            {
                Coordinate tempForce = tempForces[thread][k][p];
                tempSum = tempSum + tempForce;
            }

            if(k == 2) 
            {
                if(tempSum.getMagnitude() > 0.0)
                {
                    ExcVolIsActing = true; 
                    ExcVolForce = tempSum;
                }

            }

            if(k == 6 || k == 7) 
            {
                // Handling the case in which the cortical boundary closely approaches the nuclear boundary and the contractile and cortical area forces are locally reduced to zero in these regions.
                if(!ExcVolIsActing) 
                {
                    totalForce = totalForce + tempSum;
                }
                
                else
                {
                    
                    if(ExcVolForce.getMagnitude() < tempSum.getMagnitude())
                    {
                        
                        totalForce = totalForce + tempSum - ExcVolForce;
                        for(int threadEV = 0; threadEV < omp_get_max_threads(); threadEV++)
                        {
                            
                            tempForces[threadEV][2][p] = Coordinate {0.0, 0.0, 0.0};
                        }
                    }

                    for(int threadExcluded = 0; threadExcluded < omp_get_max_threads(); threadExcluded++) 
                    {
                        tempForces[threadExcluded][k][p] = Coordinate {0.0, 0.0, 0.0}; 
                    }
                }

            }

            else 
            {
                totalForce = totalForce + tempSum;

            }

        } 

        pinfo.addForceByIndex(p, totalForce);

        pinfo.addStressByIndex(p, totalStress, 0);
    }
    
    double totalLeadingEdgeForce = 0.0;

    return totalLeadingEdgeForce;
}

void Simulation::writeXYZFile(std::vector<std::string>& xyzArray)
{
    std::ostringstream ss; 
    ss << pinfo.getNumParticles() << "\n"; 
    ss << "t=" << currSimTime << "\n"; 

    for(int p = 0; p < pinfo.getNumParticles(); p++)
    {
        Coordinate tempCoordinate = pinfo.getPosByIndex(p);

        int output_type = pinfo.getBeadToFilamentByIndex(p);

        if(CBN_2 != -1 && (p == CBN_1 || p == CBN_2 || p == CBM_1 || p == CBM_2))
        {
            output_type = 6;
            
        }

        ss << output_type << " " << pinfo.getTagAtIndex(p) << " " << tempCoordinate.x << " " << tempCoordinate.y << " " << tempCoordinate.z << " ";

        for(int k = 0; k < num_force_types; k++)
        {
            Coordinate tmpSpecificForce = Coordinate (0.0, 0.0, 0.0);

            for(int thread = 0; thread < omp_get_max_threads(); thread++)
            {
                tmpSpecificForce = tmpSpecificForce + tempForces[thread][k][p];
            }

            ss << tmpSpecificForce.x << " " << tmpSpecificForce.y << " " << tmpSpecificForce.z << " ";
        }

        Coordinate tempForce = pinfo.getForceByIndex(p);
        ss << tempForce.x << " " << tempForce.y << " " << tempForce.z;

        ss << "\n";
    }

    xyzArray.push_back(ss.str());
}


void Simulation::writeBondFile(std::vector<std::string>& bondArray)
{
    std::ostringstream ss; 

    int numBonds = pinfo.getNumBonds();

    ss << numBonds << "\n";
    ss << "t=" << currSimTime << "\n";
    for(int i = 0; i < numBonds; i++)
    {
        Bond tempBond = pinfo.getBondByIndex(i);
        ss << tempBond.bondTypeIndex << " " << tempBond.i << " " << tempBond.j << "\n";
    }

    ss << "\n";
    bondArray.push_back(ss.str());
}

double Simulation::HeavisideTheta(double val)
{
    if(val > 0.0)
    {
        return 1.0;
    }
    else if(val < 0.0)
    {
        return 0.0;
    }
    else
    {
        return 0.5;
    }
}

double Simulation::DiracDelta(double val)
{
    return 0.0;
}

double Simulation::AreaConservationDerivativeXb(Coordinate a, Coordinate b, Coordinate c)
{
    double centerX = r_nuc_com.x;
    double centerY = r_nuc_com.y;
    
    double xa = a.x - centerX;
    double ya = a.y - centerY;
    double xb = b.x - centerX;
    double yb = b.y - centerY;
    double xc = c.x - centerX;
    double yc = c.y - centerY;
   
   double first_term =  xb*ya - xa*yb;
   double prefactor_first_term = HeavisideTheta(first_term)*2 - 1;
   

   double second_term = xc*yb - xb*yc;
   double prefactor_second_term = HeavisideTheta(second_term)*2 - 1;

   return ya*abs(xb*ya - xa*yb)*DiracDelta(xb*ya - xa*yb)
   - yc*abs(xc*yb - xb*yc)*DiracDelta(xc*yb - xb*yc)
   + prefactor_first_term*(ya*(-0.5 + HeavisideTheta(xb*ya - xa*yb)))
   - prefactor_second_term*(yc*(-0.5 + HeavisideTheta(xc*yb - xb*yc)));
}

double Simulation::AreaConservationDerivativeYb(Coordinate a, Coordinate b, Coordinate c)
{
    double centerX = r_nuc_com.x;
    double centerY = r_nuc_com.y;
    
    double xa = a.x - centerX;
    double ya = a.y - centerY;
    double xb = b.x - centerX;
    double yb = b.y - centerY;
    double xc = c.x - centerX;
    double yc = c.y - centerY;
   
   double first_term =  xb*ya - xa*yb;
   double prefactor_first_term = HeavisideTheta(first_term)*2 - 1;
   

   double second_term = xc*yb - xb*yc;
   double prefactor_second_term = HeavisideTheta(second_term)*2 - 1;
   
   return -xa*abs(xb*ya - xa*yb)*DiracDelta(xb*ya - xa*yb)
   + xc*abs(xc*yb - xb*yc)*DiracDelta(xc*yb - xb*yc)
   - prefactor_first_term*(xa*(-0.5 + 1.*HeavisideTheta(xb*ya - xa*yb)))
   + prefactor_second_term*(xc*(-0.5 + 1.*HeavisideTheta(xc*yb - xb*yc)));
   
}


double Simulation::MeasureAreaNucleus()
{
    double Area = 0.0;

    for(int p = 0; p < numNucBeads; p++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(p);
        
        int currTag = pinfo.getTagAtIndex(p);

        std::vector <struct Angle> tempAngleList = pinfo.getAngleTagsByIndex(p);
        for(int la = 0; la < tempAngleList.size(); la++)
        {
            if(currTag == tempAngleList[la].i)
            {
                Angle ang = pinfo.getAngle(tempAngleList[la].angleTypeIndex);
                
                Coordinate a = pinfo.getPos(ang.i);                
                Coordinate b = pinfo.getPos(ang.j);
                
                Area += MeasureAreaTriangle(a, b);
            }
        }
        
	}

    return Area;

}

double Simulation::MeasureAreaCM()
{
    double Area = 0.0;

    for(int p = numNucBeads; p < pinfo.getNumParticles(); p++)
    {
        Coordinate currPosition = pinfo.getPosByIndex(p);

        int currTag = pinfo.getTagAtIndex(p);

        std::vector <struct Angle> tempAngleList = pinfo.getAngleTagsByIndex(p);
        for(int la = 0; la < tempAngleList.size(); la++)
        {
            if(currTag == tempAngleList[la].i)
            {
                Angle ang = pinfo.getAngle(tempAngleList[la].angleTypeIndex);

                Coordinate a = pinfo.getPos(ang.i);
                Coordinate b = pinfo.getPos(ang.j);

                Area += MeasureAreaTriangle(a, b);
            }
        }

	}

    return Area;

}

double Simulation::MeasureAreaTriangle(Coordinate a, Coordinate b)
{
    double Area = 0.0;

    Coordinate aUnit = a.getUnitCoord();

    Coordinate normAB = {-(b.y-a.y), (b.x-a.x), 0.0};
    normAB = normAB.getUnitCoord();

    double dotProd = aUnit*normAB;

    int sign = 1;

    if(dotProd > 0.0)
    {
        sign = 1;
    }
    else if(dotProd < 0.0)
    {
        sign = -1;
    }
    else
    {

        sign = 0;
    }

    double nonSignedArea = fabs(0.5*(-b.x*a.y + a.x*b.y));

    double signedArea = sign * nonSignedArea;

    Area = signedArea;

    return Area;
}


void Simulation::run()
{
    std::vector<std::string> xyzArray;
    std::vector<std::string> bondArray;

    bool verbose = false;
    double t1, t2;
    t1 = omp_get_wtime();
    double initialSimTime = currSimTime;

    omp_set_num_threads(numThreads);

    std::vector<Coordinate> tempVec;
    std::vector<std::vector<Coordinate>> tempVecVec;

    for(int j = 0; j < pinfo.getNumParticles(); j++)
    {
        Coordinate c = {0.0, 0.0, 0.0};
        tempVec.push_back(c);
    }

    for(int k = 0; k < num_force_types; k++)
    {
        tempVecVec.push_back(tempVec);
    }

    for(int i = 0; i < omp_get_max_threads(); i++)
    {
        tempForces.push_back(tempVecVec);
        tempStresses.push_back(tempVec);
    }

    double tempCrosslinkLifetime = 1.0;
    int saveInterval = 70*snapshotStep; 

    std::ofstream outputXYZ("out-" + outputFileName + ".xyz");
    std::ofstream outputBND("out-" + outputFileName + ".bnd");
	
	Coordinate r_nuc_com_previous = r_nuc_com;
    double previousTime = 0.0;
    for(int mainLoopVar = 0; mainLoopVar < numSteps; mainLoopVar++)
    {
        if(verbose == true)
            cout << "before putting in grid" << endl;

        if(mainLoopVar % updateGridStep == 0)
        {
            putAllParticlesInGrid();
        }

        double leadingEdgeForce = calcForces();
        moveParticles();
        currSimTime += dt;

        if(verbose == true)
            cout << "before snapshot" << endl;

        if(mainLoopVar % snapshotStep == 0)
        {
            std::cout << fixed << currSimTime << " / " << (totalSimulationTime + initialSimTime) << endl;

            writeXYZFile(xyzArray);
            writeBondFile(bondArray);

            // Repolarization of the cell leading edge beads across the rest of the cell cortex beads.
            if(changeDirection)
			{
				if((currSimTime - previousTime) > 16000.0){
					std::mt19937& threadEngine = generators[omp_get_thread_num()];
					double rnd = (dis(threadEngine) - 0.5)*2.0; 

					shiftBeadTypes(int(rnd*numLeadingEdgeBeads)); // For the case of narrow repolarization in which L_Repol=L_Front; L_Front is given as the number of leading edge beads in these simulations.

					previousTime = currSimTime;
				}
				r_nuc_com_previous = r_nuc_com;
			}
        }

        if(mainLoopVar % saveInterval == 0)
        {
            outputXYZ << std::setprecision(16);
            outputBND << std::setprecision(16);

            for (const auto& line : xyzArray) {
                outputXYZ << line;
            }

            for (const auto& line : bondArray) {
                outputBND << line;
            }

            xyzArray.clear();
            bondArray.clear();
        }
    }
	
    outputXYZ.close();
    outputBND.close();


    t2 = omp_get_wtime();
    double seconds = t2 - t1;
    std::cout << "Time to run: " << seconds << endl;
    cout << seconds << endl;
}