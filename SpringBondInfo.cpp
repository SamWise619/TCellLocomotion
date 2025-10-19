// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#include "SpringBondInfo.h"

using namespace std;

SpringBondInfo::SpringBondInfo()
{
	k = 0.0;
	eqDist = 0.0;
}

SpringBondInfo::SpringBondInfo(double newK, double newEqDist)
{
	k = newK;
	eqDist = newEqDist;
}

double SpringBondInfo::calcForce(double dist)
{
	return -k * (dist - eqDist);
}

double SpringBondInfo::getEqDist()
{
    return eqDist;
}