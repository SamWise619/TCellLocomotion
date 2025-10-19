// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#ifndef SPRINGBONDINFO_H
#define SPRINGBONDINFO_H

#include "BondInfo.h"

class SpringBondInfo:public BondInfo
{
	private:
		double k;
		double eqDist;
	public:
		SpringBondInfo();
		SpringBondInfo(double newK, double newEqDist);
		double calcForce(double dist);
        double getEqDist();
};

#endif