// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#ifndef BONDINFO_H
#define BONDINFO_H

class BondInfo
{
	public:
		virtual double calcForce(double dist)
		{ return 0.0; }
        virtual double getEqDist()
        { return 0.0; }
};

#endif