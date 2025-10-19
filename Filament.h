// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#ifndef FILAMENT_H
#define FILAMENT_H

#include <vector>
#include <queue>
#include <iostream>
#include "Coordinate.h"
#include "ParticleInfo.h"

class Filament
{
	private:
		std::deque <int> p_tags;
	public:
		int getNumParticles();
		int getTagAtIndex(int pos);
		int addParticleFront(int newTag);
		int addParticleBack(int newTag);
        int getIndexOfTag(int tempTag);

};

#endif