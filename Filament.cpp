// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#include "Filament.h"

using namespace std;

int Filament::getNumParticles()
{
	return p_tags.size();
}

int Filament::addParticleFront(int newTag)
{
	p_tags.push_front(newTag);
		
	if(p_tags.size() > 1)
	{
		int oldFront = p_tags[1];
		return oldFront;
	}
	else
	{
		return -1;
	}
}

int Filament::addParticleBack(int newTag)
{
	p_tags.push_back(newTag);
	
	if(p_tags.size() > 1)
	{
		int oldBack = p_tags[p_tags.size()-2];
		return oldBack;
	}
	else
	{
		return -1;
	}	
}

int Filament::getTagAtIndex(int index)
{
	if(index < p_tags.size())
	{
		return p_tags[index];
	}
	else
	{
		return -1;
	}
}

int Filament::getIndexOfTag(int tempTag)
{
    for(int i = 0; i < p_tags.size(); i++)
    {
        if(p_tags[i] == tempTag)
            return i;
    }
    
   return -1;
}