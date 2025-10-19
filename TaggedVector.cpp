// Sami Alawadhi, David M. Rutkowski (2025) 
// Vavylonis Group
// Department of Physics, Lehigh University
// If you use any part of the scripts in this package, please cite:
// TBD

#include "TaggedVector.h"
#include <iostream>

using namespace std;

TaggedVector::TaggedVector()
{
	numCurrEntries = 0;
    maxHistoricalTag = 0;
}

int TaggedVector::add()
{
    int newTag;
    
    if(q_recycledTags.size() > 0)
    {
        newTag = q_recycledTags.front();
        q_recycledTags.pop();
    }
    else
    {
        // if recycling tags use numCurrEntires here
        newTag = maxHistoricalTag;
    }
    
    if(newTag < v_tag.size())
        v_tag[newTag] = numCurrEntries;
    else
        v_tag.push_back(numCurrEntries);
    
    if(numCurrEntries < v_rtag.size())
    {
        v_rtag[numCurrEntries] = newTag;
    }
    else
        v_rtag.push_back(newTag);
    
    numCurrEntries++;
    maxHistoricalTag++;
    
    if(maxHistoricalTag == 0)
    {
        throw std::runtime_error("taggedVector ran out of tags!");
    }
    
    return newTag;
}

void TaggedVector::remove(int oldTag)
{
    if(oldTag < v_tag.size())
    {
        int tempIndex = v_tag[oldTag];
        
        if(tempIndex > -1)
        {
            numCurrEntries--;
            
            // swappedTag is the last tag
            int swappedTag = v_rtag[numCurrEntries];
            v_tag[swappedTag] = tempIndex;
            v_rtag[tempIndex] = swappedTag;
            
            v_tag[oldTag] = -1;
            v_rtag[numCurrEntries] = -1;
        }
        else
        {
            throw std::runtime_error("Attempted to remove tag that didnt exist: " + std::to_string(oldTag));
        }
    }
}

int TaggedVector::getNumEntries()
{
    return numCurrEntries;
}

int TaggedVector::getTagAtIndex(int index)
{
    if(index < 0 || index >= numCurrEntries)
    {
        return -2;
    }
    
    return v_rtag[index];
}

int TaggedVector::getIndexOfTag(int tag)
{
    if(tag < 0 || tag >= v_tag.size())
    {
        return -2;
    }
    
    return v_tag[tag];
}

int TaggedVector::getCurrMaxTag()
{
    return maxHistoricalTag-1;
}

void TaggedVector::setTagAtPos(int pos, int newTag)
{
    if(pos < 0 || pos >= v_rtag.size())
        throw std::runtime_error("Cannot change tag at position " + std::to_string(pos) + " out of bounds of 0-" + std::to_string(v_rtag.size()));

    int oldTag = v_rtag[pos];
    v_rtag[pos] = newTag;
    
    if(newTag >= v_tag.size())
        v_tag.resize(newTag+1);    
    v_tag[newTag] = pos;
    
    if(maxHistoricalTag < newTag)
        maxHistoricalTag = newTag+1;
    else if(maxHistoricalTag-1 == oldTag)
        maxHistoricalTag--;
}

bool TaggedVector::checkTagExists(int tag)
{
    if(tag < v_tag.size())
    {
        int tempIndex = v_tag[tag];
        
        if(tempIndex > -1)
        {
            return true;
        }
    }
    
    return false;
}