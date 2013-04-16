#include "CellData.h"

//Constructor
CellData::CellData() : 
	mDfltValue(0.0), mMax(0.0, 0.0) {
	//Nothing
}

//Copy Constructor
CellData::CellData(const CellData& orig) :
mDfltValue(orig.mDfltValue) {
	mData = orig.mData;
	mMax = orig.mMax;
}

//Destructor
CellData::~CellData() {
	//Nothing
}

//Returns the height data in a vector
std::vector<double>& CellData::data() {
	return mData;
}

//Setter  constructor
CellData& CellData::operator=(const CellData& orig) {
	if (this == &orig) {
		return *this;
	}
	
	mDfltValue = orig.mDfltValue;
	mData = orig.mData;
	mMax = orig.mMax;
	return *this;
}
//Initialize values in the terrain heights
void CellData::initialize(double dfltValue) {
	mDfltValue = dfltValue;
	mMax[0] = theCellSize*theDim[0];  //x-component
    mMax[1] = theCellSize*theDim[2];  //z-component
	mData.resize(theDim[0]*theDim[2], false);
    std::fill(mData.begin(), mData.end(), mDfltValue);
}


double& CellData::operator()(int i, int k) {
	static double dflt = 0;
	dflt = mDfltValue;
	
	//Default - OUT OF BOUNDS Case
	if (i< 0 || k<0 || 
       i > theDim[0]-1 || 
       k > theDim[2]-1) return dflt;

	int x = i;
	int z = k * theDim[2];

	return mData[x+z];
}

const double CellData::operator()(int i, int k) const {
	static double dflt = 0;
	dflt = mDfltValue;
	
	//Default - OUT OF BOUNDS Case
	if (i < 0 || k < 0 || 
       i > theDim[0]-1 || 
       k > theDim[2]-1) return dflt;

	int x = i;
	int z = k * theDim[2];

	return mData[x+z];

}