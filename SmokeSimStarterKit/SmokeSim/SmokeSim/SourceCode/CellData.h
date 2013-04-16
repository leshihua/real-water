#pragma once

#pragma warning(disable: 4244 4267 4996)
#include <vector>
#include "vec.h"
#include "Constants.h"

class CellData
{
public:
	CellData(void);
	CellData(const CellData& orig);
	~CellData(void);
	CellData& operator=(const CellData& orig);
	std::vector<double>& data();
	virtual double& operator()(int i, int k);
	virtual const double operator()(int i, int k) const;
	//virtual double interpolate(const vec3& pt);
	virtual void initialize(double dftlValue);


protected:
	vec3 worldToSelf(const vec3& pt) const;
	double mDfltValue;
	vec2 mMax;
	std::vector<double> mData;

};

