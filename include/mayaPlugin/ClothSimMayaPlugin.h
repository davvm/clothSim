#ifndef ClothSimMayaPlugin_H
#define ClothSimMayaPlugin_H

#include <memory>

#include <maya/MPxNode.h> 
#include <maya/MVectorArray.h> 

#include "simLib\ClothMesh.h"

class ClothSimMayaPlugin : public MPxNode
{
public:

	ClothSimMayaPlugin();
	virtual ~ClothSimMayaPlugin();

	virtual MStatus compute(const MPlug& plug, MDataBlock& data);

	static void* creator();
	static MStatus initialize();

	static MTypeId id;

	static MObject g_aTime;
	static MObject g_aOutputMesh;

private:

	MObject createMesh(const MTime& time, MObject& outData, MStatus& stat);

	std::auto_ptr< ClothMesh<float> > m_simMesh;
	double m_prevTime;

};

#endif // ClothSimMayaPlugin_H