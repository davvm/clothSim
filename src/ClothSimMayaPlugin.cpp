#include "ClothSimMayaPlugin.h"

#include "BasicCGSolver.h"
#include "DirectSolver.h"

#include <maya/MFnVectorArrayData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnIntArrayData.h>

#include <maya/MGlobal.h>
#include <maya/MTime.h>
#include <maya/MIOStream.h>
#include <string.h>
#include <math.h>

#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnPlugin.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnMesh.h>

#include <maya/MFloatPointArray.h>
#include <maya/MString.h> 
#include <maya/MTypeId.h> 
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

MTypeId ClothSimMayaPlugin::id(0x8f17c);

MObject ClothSimMayaPlugin::g_aTime;
MObject ClothSimMayaPlugin::g_aOutputMesh;


#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
		}

ClothSimMayaPlugin::ClothSimMayaPlugin() : m_prevTime( 1.0/24 )
{
}

ClothSimMayaPlugin::~ClothSimMayaPlugin()
{
}


MStatus ClothSimMayaPlugin::initialize()
{
	MFnUnitAttribute unitAttr;
	MFnTypedAttribute typedAttr;

	MStatus returnStatus;

	g_aTime = unitAttr.create("time", "tm",	MFnUnitAttribute::kTime, 0.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating time attribute\n");

	g_aOutputMesh = typedAttr.create("outputMesh", "out", MFnData::kMesh, &returnStatus);
	McheckErr(returnStatus, "ERROR creating output attribute\n");
	typedAttr.setStorable(false);

	returnStatus = addAttribute(g_aTime);
	McheckErr(returnStatus, "ERROR adding time attribute\n");

	returnStatus = addAttribute(g_aOutputMesh);
	McheckErr(returnStatus, "ERROR adding outputMesh attribute\n");

	returnStatus = attributeAffects(g_aTime, g_aOutputMesh);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	return MS::kSuccess;
}

MStatus ClothSimMayaPlugin::compute(const MPlug& plug, MDataBlock& data)
{
	MStatus returnStatus;

	if (plug == g_aOutputMesh)
	{
		MDataHandle timeData = data.inputValue(g_aTime, &returnStatus);
		McheckErr(returnStatus, "Error getting time data handle\n");
		MTime time = timeData.asTime();

		/* Get output object */

		MDataHandle outputHandle = data.outputValue(g_aOutputMesh, &returnStatus);
		McheckErr(returnStatus, "ERROR getting polygon data handle\n");

		MFnMeshData dataCreator;
		MObject newOutputData = dataCreator.create(&returnStatus);
		McheckErr(returnStatus, "ERROR creating outputData");

		createMesh(time, newOutputData, returnStatus);
		if (!returnStatus)
		{
			std::cerr << "ERROR creating new Cube: " << returnStatus.errorString() << std::endl;
			return returnStatus;
		}

		outputHandle.set(newOutputData);
		data.setClean(plug);
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}

void* ClothSimMayaPlugin::creator()
{
	return new ClothSimMayaPlugin;
}

MObject ClothSimMayaPlugin::createMesh(const MTime& time,
	MObject& outData,
	MStatus& stat)

{
	double t = time.as(MTime::kSeconds);
	if (t <= 2.0 / 24)
	{
		m_simMesh.reset(0);
	}

	if (!m_simMesh.get())
	{
		int nx = 20;
		int ny = 20;

		Eigen::VectorXf v((nx+1) * (ny+1) * 3);
		Eigen::VectorXf x((nx + 1) * (ny + 1) * 3);
		Eigen::VectorXf uv((nx + 1) * (ny + 1) * 2);

		Eigen::Vector2f xzCoords(0, 0);
		for (int i = 0; i <= nx; ++i)
		{
			for (int j = 0; j <= ny; ++j)
			{
				int base = i + (nx+1) * j;
				uv[2 * base + 0] = (float)i / nx;
				uv[2 * base + 1] = (float)j / ny;

				x[3 * base + 0] = xzCoords[0];
				x[3 * base + 1] = (float)j / ny;
				x[3 * base + 2] = xzCoords[1];

				v[3 * base + 0] = v[3 * base + 1] = v[3 * base + 2] = 0;
			}
			float theta(i * 3.5f / nx);
			xzCoords += Eigen::Vector2f(cos(theta) / nx, sin(theta) / nx);
		}

		std::vector<int> triangleInds;
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < ny; ++j)
			{
				int base = i + (nx + 1) * j;

				triangleInds.push_back(base + 0);
				triangleInds.push_back(base + 1);
				triangleInds.push_back(base + (nx + 1));

				triangleInds.push_back(base + 1);
				triangleInds.push_back(base + (nx + 2));
				triangleInds.push_back(base + (nx + 1));

			}
		}

		m_simMesh.reset( new ClothMesh<float>(
			x, v, uv, triangleInds,
			0.01f, 10000.0f, 10000.0f,
			0.001f, 100.0f, 100.0f,
			1.0f
			));
	}

	if (t >= m_prevTime)
	{
		//BasicCGSolver<float> solver;
		DirectSolver<float> solver;
		try
		{
			m_simMesh->advance(float(t - m_prevTime), solver);
		}
		catch (const std::exception &e)
		{
			std::cerr << e.what() << std::endl;
			stat = MStatus::kFailure;
			return MObject();
		}
		catch (...)
		{
			std::cerr << "unknown exception" << std::endl;
			stat = MStatus::kFailure;
			return MObject();
		}
	}

	m_prevTime = t;

	MFloatPointArray points;
	for (int i = 0; i < m_simMesh->x().size(); i += 3)
	{
		MFloatPoint p(m_simMesh->x()[i], m_simMesh->x()[i + 1], m_simMesh->x()[i + 2]);
		points.append(p);
	}

	MFnMesh meshFS;
	MIntArray faceCounts((int)m_simMesh->triangleIndices().size()/3, 3);
	MIntArray faceConnects;
	for (unsigned i = 0; i < m_simMesh->triangleIndices().size(); ++i)
	{
		faceConnects.append(m_simMesh->triangleIndices()[i]);
	}

	MObject newMesh = meshFS.create((int)m_simMesh->x().size() / 3, (int)m_simMesh->triangleIndices().size() / 3, points, faceCounts, faceConnects, outData, &stat);

	return newMesh;
}

MStatus initializePlugin(MObject obj)
{
	MStatus   status;
	MFnPlugin plugin(obj, PLUGIN_COMPANY, "6.0", "Any");

	status = plugin.registerNode("ClothSimMayaPlugin", ClothSimMayaPlugin::id, ClothSimMayaPlugin::creator, ClothSimMayaPlugin::initialize);
	if (!status)
	{
		status.perror("registerNode");
		return(status);
	}

	return(status);
}

MStatus uninitializePlugin(MObject obj)
{
	MStatus   status;
	MFnPlugin plugin(obj);

	status = plugin.deregisterNode(ClothSimMayaPlugin::id);
	if (!status)
	{
		status.perror("deregisterNode");
		return(status);
	}

	return(status);
}
