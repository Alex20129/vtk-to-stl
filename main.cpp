#include <iostream>
#include <string>

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkResampleToImage.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageConstantPad.h>
#include <vtkMarchingCubes.h>
#include <vtkSTLWriter.h>

enum symmetryAxis
{
	XAxis,
	YAxis,
	ZAxis,
};

void MakeSymmetrical(vtkImageData *imageData, symmetryAxis axis)
{
	vtkDataArray *scalars = imageData->GetPointData()->GetScalars("element_states_averaged_at_nodes");
	int dims[3];
	imageData->GetDimensions(dims);
	double bounds[6];
	imageData->GetBounds(bounds);
	int symAxis = axis;
	int iMax = dims[0];
	int jMax = dims[1];
	int kMax = dims[2];
	for (int k = 0; k < kMax; ++k)
	{
		for (int j = 0; j < jMax; ++j)
		{
			for (int i = 0; i < iMax; ++i)
			{
				double coord[3];
				int ijk[3] = {i, j, k};
				vtkIdType pointId = imageData->ComputePointId(ijk);
				imageData->GetPoint(pointId, coord);
				if (coord[symAxis] < 0)
				{
					continue;
				}
				double symCoord[3] = {coord[0], coord[1], coord[2]};
				symCoord[symAxis] = -coord[symAxis];
				int symI = i, symJ = j, symK = k;
				double spacing[3];
				imageData->GetSpacing(spacing);
				double origin[3];
				imageData->GetOrigin(origin);
				if (symAxis == XAxis)
				{
					symI = round((symCoord[XAxis] - origin[0]) / spacing[0]);
				}
				else if (symAxis == YAxis)
				{
					symJ = round((symCoord[YAxis] - origin[1]) / spacing[1]);
				}
				else // symAxis == ZAxis
				{
					symK = round((symCoord[ZAxis] - origin[2]) / spacing[2]);
				}
				if (symI < 0 || symI >= iMax || symJ < 0 || symJ >= jMax || symK < 0 || symK >= kMax)
				{
					continue;
				}
				int symijk[3] = {symI, symJ, symK};
				vtkIdType symPointId = imageData->ComputePointId(symijk);
				double value = scalars->GetComponent(pointId, 0);
				double symValue = scalars->GetComponent(symPointId, 0);
				double avgValue = (value + symValue) / 2.0;
				scalars->SetComponent(pointId, 0, avgValue);
				scalars->SetComponent(symPointId, 0, avgValue);
			}
		}
	}
	scalars->Modified();
}

int main(int argc, char *argv[])
{
	double threshold=0.5, resolution=128;
	double gauss_radius=0.0, gauss_deviation=0.0;
	bool XsymmetryFlag=false, YsymmetryFlag=false, ZsymmetryFlag=false;

	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0] << " <input.vtk> <output.stl> [threshold=0.5] [resolution=128] [gauss_radius=0.0] [gauss_deviation=0.0] [symmetry]\n";
		return 1;
	}

	std::string inputFilename=argv[1];
	std::string outputFilename=argv[2];
	if (argc > 3)
	{
		threshold=std::stod(argv[3]);
	}
	if (argc > 4)
	{
		resolution=std::stod(argv[4]);
	}
	if (argc > 5)
	{
		gauss_radius=std::stod(argv[5]);
	}
	if (argc > 6)
	{
		gauss_deviation=std::stod(argv[6]);
	}
	if (argc > 7)
	{
		std::string symmetryArg = std::string(argv[7]);
		for (char c : symmetryArg)
		{
			if (c == 'x' || c == 'X')
			{
				XsymmetryFlag = true;
			}
			if (c == 'y' || c == 'Y')
			{
				YsymmetryFlag = true;
			}
			if (c == 'z' || c == 'Z')
			{
				ZsymmetryFlag = true;
			}
		}
	}

	if(gauss_deviation<=0.0)
	{
		gauss_deviation=gauss_radius=0.0;
	}
	if(gauss_radius<=0.0)
	{
		gauss_deviation=gauss_radius=0.0;
	}

	vtkUnstructuredGridReader *vtkGridReader=vtkUnstructuredGridReader::New();
	vtkGridReader->SetFileName(inputFilename.c_str());

	vtkGridReader->SetScalarsName("element_states_averaged_at_nodes");
	vtkGridReader->Update();

	if (!vtkGridReader->GetOutput())
	{
		std::cerr << "Error: Could not read VTK file: " << inputFilename << std::endl;
		return 1;
	}

	vtkUnstructuredGrid *grid=vtkGridReader->GetOutput();
	std::cout << *grid;

	double gridBounds[6];
	grid->GetBounds(gridBounds);

	double xSize=gridBounds[1] - gridBounds[0];
	double ySize=gridBounds[3] - gridBounds[2];
	double zSize=gridBounds[5] - gridBounds[4];
	std::cout << "Model bounds:\n"
		<< "  X: [" << gridBounds[0] << ", " << gridBounds[1] << "] (size: " << xSize << ")\n"
		<< "  Y: [" << gridBounds[2] << ", " << gridBounds[3] << "] (size: " << ySize << ")\n"
		<< "  Z: [" << gridBounds[4] << ", " << gridBounds[5] << "] (size: " << zSize << ")\n";

	double voxelSize=std::min({xSize, ySize, zSize});
	voxelSize /= resolution;

	int nx=std::max(2, (int)std::round(xSize / voxelSize));
	int ny=std::max(2, (int)std::round(ySize / voxelSize));
	int nz=std::max(2, (int)std::round(zSize / voxelSize));
	std::cout << "Sampling dimensions: " << nx << "x" << ny << "x" << nz << std::endl;

	vtkPointData *pointData=grid->GetPointData();
	if (!pointData->GetScalars("element_states_averaged_at_nodes"))
	{
		std::cerr << "Error: Scalar field 'element_states_averaged_at_nodes' not found in VTK file.\n";
		std::cerr << "Available scalar fields in POINT_DATA:\n";
		for (int i=0; i < pointData->GetNumberOfArrays(); ++i)
		{
			const char *arrayName=pointData->GetArrayName(i);
			if (arrayName)
			{
				std::cerr << " - " << arrayName << std::endl;
			}
		}
		return 1;
	}

	vtkImageData *imageData;

	vtkResampleToImage *resampler=vtkResampleToImage::New();
	resampler->SetInputDataObject(grid);
	resampler->SetSamplingDimensions(nx, ny, nz);
	resampler->SetUseInputBounds(true);
	resampler->Update();
	imageData=resampler->GetOutput();

	imageData->GetPointData()->SetActiveScalars("element_states_averaged_at_nodes");

	if(XsymmetryFlag)
	{
		MakeSymmetrical(imageData, XAxis);
	}
	if(YsymmetryFlag)
	{
		MakeSymmetrical(imageData, YAxis);
	}
	if(ZsymmetryFlag)
	{
		MakeSymmetrical(imageData, ZAxis);
	}

	double padding=1.0+gauss_radius+gauss_deviation;
	vtkImageConstantPad *padder=vtkImageConstantPad::New();
	padder->SetInputData(imageData);
	padder->SetConstant(0.0);
	padder->SetOutputWholeExtent(-padding, nx+padding, -padding, ny+padding, -padding, nz+padding);
	padder->Update();
	imageData=padder->GetOutput();

	if(gauss_radius>0.0 && gauss_deviation>0.0)
	{
		vtkImageGaussianSmooth *gaussianSmooth=vtkImageGaussianSmooth::New();
		gaussianSmooth->SetInputData(imageData);
		gaussianSmooth->SetRadiusFactor(gauss_radius);
		gaussianSmooth->SetStandardDeviation(gauss_deviation, gauss_deviation, gauss_deviation);
		gaussianSmooth->Update();
		imageData=gaussianSmooth->GetOutput();
	}

	std::cout << *imageData;

	vtkMarchingCubes *marchingCubes=vtkMarchingCubes::New();
	marchingCubes->SetInputData(imageData);
	marchingCubes->SetNumberOfContours(1);
	marchingCubes->SetValue(0, threshold);
	marchingCubes->SetComputeNormals(true);
	marchingCubes->SetComputeGradients(true);
	marchingCubes->SetComputeScalars(true);
	marchingCubes->Update();

	vtkSTLWriter *stlWriter=vtkSTLWriter::New();
	stlWriter->SetInputData(marchingCubes->GetOutput());
	stlWriter->SetFileName(outputFilename.c_str());
	stlWriter->SetFileTypeToASCII();
	stlWriter->Write();

	if (stlWriter->GetErrorCode())
	{
		std::cerr << "Failed to generate STL file: " << outputFilename << std::endl;
		return 1;
	}
	else
	{
		std::cout << "Successfully generated STL file: " << outputFilename << std::endl;
	}

	return(0);
}
