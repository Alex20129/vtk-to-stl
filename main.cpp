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

int main(int argc, char *argv[])
{
	double threshold=0.5, resolution=128;
	double gauss_radius=0.0, gauss_deviation=0.0;

	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0] << " <input.vtk> <output.stl> [threshold=0.5] [resolution=128] [gauss_radius=0.0] [gauss_deviation=0.0]\n";
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
