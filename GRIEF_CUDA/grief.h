#include <sys/time.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d.hpp>
#include <stdio.h>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include <opencv2/core/core.hpp>
#include <opencv2/flann/miniflann.hpp>
#include "opencv2/xfeatures2d.hpp"
#include <Eigen/Dense>
#include "../DE/DE.h"
//#include "matplotlibcpp.h"
namespace cv{

namespace xfeatures2d
{

typedef int arr64[64];
typedef int arr4[4];
typedef int arr2[2];
float evaluation(Eigen::MatrixXd individual);

class CV_EXPORTS_W GriefDescriptorExtractor : public Feature2D
{
public:
	CV_WRAP static Ptr<GriefDescriptorExtractor> create( int bytes = 32, bool use_orientation = false, EvalFunction evaluation = evaluation, 
				int N_pop = 0, float cr = 0.5, float jr = 0.3, float F = 0.5, int mutation_algorithm=RAND_1, int crossover_algorithm=BIN);
	CV_WRAP virtual void getInd( );
	CV_WRAP virtual void setInd(Eigen::MatrixXd new_individual);
	CV_WRAP virtual void evolve(uint ng);
	CV_WRAP virtual float get_b_fit();
	CV_WRAP virtual std::vector<float> gbfit();
	CV_WRAP virtual void compute(InputArray image, std::vector<KeyPoint>& keypoints, cuda::GpuMat& descriptors);
	//CV_WRAP virtual void plot_convergence();
};

/*
 * Grief Descriptor
 */
class GriefDescriptorExtractorImpl : public GriefDescriptorExtractor, public DE::DE
{
public:
	enum { PATCH_SIZE = 48, KERNEL_SIZE = 9 };

	// bytes is a length of descriptor in bytes. It can be equal 16, 32 or 64 bytes.
	GriefDescriptorExtractorImpl( int bytes = 32, bool use_orientation = false, EvalFunction evaluation = evaluation, 
											  int N_pop = 0, float cr = 0.9, float jr = 0.2, float F = 0.8, int mutation_algorithm=RAND_1, int crossover_algorithm=BIN);
											  
	int load(int mat[512][4], std::string fileName);
	virtual void read( const FileNode& ) CV_OVERRIDE;
	virtual void write( FileStorage& ) const CV_OVERRIDE;
	
	virtual int descriptorSize() const CV_OVERRIDE;
	virtual int descriptorType() const CV_OVERRIDE;
	virtual int defaultNorm() const CV_OVERRIDE;

	virtual void compute(InputArray image, std::vector<KeyPoint>& keypoints, cuda::GpuMat& descriptors) CV_OVERRIDE;
	virtual void getInd() CV_OVERRIDE;
	virtual void evolve(uint ng) CV_OVERRIDE;
	virtual void setInd(Eigen::MatrixXd new_individual) CV_OVERRIDE;
	virtual float get_b_fit() CV_OVERRIDE;
	virtual std::vector<float> gbfit() CV_OVERRIDE;
	//virtual void plot_convergence() CV_OVERRIDE;
	void pixelTests64(InputArray sum, const std::vector<KeyPoint>& keypoints,cuda::GpuMat& descriptors, bool use_orientation, int individual[512][4]);
	
protected:
	std::vector<float> bfit;
	int N_pop;
	//typedef void(*PixelTestFn)(InputArray, const std::vector<KeyPoint>&, OutputArray, bool use_orientation, int individual[512][4]);
	int individual[512][4];
	int bytes_;
	bool use_orientation_;
	
	//PixelTestFn test_fn_;
};
}
}
