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


namespace cv{
	namespace xfeatures2d
	{
		typedef int arr64[64];
		typedef int arr4[4];
		typedef int arr2[2];

		// #if CURRENT_TO_RAND||RAND_TO_BEST_MOD
		// std::vector<double> evaluation(Eigen::MatrixXd individual);
		// #else
		// float evaluation(Eigen::MatrixXd individual);
		// #endif

		#if ( CURRENT_TO_RAND || RAND_TO_BEST_MOD ) && BIN_CROSS_GENE
			void evaluation(Eigen::MatrixXd individual, std::vector<double> &fit, std::vector<float> &gene_fit_vec);
		#elif CURRENT_TO_RAND || RAND_TO_BEST_MOD
			std::vector<double> evaluation(Eigen::MatrixXd individual);
		#elif BIN_CROSS_GENE
			void evaluation(Eigen::MatrixXd individual, float &fit, std::vector<float> &gene_fit_vec);
		#else
			float evaluation(Eigen::MatrixXd individual);
		#endif


		class CV_EXPORTS_W GriefDescriptorExtractor : public Feature2D
		{
			public:
				CV_WRAP static Ptr<GriefDescriptorExtractor> create( int bytes = 32, bool use_orientation = false, EvalFunction evaluation = evaluation, 
							int N_pop = 0, int K=10, float cr = 0.8, float jr = 0.3, float F = 0.8, int mutation_algorithm=RAND_TO_BEST, int crossover_algorithm=BIN_G);
				CV_WRAP virtual void getInd( );
				CV_WRAP virtual void setInd(Eigen::MatrixXd new_individual);
				CV_WRAP virtual void evolve(uint ng);
				CV_WRAP virtual float get_b_fit();
				CV_WRAP virtual std::vector<float> get_change_percentage(uint ng);
				CV_WRAP virtual std::vector<float> gbfit();
				CV_WRAP virtual Eigen::MatrixXd get_best_indv();
				CV_WRAP virtual void compute(InputArray image, std::vector<KeyPoint>& keypoints, cuda::GpuMat& descriptors);
				//CV_WRAP virtual void plot_convergence();
		};

		/*
		* Grief Descriptor
		*/
		class GriefDescriptorExtractorImpl : public GriefDescriptorExtractor, public DE::DE
		{
			public:
				enum { PATCH_SIZE = 100, KERNEL_SIZE = 9 };

				// bytes is a length of descriptor in bytes. It can be equal 16, 32 or 64 bytes.
				GriefDescriptorExtractorImpl( int bytes = 32, bool use_orientation = false, EvalFunction evaluation = evaluation, 
														int N_pop = 0, int K=10, float cr = 0.8, float jr = 0.3, float F = 0.8, int mutation_algorithm=RAND_TO_BEST, int crossover_algorithm=BIN_G);
														
				int load(std::string fileName);
				virtual void read( const FileNode& ) CV_OVERRIDE;
				virtual void write( FileStorage& ) const CV_OVERRIDE;
				
				virtual int descriptorSize() const CV_OVERRIDE;
				virtual int descriptorType() const CV_OVERRIDE;
				virtual int defaultNorm() const CV_OVERRIDE;
				virtual std::vector<float> get_change_percentage(uint ng) CV_OVERRIDE;
				virtual void compute(InputArray image, std::vector<KeyPoint>& keypoints, cuda::GpuMat& descriptors) CV_OVERRIDE;
				virtual void getInd() CV_OVERRIDE;
				virtual void evolve(uint ng) CV_OVERRIDE;
				virtual void setInd(Eigen::MatrixXd new_individual) CV_OVERRIDE;
				virtual float get_b_fit() CV_OVERRIDE;
				virtual std::vector<float> gbfit() CV_OVERRIDE;
				//virtual void plot_convergence() CV_OVERRIDE;
				void pixelTests64(InputArray sum, const std::vector<KeyPoint>& keypoints,cuda::GpuMat& descriptors, bool use_orientation);

				virtual Eigen::MatrixXd get_best_indv() CV_OVERRIDE;
				
			protected:
				std::vector<float> bfit;
				int N_pop;
				std::vector<float> change_percentage;
				//typedef void(*PixelTestFn)(InputArray, const std::vector<KeyPoint>&, OutputArray, bool use_orientation, int individual[512][4]);
				int individual[512][4];
				int bytes_;
				bool use_orientation_;
				int exp = 1;
				//PixelTestFn test_fn_;
		};
	}
}
