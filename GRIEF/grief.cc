/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000-2008, Intel Corporation, all rights reserved.
// Copyright (C) 2009-2010, Willow Garage Inc., all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/
#include "grief.h"
#include "precomp.hpp"
#include <algorithm>
#include <vector>

#include <iostream>
#include <iomanip>
 
namespace cv
{
	namespace xfeatures2d
	{

		float evaluation(Eigen::MatrixXd individual){
			return -1;
		}

		Ptr<GriefDescriptorExtractor> GriefDescriptorExtractor::create(int bytes, bool use_orientation, EvalFunction evaluation, int N_pop, 
																	   float cr, float jr, float F, int mutation_algorithm, int crossover_algorithm)
		{
			return makePtr<GriefDescriptorExtractorImpl>(bytes, use_orientation, evaluation, N_pop, cr, jr, F, mutation_algorithm, crossover_algorithm);
		}

		std::vector<float> GriefDescriptorExtractorImpl::gbfit(){
			return bfit;
		}

		std::vector<float> GriefDescriptorExtractor::gbfit(){
			return std::vector<float>{0,0};
		}
		
		Eigen::MatrixXd GriefDescriptorExtractorImpl::get_best_indv(){
			return get_best_ind();
		}

		Eigen::MatrixXd GriefDescriptorExtractor::get_best_indv(){
			return Eigen::MatrixXd(1,1);
		}
#include <unistd.h>
		int GriefDescriptorExtractorImpl::load(std::string fileName) {
			std::string CURRENT_DIR = get_current_dir_name();
			using namespace std;
			ifstream file(CURRENT_DIR +"/../GRIEF_CUDA/" + fileName);

			std::string line;
			uint16_t i = 0, j = 0;
			bool successful=false;
			std::string cell;
			//std::cout << CURRENT_DIR +"../GRIEF_CUDA/" + fileName;
			while (std::getline(file, line)) {
				//std::cout << line;
				std::vector<int> v;
				istringstream is(line);
				while (std::getline(is, cell, ' ')) {
					//std::cout << std::stoi(cell);
					individual[i][j] = std::stoi(cell);
					j++;
				}
				successful=true;
				i++;
				j = 0;
			}

			file.close();
			return successful;
		}

		float GriefDescriptorExtractorImpl::get_b_fit(){
			return get_best_fit();
		}

		float GriefDescriptorExtractor::get_b_fit(){
			return 1;
		}

		inline int smoothedSum(const Mat& sum, const KeyPoint& pt, int y, int x)
		{
			static const int HALF_KERNEL = GriefDescriptorExtractorImpl::KERNEL_SIZE / 2;

			int img_y = (int)(pt.pt.y + 0.5) + y;
			int img_x = (int)(pt.pt.x + 0.5) + x;
			return   sum.at<int>(img_y + HALF_KERNEL + 1, img_x + HALF_KERNEL + 1)
				- sum.at<int>(img_y + HALF_KERNEL + 1, img_x - HALF_KERNEL)
				- sum.at<int>(img_y - HALF_KERNEL, img_x + HALF_KERNEL + 1)
				+ sum.at<int>(img_y - HALF_KERNEL, img_x - HALF_KERNEL);
		}

		static void pixelTests16(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, int individual[512][4])
		{
			Matx21f R;
			Mat sum = _sum.getMat(), descriptors = _descriptors.getMat();
			for (size_t i = 0; i < keypoints.size(); ++i)
			{
				uchar* desc = descriptors.ptr(static_cast<int>(i));
				const KeyPoint& pt = keypoints[i];
				if ( use_orientation )
				{
					float angle = pt.angle;
					angle *= (float)(CV_PI/180.f);
					R(0,0) = sin(angle);
					R(1,0) = cos(angle);
				}

#include "generated_16.i"
			}
		}

		static void pixelTests32(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, int individual[512][4])
		{
			Matx21f R;
			Mat sum = _sum.getMat(), descriptors = _descriptors.getMat();
			for (size_t i = 0; i < keypoints.size(); ++i)
			{
				uchar* desc = descriptors.ptr(static_cast<int>(i));
				const KeyPoint& pt = keypoints[i];
				if ( use_orientation )
				{
				float angle = pt.angle;
				angle *= (float)(CV_PI / 180.f);
				R(0,0) = sin(angle);
				R(1,0) = cos(angle);
				}

#include "generated_32.i"
			}
		}

		static void pixelTests64(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, int individual[512][4])
		{
			Matx21f R;
			Mat sum = _sum.getMat(), descriptors = _descriptors.getMat();
			for (size_t i = 0; i < keypoints.size(); ++i)
			{
				uchar* desc = descriptors.ptr(static_cast<int>(i));
				const KeyPoint& pt = keypoints[i];
				if ( use_orientation )
				{
					float angle = pt.angle;
					angle *= (float)(CV_PI/180.f);
					R(0,0) = sin(angle);
					R(1,0) = cos(angle);
				}

#include "generated_64.i"
			}
		}


		void GriefDescriptorExtractorImpl::evolve(uint ng){
			
			std::random_device rseed;
			std::mt19937 rng(rseed());
			std::uniform_real_distribution<float> distr(0,1);

			for(int g = 0; g < ng; g++){
				std::cout << "Gen " << g+1 << std::endl;

				auto start = std::chrono::high_resolution_clock::now();

				#if OPPOSITION_LEARNING
					if(g == 0 || distr(rng) <= jr){
						// std::cout << "Applying Opposition..." << std::endl;
						apply_opposition();
						continue;
					}
				#endif

				for(int i = 0; i < N_pop; i++){
					
					// std::cout << "[ Applying DE mutation and crossover operators ]" << std::endl;
					
					mutate(i);
					crossover(i);
					if(is_infeasible())
						repair(i);
					selection(i);
				}
				auto finish = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> elapsed = finish - start;
				std::cout << "Gen " << g+1 << ": Elapsed time: " << elapsed.count() << " ms." << std::endl;
				
			}
		}

		void GriefDescriptorExtractor::evolve(uint ng){
			
		}

		GriefDescriptorExtractorImpl::GriefDescriptorExtractorImpl( int bytes, bool use_orientation, EvalFunction evaluation, 
																	int N_pop, float cr, float jr, float F, int mutation_algorithm, int crossover_algorithm) :
			bytes_(bytes), test_fn_(NULL), 
			DE(N_pop, std::vector<int>{bytes*8, 4}, cr, jr, evaluation, F, MAXIMIZATION, std::vector<int>{-24, 24}, mutation_algorithm, crossover_algorithm)
		{
			this->N_pop = N_pop;
			this->jr = jr;

			load("test_pairs.brief");
			use_orientation_ = use_orientation;

			switch (bytes)
			{
				case 16:
					test_fn_ = pixelTests16;
					break;
				case 32:
					test_fn_ = pixelTests32;
					break;
				case 64:
					test_fn_ = pixelTests64;
					break;
				default:
					CV_Error(Error::StsBadArg, "bytes must be 16, 32, or 64");
			}
		}

		int GriefDescriptorExtractorImpl::descriptorSize() const
		{
			return bytes_;
		}

		int GriefDescriptorExtractorImpl::descriptorType() const
		{
			return CV_8UC1;
		}

		int GriefDescriptorExtractorImpl::defaultNorm() const
		{
			return NORM_HAMMING;
		}

		void GriefDescriptorExtractorImpl::read( const FileNode& fn)
		{
			int dSize = fn["descriptorSize"];
			switch (dSize)
			{
				case 16:
					test_fn_ = pixelTests16;
					break;
				case 32:
					test_fn_ = pixelTests32;
					break;
				case 64:
					test_fn_ = pixelTests64;
					break;
				default:
					CV_Error(Error::StsBadArg, "descriptorSize must be 16, 32, or 64");
			}
			bytes_ = dSize;
		}

		void GriefDescriptorExtractorImpl::write( FileStorage& fs) const
		{
			fs << "descriptorSize" << bytes_;
		}

		void GriefDescriptorExtractorImpl::getInd(){
			for(int i = 0; i < 512; i++){
				for(int j=0; j < 4; j++)
					std::cout << individual[i][j] << " ";
				std::cout << std::endl;
			}

		}

		void GriefDescriptorExtractor::getInd(){
			
		}

		void GriefDescriptorExtractorImpl::compute(InputArray image, std::vector<KeyPoint>& keypoints, OutputArray descriptors)
		{
			// Construct integral image for fast smoothing (box filter)
			Mat sum;

			Mat grayImage = image.getMat();
			if( image.type() != CV_8U ) cvtColor( image, grayImage, COLOR_BGR2GRAY );

			///TODO allow the user to pass in a precomputed integral image
			//if(image.type() == CV_32S)
			//  sum = image;
			//else

			integral( grayImage, sum, CV_32S);

			//Remove keypoints very close to the border
			KeyPointsFilter::runByImageBorder(keypoints, image.size(), PATCH_SIZE/2 + KERNEL_SIZE/2);

			descriptors.create((int)keypoints.size(), bytes_, CV_8U);
			descriptors.setTo(Scalar::all(0));
			test_fn_(sum, keypoints, descriptors, use_orientation_, individual);
		}

		void GriefDescriptorExtractorImpl::setInd(Eigen::MatrixXd new_individual){
			// load("test_pairs.brief");
			// std::cout << individual[0][0] << individual[0][1] << individual[0][2] << individual[0][3]  << std::endl;
			for(int i = 0; i < bytes_*8; i++){
				for(int j = 0; j < 4; j++){
					individual[i][j] = new_individual(i,j);
				}
			}
			//std::cout << individual[0][0] << individual[0][1] << individual[0][2] << individual[0][3]  << std::endl;
		}

		void GriefDescriptorExtractor::setInd(Eigen::MatrixXd new_individual){
		
		}
	} //namespace xfeatures2d
} // namespace cv