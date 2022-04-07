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

__device__ int smoothedSum(cuda::PtrStepSz<int> sum, int x, int y, int _x, int _y)
{
	static const int HALF_KERNEL = GriefDescriptorExtractorImpl::KERNEL_SIZE / 2;
	//printf("ajsdihji%d\n", sum(1,1));
	//printf("passou aqui%d\n", x);
	int img_y = (int)(y + 0.5) + _y;
	int img_x = (int)(x + 0.5) + _x;
	return   sum(img_y + HALF_KERNEL + 1, img_x + HALF_KERNEL + 1)
		   - sum(img_y + HALF_KERNEL + 1, img_x - HALF_KERNEL)
		   - sum(img_y - HALF_KERNEL, img_x + HALF_KERNEL + 1)
		   + sum(img_y - HALF_KERNEL, img_x - HALF_KERNEL);
	//printf("%d e %d\n", individual[0][0], individual[0][0]);
	//printf("passou aqui%d\n", individual[0][0]);
}

__global__ void compare_results(uchar* desc, arr2 * result){
	//printf("%d\n", result[0][0]);
	desc[blockIdx.x] += (result[threadIdx.x][0] <  result[threadIdx.x][1])<< (7 - threadIdx.x);
	//printf("%d\n", desc[blockIdx.x]);
}

static void pixelTests16(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, int individual[512][4])
{
	Matx21f R;
	Mat sum = _sum.getMat(), descriptors = _descriptors.getMat();
	int result;
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

//#include "generated_16.i"
	}
}

static void pixelTests32(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, int individual[512][4])
{
	Matx21f R;
	
	Mat sum = _sum.getMat(), descriptors = _descriptors.getMat();
	exit(-1);
	int result;
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

//#include "generated_32.i"
	}
}

__global__ static void pixelTests64_kernel(cuda::PtrStepSz<int> sum, float* x, float* y,cuda::PtrStepSz<uchar>  descriptors, bool* use_orientation, arr4* individual){
	//Matx21f R;
	//arr2* result_child = nullptr;
	//arr4* individual_child = nullptr;
	//cudaMalloc(&result_child, (sizeof(int) *512*2));
	//cudaMalloc(&result_child2, (sizeof(int) *512*2));
	//cudaMemcpyAsync(result_child, result, (sizeof(int)*512*2), cudaMemcpyDeviceToDevice);

	//cudaMalloc(&individual_child, (sizeof(int) *512*4));
	//cudaMemcpyAsync(individual_child, individual, (sizeof(int)*512*4), cudaMemcpyDeviceToDevice);
	
	
	//cudaMalloc(&desc_child, (sizeof(uchar) *64));
	//printf("%d\n", individual[0][0]);

	//KeyPoint& pt = ;
	//if ( use_orientation )
	//{
	//  //float angle = pt.angle;
	//  //angle *= (float)(CV_PI/180.f);
	//  //R(0,0) = sin(angle);
	//  //R(1,0) = cos(angle);
	//}
	//printf("passou aqui%d\n", sum(1,1));
	//desc[0] = 5;
	//smoothedSum<<<512,2>>>(sum, x[blockIdx.x], y[blockIdx.x], individual_child, result_child);
	//cudaDeviceSynchronize();
	//memcpy(result, result_child, (sizeof(int)*512*2));
	//memcpy(result_child2, result, (sizeof(int)*512*2));
	//printf("%d\n", result[0][0]);
	//printf("opa%d", result[0][0]);
	//compare_results<<<64,8>>>(desc_child, result_child);
	//memcpy(desc, desc_child, (sizeof(int)*512*2));
	//cudaDeviceSynchronize();
	//desc[0] = 1;
	//printf("%d\n", desc[0]);
	descriptors(blockIdx.x,threadIdx.x) = ((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x][0], individual[8*threadIdx.x][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x][2], individual[8*threadIdx.x][3])) << 7)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 1][0], individual[8*threadIdx.x + 1][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 1][2], individual[8*threadIdx.x + 1][3])) << 6)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 2][0], individual[8*threadIdx.x + 2][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 2][2], individual[8*threadIdx.x + 2][3])) << 5)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 3][0], individual[8*threadIdx.x + 3][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 3][2], individual[8*threadIdx.x + 3][3])) << 4)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 4][0], individual[8*threadIdx.x + 4][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 4][2], individual[8*threadIdx.x + 4][3])) << 3)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 5][0], individual[8*threadIdx.x + 5][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 5][2], individual[8*threadIdx.x + 5][3])) << 2)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 6][0], individual[8*threadIdx.x + 6][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 6][2], individual[8*threadIdx.x + 6][3])) << 1)
	+((smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 7][0], individual[8*threadIdx.x + 7][1]) < smoothedSum(sum, x[blockIdx.x], y[blockIdx.x], individual[8*threadIdx.x + 7][2], individual[8*threadIdx.x + 7][3])));// << 0
//#include "generated_64.i"
}

Eigen::MatrixXd GriefDescriptorExtractorImpl::get_best_indv(){
	return get_best_ind();
}

Eigen::MatrixXd GriefDescriptorExtractor::get_best_indv(){
	return Eigen::MatrixXd(1,1);
}

void GriefDescriptorExtractorImpl::pixelTests64(InputArray sum, const std::vector<KeyPoint>& keypoints,cuda::GpuMat&  descriptors, bool use_orientation)
{
	//auto start = std::chrono::high_resolution_clock::now();
	cuda::GpuMat _sum;
	_sum.upload(sum.getMat());
	bool* _use_orientation;
	
	std::vector<KeyPoint> keypoints_aux  = keypoints;
	KeyPoint *arr_keypoints = &keypoints_aux[0], *gpu_keypoints;

	float*x, *y;
	float*_x, *_y;

	x = (float*)malloc(sizeof(float)*keypoints.size());
	y = (float*)malloc(sizeof(float)*keypoints.size());
	
	for(int i = 0; i < keypoints.size(); i++){
		x[i] = keypoints[i].pt.x;
		y[i] = keypoints[i].pt.y;
	}
	
	//printf("passou aqui%f\n", x[2]);
	//arr_keypoints = (KeyPoint*)*malloc(sizeof(KeyPoint));
	arr4* gpu_individual;
	
	//cuda::GpuMat _descriptors(keypoints.size(), 64, CV_8UC1);
	
	//cudaMalloc((void **)&_descriptors, sizeof(cv::KeyPoint));
	cudaMalloc((void **)&gpu_individual,sizeof(int)*512*4);
	
	cudaMalloc((void **)&_x,sizeof(float)*keypoints.size());
	cudaMalloc((void **)&_y,sizeof(float)*keypoints.size());
	//cudaMallocManaged(&_sum,sizeof(Mat)*keypoints.size());

	cudaMemcpy(_x, x, sizeof(float)*keypoints.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(_y, y, sizeof(float)*keypoints.size(), cudaMemcpyHostToDevice);
	

	//cudaMemcpy(_use_orientation, &use_orientation, sizeof(bool), cudaMemcpyHostToDevice);
	
	cudaMemcpy(gpu_individual, individual, sizeof(int)*512*4, cudaMemcpyHostToDevice);
	//cudaMemcpy(_descriptors, &aux, sizeof(cuda::GpuMat), cudaMemcpyHostToDevice);
	
	
	
	
	//cudaMemcpy(_sum, &aux, sizeof(KeyPoint)*keypoints.size(), cudaMemcpyHostToDevice);
	pixelTests64_kernel<<<keypoints.size(),64>>>(_sum, _x, _y,descriptors, _use_orientation, gpu_individual);
	cudaDeviceSynchronize();
	
	cudaFree(gpu_individual); cudaFree(_x); cudaFree(_y);
	//std::cout << descriptors.size() << std::endl;
	//cudaMemcpy(aux, _descriptors, sizeof(Mat), cudaMemcpyDeviceToHost);

	//Mat a;
	//_descriptors.download(a);
	//descriptors.assign(a);
	//std::cout << a << std::endl;
	//std::cout << descriptors.getMat().size() << std::endl;
	//auto finish = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double, std::micro> elapsed = finish - start;
	//std::cout << "elapsed " << elapsed.count() << std::endl;
	
	//exit(-1);
	
}



void GriefDescriptorExtractorImpl::evolve(uint ng){
	
	for(int g = 0; g < ng; g++){
		
		auto start = std::chrono::high_resolution_clock::now();
		for(int i = 0; i < N_pop; i++){
			
			mutate(i);
			
			crossover(i);
			
			if(is_infeasible())
				repair(i);
			
			selection(i);
			//std::cout << i << std::endl;

		}//exit(-1);
		
		std::cout <<  get_best_fit() << std::endl;
		
		bfit.emplace_back(get_best_fit());
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> elapsed = finish - start;
		std::cout << "Gen " << g+1 << ": Elapsed time: " << elapsed.count() << " ms." << std::endl;
		
		
	}
}

void GriefDescriptorExtractor::evolve(uint ng){
	
}

std::vector<float> GriefDescriptorExtractorImpl::gbfit(){
	return bfit;
}

std::vector<float> GriefDescriptorExtractor::gbfit(){
	return std::vector<float>{0,0};
}



//void GriefDescriptorExtractor::plot_convergence(){
//}

GriefDescriptorExtractorImpl::GriefDescriptorExtractorImpl( int bytes, bool use_orientation, EvalFunction evaluation, 
															int N_pop, float cr, float jr, float F, int mutation_algorithm, int crossover_algorithm) :
	bytes_(bytes), 
	DE(N_pop, std::vector<int>{bytes*8, 4}, cr, jr, evaluation, F, MAXIMIZATION, std::vector<int>{-24, 24}, mutation_algorithm, crossover_algorithm)
{
	this->N_pop = N_pop;
	this->jr = jr;

	load("test_pairs.brief");
	
	use_orientation_ = use_orientation;
	switch (bytes)
	{
		case 16:
			//test_fn_ = pixelTests16;
			break;
		case 32:
			//test_fn_ = pixelTests32;
			break;
		case 64:
			//test_fn_ = pixelTests64;
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
			//test_fn_ = pixelTests16;
			break;
		case 32:
			//test_fn_ = pixelTests32;
			break;
		case 64:
			//test_fn_ = pixelTests64;
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

void GriefDescriptorExtractorImpl::compute(InputArray image,
										   std::vector<KeyPoint>& keypoints,
										   cuda::GpuMat& descriptors)
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
	pixelTests64(sum, keypoints, descriptors, use_orientation_);
}

void GriefDescriptorExtractor::compute(InputArray image,
										   std::vector<KeyPoint>& keypoints,
										   cuda::GpuMat& descriptors)
{
}

void GriefDescriptorExtractorImpl::setInd(Eigen::MatrixXd new_individual){
	load("test_pairs.brief");
	std::cout << individual[0][0] << individual[0][1] << individual[0][2] << individual[0][3]  << std::endl;
	//for(int i = 0; i < bytes_*8; i++){
	//	for(int j=0; j<4; j++){
	//		individual[i][j] = new_individual(i,j);
	//	}
	//}
}

void GriefDescriptorExtractor::setInd(Eigen::MatrixXd new_individual){
}



}
} // namespace cv