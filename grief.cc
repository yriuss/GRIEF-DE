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



Ptr<GriefDescriptorExtractor> GriefDescriptorExtractor::create( int bytes, bool use_orientation )
{
    return makePtr<GriefDescriptorExtractorImpl>(bytes, use_orientation );
}

int GriefDescriptorExtractorImpl::load(IndMat& mat, std::string fileName) {
    
    using namespace std;
    ifstream file(fileName);

    std::string line;
    uint16_t i = 0, j = 0;
    bool successful=false;
    std::string cell;
    while (std::getline(file, line)) {
        //std::cout << line;
        Vector v;
        istringstream is(line);
        while (std::getline(is, cell, ' ')) {
            mat[i][j] = std::stoi(cell);
            j++;
        }
        successful=true;
        i++;
        j = 0;
    }
    
    file.close();
    return successful;
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

static void pixelTests16(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, IndMat& individual)
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

static void pixelTests32(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, IndMat& individual)
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

static void pixelTests64(InputArray _sum, const std::vector<KeyPoint>& keypoints, OutputArray _descriptors, bool use_orientation, IndMat& individual)
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

GriefDescriptorExtractorImpl::GriefDescriptorExtractorImpl(int bytes, bool use_orientation) :
    bytes_(bytes), test_fn_(NULL), individual(bytes*8, Vector(4))
{

    load(individual, "test_pairs.brief");
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
    std:: cout << individual.size() << " " << individual[0].size();
    for(int i = 0; i < individual.size(); i++){
        for(int j=0; j < individual[i].size(); j++)
            std::cout << individual[i][j] << " ";
        std::cout << std::endl;
    }

}

void GriefDescriptorExtractor::getInd(){
    
}

void GriefDescriptorExtractorImpl::compute(InputArray image,
                                           std::vector<KeyPoint>& keypoints,
                                           OutputArray descriptors)
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

}
} // namespace cv