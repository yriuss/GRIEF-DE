#include <sys/time.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d.hpp>
#include <stdio.h>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include "grief.h"

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"

#include <iostream>

#define CROSSCHECK true 




float distance_factor = 1.0;

//feature matching - this can combine 'ratio' and 'cross-check' 
void distinctiveMatch(const cv::Mat& descriptors1, const cv::Mat& descriptors2, std::vector<cv::DMatch>& matches, bool crossCheck=false)
{
	cv::Ptr<cv::DescriptorMatcher> descriptorMatcher;
	std::vector<std::vector<cv::DMatch> > allMatches1to2, allMatches2to1;

	descriptorMatcher = new cv::BFMatcher(cv::NORM_HAMMING, false);
	descriptorMatcher->knnMatch(descriptors1, descriptors2, allMatches1to2, 2);

	if (!crossCheck)
	{
		for(unsigned int i=0; i < allMatches1to2.size(); i++)
		{
			if (allMatches1to2[i].size() == 2)
			{ 
				if (allMatches1to2[i][0].distance < allMatches1to2[i][1].distance * distance_factor)
				{
					cv::DMatch match = cv::DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
					matches.push_back(match);
				}
			}         
			else if (allMatches1to2[i].size() == 1)
			{
				cv::DMatch match = cv::DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
				matches.push_back(match);
				std::cout << "ERROR" << std::endl;
			}

		}
	}
	else
	{
		descriptorMatcher->knnMatch(descriptors2, descriptors1, allMatches2to1, 2);
		for(unsigned int i=0; i < allMatches1to2.size(); i++)
		{
			if (allMatches1to2[i].size() == 2)
			{ 
				if (allMatches2to1[allMatches1to2[i][0].trainIdx].size() == 2)
				{
					if (allMatches1to2[i][0].distance < allMatches1to2[i][1].distance * distance_factor && allMatches2to1[allMatches1to2[i][0].trainIdx][0].distance < allMatches2to1[allMatches1to2[i][0].trainIdx][1].distance * distance_factor && allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
					{
						cv::DMatch match = cv::DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
						matches.push_back(match);
					}
				}
				else if (allMatches2to1[allMatches1to2[i][0].trainIdx].size() == 1)
					if (allMatches1to2[i][0].distance  < allMatches1to2[i][1].distance * distance_factor && allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
					{
						cv::DMatch match = cv::DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
						matches.push_back(match);
						std::cout << "ERROR" << std::endl;
					} 
			}
			else if (allMatches2to1[allMatches1to2[i][0].trainIdx].size() == 2)
			{
				if (allMatches2to1[allMatches1to2[i][0].trainIdx][0].distance < allMatches2to1[allMatches1to2[i][0].trainIdx][1].distance * distance_factor && allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
				{
					cv::DMatch match = cv::DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
					matches.push_back(match);
					std::cout << "ERROR" << std::endl;
				}
			}
			else if (allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
			{
				cv::DMatch match = cv::DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
				matches.push_back(match);
				std::cout << "ERROR" << std::endl;
			} 

		}

	}

}

int main( int argc, char* argv[] )
{
    cv::Mat img[2];
    cv::Mat descriptors[2];
	std::vector<cv::KeyPoint> keypoints[2];
    std::string filename = "test_pairs.brief";

	
    cv::Ptr<cv::xfeatures2d::GriefDescriptorExtractor> descriptor = cv::xfeatures2d::GriefDescriptorExtractor::create(64);
    cv::Ptr<cv::xfeatures2d::StarDetector>detector = cv::xfeatures2d::StarDetector::create(45,0,10,8,5);
    
    //descriptor->getInd();

    img[0] = imread("0.bmp", cv::IMREAD_GRAYSCALE);
    img[1] = imread("1.bmp", cv::IMREAD_GRAYSCALE);


	auto start = std::chrono::high_resolution_clock::now();
	
    detector->detect(img[0], keypoints[0]);
    descriptor->compute(img[0], keypoints[0], descriptors[0]);

    detector->detect(img[1], keypoints[1]);
    descriptor->compute(img[1], keypoints[1], descriptors[1]);

	


    std::vector<cv::DMatch> matches, inliers_matches,working_matches;
    if (descriptors[0].rows*descriptors[1].rows > 0) distinctiveMatch(descriptors[0], descriptors[1], matches, CROSSCHECK);

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	std::cout << "Elapsed Time: " << elapsed.count() << " milliseconds" << std::endl;
	

    cv::Mat imA,imB,img_matches,img_matches_transposed;

    cv::transpose(img[0], imA);
	cv::transpose(img[1], imB);
	
	
	cv::namedWindow("matches", 1);
	
	cv::drawMatches(img[0], keypoints[0], img[1], keypoints[1], matches, img_matches, cv::Scalar(0, 0, 255), cv::Scalar(0, 0, 255), std::vector<char>(), cv::DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	




	cv::transpose(img_matches,img_matches_transposed);
	cv::imshow("matches", img_matches);
	cv::waitKey(0);

    //std::vector<std::vector<int> > individuo(512);
    //for ( int i = 0 ; i < 512 ; i++ )
    //   individuo[i].resize(2);
    //std::vector<int> individuo(2);
    //load(individuo, filename);

	

    return 0;
}