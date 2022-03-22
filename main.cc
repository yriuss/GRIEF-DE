#include "matplotlibcpp.h"
#include "DE/DE.h"
#include <stdio.h>
#include<cmath>
#include<Eigen/Dense>
#include <vector>
#include<iostream>
#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>


#include "grief.h"



#define CROSSCHECK true 
#define VERTICAL_LIMIT 100
#define MAX_SEASONS 100
#define MAX_LOCATIONS 1000
#define WINDOW_SIZE 48 

char fileInfo[1000];
int numExchange = 10;
int runs = 0;
bool save=false;
bool draw=false; 
bool matchFail=false; 
using namespace std;
using namespace cv;
unsigned int n;
float distance_factor = 1.0;
int griefDescriptorLength= 512;
char dataset[1000];
char season[1000][1000]; 


int* offsetX;

typedef struct{
	int id;
	int value;
}TRating;

TRating griefRating[1024];



int numLocations = 0;
int numDisplacements = 0;
int numSeasons = 0;

typedef std::vector<Mat> Vector;
typedef std::vector<Vector> ImgMat;


int load(Eigen::MatrixXd& mat, std::string fileName) {
    
    using namespace std;
    ifstream file(fileName);

    std::string line;
    uint16_t i = 0, j = 0;
    bool successful=false;
    std::string cell;
    while (std::getline(file, line)) {
        //std::cout << line;
        std::vector<int> v;
        istringstream is(line);
        while (std::getline(is, cell, ' ')) {
			
            mat(i,j) = std::stoi(cell);
            j++;
        }
        successful=true;
        i++;
        j = 0;
    }
    
    file.close();
    return successful;
}

int x,y; 
FILE *displacements;

std::string CURRENT_DIR = get_current_dir_name();




void distinctiveMatch(const Mat& descriptors1, const Mat& descriptors2, vector<DMatch>& matches, bool crossCheck=false)
{
	Ptr<DescriptorMatcher> descriptorMatcher;
	vector<vector<DMatch> > allMatches1to2, allMatches2to1;

	descriptorMatcher = new BFMatcher(cv::NORM_HAMMING, false);
	descriptorMatcher->knnMatch(descriptors1, descriptors2, allMatches1to2, 2);

	if (!crossCheck)
	{
		for(unsigned int i=0; i < allMatches1to2.size(); i++)
		{
			if (allMatches1to2[i].size() == 2)
			{ 
				if (allMatches1to2[i][0].distance < allMatches1to2[i][1].distance * distance_factor)
				{
					DMatch match = DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
					matches.push_back(match);
				}
			}         
			else if (allMatches1to2[i].size() == 1)
			{
				DMatch match = DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
				matches.push_back(match);
				cout << "ERROR" << endl;
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
						DMatch match = DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
						matches.push_back(match);
					}
				}
				else if (allMatches2to1[allMatches1to2[i][0].trainIdx].size() == 1)
					if (allMatches1to2[i][0].distance  < allMatches1to2[i][1].distance * distance_factor && allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
					{
						DMatch match = DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
						matches.push_back(match);
						cout << "ERROR" << endl;
					} 
			}
			else if (allMatches2to1[allMatches1to2[i][0].trainIdx].size() == 2)
			{
				if (allMatches2to1[allMatches1to2[i][0].trainIdx][0].distance < allMatches2to1[allMatches1to2[i][0].trainIdx][1].distance * distance_factor && allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
				{
					DMatch match = DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
					matches.push_back(match);
					cout << "ERROR" << endl;
				}
			}
			else if (allMatches1to2[i][0].trainIdx == allMatches2to1[allMatches1to2[i][0].trainIdx][0].queryIdx)
			{
				DMatch match = DMatch(allMatches1to2[i][0].queryIdx, allMatches1to2[i][0].trainIdx, allMatches1to2[i][0].distance);
				matches.push_back(match);
				cout << "ERROR" << endl;
			} 

		}

	}

}


Mat dataset_imgs[100][100];

namespace plt = matplotlibcpp;

void plot_convergence(std::vector<int> x,std::vector<int> y){
    plt::plot(x,y);
    plt::show();
}


float eval(Eigen::MatrixXd individual){
	Ptr<cv::xfeatures2d::StarDetector>detector = cv::xfeatures2d::StarDetector::create(45,0,10,8,5);
	cv::Ptr<cv::xfeatures2d::GriefDescriptorExtractor> descriptor = cv::xfeatures2d::GriefDescriptorExtractor::create(64);
	descriptor->setInd(individual);
	
	std::vector<cv::DMatch> matches;
	

	int matchingTests = 0;
	int matchingFailures = 0;


	
	
	int i1,i2;
	
	bool supervised = false;
	

	
	for (int location = 0;location<numLocations;location++){
			
		// detecting keypoints and generating descriptors
		Mat descriptors[numSeasons];
		vector<KeyPoint> keypoints[numSeasons];
		KeyPoint kp;
		Mat dp;
		
		
		for (int i = 0;i<numSeasons;i++){
			sprintf(fileInfo,"%s/season_%02i/spgrid_regions_%09i.txt","GRIEF-datasets/michigan",i,location);
			detector->detect(dataset_imgs[i][location], keypoints[i]);
			descriptor->compute(dataset_imgs[i][location], keypoints[i], descriptors[i]);
		}
		
		
		

		// matching the extracted features
		for (int ik = 0;ik<numSeasons;ik++){
			for (int jk = ik+1;jk<numSeasons;jk++){
				matches.clear();
				/*if not empty*/
				if (descriptors[ik].rows*descriptors[jk].rows > 0) distinctiveMatch(descriptors[ik], descriptors[jk], matches, CROSSCHECK);

				/*are there any tentative correspondences ?*/
				int sumDev = 0;
				int numPoints = 0;

				int histMax = 0;
				int auxMax=0;
				int manualDir = 0; 
				int histDir = 0;
				int numBins = 100; 
				int granularity = 20;
				int maxS = 0;
				int domDir = 0;
				int histogram[numBins];
				int bestHistogram[numBins];
				vector<unsigned char> mask;
				if (matches.size() > 0){
					//histogram assembly
					int histogram[numBins];
					int bestHistogram[numBins];
					memset(bestHistogram,0,sizeof(int)*numBins);
					for (int s = 0;s<granularity;s++){
						memset(histogram,0,sizeof(int)*numBins);
						for( size_t i = 0; i < matches.size(); i++ )
						{
							int i1 = matches[i].queryIdx;
							int i2 = matches[i].trainIdx;
							if ((fabs(keypoints[ik][i1].pt.y-keypoints[jk][i2].pt.y))<VERTICAL_LIMIT){
								int devx = (int)(keypoints[ik][i1].pt.x-keypoints[jk][i2].pt.x + numBins/2*granularity);
								int index = (devx+s)/granularity;
								if (index > -1 && index < numBins) histogram[index]++;
							}
						}
						for (int i = 0;i<numBins;i++){
							if (histMax < histogram[i]){
								histMax = histogram[i];
								maxS = s;
								domDir = i;
								memcpy(bestHistogram,histogram,sizeof(int)*numBins);
							}
						}
					}

					
					for (int i =0;i<numBins;i++){
						if (auxMax < bestHistogram[i] && bestHistogram[i] != histMax){
							auxMax = bestHistogram[i];
						}
					}

					/*calculate dominant direction*/
					for( size_t i = 0; i < matches.size(); i++ )
					{
						int i1 = matches[i].queryIdx;
						int i2 = matches[i].trainIdx;
						if ((int)((keypoints[ik][i1].pt.x-keypoints[jk][i2].pt.x+numBins/2*granularity+maxS)/granularity)==domDir && fabs(keypoints[ik][i1].pt.y-keypoints[jk][i2].pt.y)<VERTICAL_LIMIT)
						{
							sumDev += keypoints[ik][i1].pt.x-keypoints[jk][i2].pt.x;
							numPoints++;
						}
					}
					histDir = (sumDev/numPoints);
					manualDir = offsetX[location+ik*numLocations] - offsetX[location+jk*numLocations];
					if (fabs(manualDir - histDir) > 35) matchFail = true; else matchFail = false;
					
					float realDir = histDir;
					int strength = 1;
					//if (matchFail) strength = 100;
					if (matchFail && supervised) realDir = manualDir;

					if (matchFail) matchingFailures++;
					matchingTests++;
					//if (histMax > 0) printf("\nDirection histogram %i %i %i\n",-(sumDev/histMax),histMax,auxMax); else printf("\nDirection histogram 1000 0 0\n");
				}else{
					matchFail = true;
				}
				
				//end drawing
			}
		}
	}
	
	
	std::cout << "error is " << (float)100*matchingFailures/matchingTests << std::endl;
    return matchingFailures/matchingTests;
}

int main(int argc, char ** argv){
	char filename[100];
	
	
	bool supervised = false;
	Mat tmpIm;
	int detectorThreshold = 0;
	distance_factor = 1.0;
	// process command line args
	if (argc > 2 && strcmp(argv[2],"draw")==0) draw = true;
	if (argc > 2 && strcmp(argv[2],"save")==0) save = true;
		
	/*load dataset parameters, check dataset consistency*/
	/*check the number of seasons and check for existance of the displacement files*/
	auto start = std::chrono::high_resolution_clock::now();
	do{
		sprintf(filename, (CURRENT_DIR + "/%s/season_%02i/%09i.bmp").c_str(), argv[1],numSeasons,numLocations);
		tmpIm =  imread(filename, cv::IMREAD_GRAYSCALE);
		if (tmpIm.data != NULL)
		{
			sprintf(filename,"%s/season_%02i/displacements.txt",argv[1],numSeasons);
			if (fopen(filename,"r") != NULL) numDisplacements++;
			x = tmpIm.cols;
			y = tmpIm.rows;
			numSeasons++;
		}
	}	
	while (numSeasons < MAX_SEASONS && tmpIm.data != NULL);
	
	
	if (numDisplacements > 0 && numDisplacements < numSeasons){
		printf("WARNING: Dataset is annotated only partially (check if the ""displacements.txt"" files exist in every ""season_nn"" directory). Ignoring hand annotation \n");
	}

	if(numSeasons == 0){
		std::cout << "Error: no images in path!" << std::endl;
		exit(EXIT_FAILURE);
	}
	supervised = (numDisplacements == numSeasons);
	
	/*check the number of locations*/
	do{
		sprintf(filename, (CURRENT_DIR + "/%s/season_%02i/%09i.bmp").c_str(), argv[1],0,numLocations++);
		tmpIm =  imread(filename, cv::IMREAD_GRAYSCALE);
	}while (numLocations < MAX_LOCATIONS && tmpIm.data != NULL);

	numLocations--;

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = finish - start;
	printf("Dataset: %ix%i images from %i seasons and %i places, annotated %i, loadTime %f\n",x,y,numSeasons,numLocations,numDisplacements, elapsed.count());

	
	
	
	/*check the dataset consistency*/
	start = std::chrono::high_resolution_clock::now();

	//rows are the num of seasons and columns are the num of locations
	offsetX = (int*)malloc(sizeof(int)*numSeasons*numLocations);
	int offsetY[numSeasons*numLocations];
	
	for (int i=0;i<numSeasons;i++)
	{
		for (int j=0;j<numLocations;j++)
		{
			sprintf(filename,(CURRENT_DIR + "/%s/season_%02i/%09i.bmp").c_str(), argv[1],i,j);
			dataset_imgs[i][j] =  imread(filename, cv::IMREAD_GRAYSCALE);
			
			if (dataset_imgs[i][j].empty()) {
				fprintf(stderr,"ERROR: Image %s could not be loaded. \n",filename);
				exit(EXIT_FAILURE);
			}		
			if (dataset_imgs[i][j].cols != x || dataset_imgs[i][j].rows != y) {
				fprintf(stderr,"ERROR: Incosistent dataset image dimensions. Image %s is %ix%i, but we expect %ix%i.\n",filename,tmpIm.cols,tmpIm.rows,x,y);
				exit(EXIT_FAILURE);
			}
		}
		if (supervised){
			sprintf(filename,"%s/season_%02i/displacements.txt",argv[1],i);
			displacements = fopen(filename,"r");
			int aX,aY,aR;
			for (int j = 0;j<numLocations;j++){
				aR = fscanf(displacements,"%i %i\n",&aX,&aY);
				offsetX[j+i*numLocations] = aX;
				offsetY[j+i*numLocations] = aY;
				if (aR != 2){
					fprintf(stderr,"ERROR: Error reading %s file on line %i.\n",filename,j);
					exit(EXIT_FAILURE);
				}
				if (abs(aX) > x || abs(aY) > y) fprintf(stderr,"WARNING: %s file on line %i indicates displacements %i %i, but image dimensions are just %i %i.\n",filename,j,aX,aY,x,y);
			}	
			fclose(displacements);
		}
	}
	finish = std::chrono::high_resolution_clock::now();
	elapsed = finish - start;
	std::cout << "Dataset consistency check OK. Time " << elapsed.count() << " ms.\n" << std::endl;
	Eigen::MatrixXd individual(512,4);

	load(individual, "test_pairs.brief");
	//UserData data = {.numSeasons=numSeasons, .numLocations=numLocations};
	//int count = 0;
	//for(int i=0; i< numSeasons;i++){
	//	for(int j =0; j<numLocations; j++){
	//		dataset_imgs[i][j] = dataset_imgs[i][j];
	//		offsetX[count] = offsetX[count];
	//		offsetY[count] = offsetY[count];
	//		count++;
	//	}
	//}


    cv::Ptr<cv::xfeatures2d::GriefDescriptorExtractor> grief_descriptor = cv::xfeatures2d::GriefDescriptorExtractor::create(64, false, eval, 30);
	
	grief_descriptor->evolve(1000);
	
	grief_descriptor->get_b_fit();
    return 0;
}