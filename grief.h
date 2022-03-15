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

#include "DE.h"

namespace cv{

namespace xfeatures2d
{


class CV_EXPORTS_W GriefDescriptorExtractor : public Feature2D
{
public:
    CV_WRAP static Ptr<GriefDescriptorExtractor> create( int bytes = 32, bool use_orientation = false );
    CV_WRAP virtual void getInd( );
};

/*
 * Grief Descriptor
 */
class GriefDescriptorExtractorImpl : public GriefDescriptorExtractor, DE
{
public:
    enum { PATCH_SIZE = 48, KERNEL_SIZE = 9 };

    // bytes is a length of descriptor in bytes. It can be equal 16, 32 or 64 bytes.
    GriefDescriptorExtractorImpl( int bytes = 32, bool use_orientation = false );
    int load(IndMat& mat, std::string fileName);
    virtual void read( const FileNode& ) CV_OVERRIDE;
    virtual void write( FileStorage& ) const CV_OVERRIDE;
    
    virtual int descriptorSize() const CV_OVERRIDE;
    virtual int descriptorType() const CV_OVERRIDE;
    virtual int defaultNorm() const CV_OVERRIDE;

    virtual void compute(InputArray image, std::vector<KeyPoint>& keypoints, OutputArray descriptors) CV_OVERRIDE;
    virtual void getInd() CV_OVERRIDE;
    

protected:
    typedef void(*PixelTestFn)(InputArray, const std::vector<KeyPoint>&, OutputArray, bool use_orientation, IndMat& individual);
    IndMat individual;
    int bytes_;
    bool use_orientation_;
    PixelTestFn test_fn_;
};
}
}
