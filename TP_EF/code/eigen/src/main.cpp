#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <visp3/core/vpImage.h>

#include "EigenFacesDB.h"

void  writeEigenFaceInImage(int m_w, int m_h, const std::string& path, vpColVector &v)
{
    double maxV = v.getMaxValue();
    double minV = v.getMinValue();
    vpImage<unsigned char> I(m_h,m_w);
    for(int i=0; i<m_h; i++)
        for(int j = 0; j<m_w; j++)
        {
            I[i][j] = ((v[i*m_w+j]-minV)/(maxV-minV))*255;
        }
    vpImageIo::write(I, path);
}

std::vector<std::string> buildPathImagesAttFaces()
{
	std::vector<std::string> v;
	for(int nbDir=1; nbDir<=40; nbDir++)
		for(int nbImage=1; nbImage<10;nbImage++)
		{
			std::ostringstream ss;
			ss << "../../../data/s" << nbDir << "/" << nbImage << ".pgm";
			v.push_back(ss.str());
		}
	return v;
}

int main()
{
    omp_set_num_threads(4);

	std::cout << "[INFO] Construction du path ..." << std::endl;
	std::vector<std::string> paths = buildPathImagesAttFaces();
	std::cout << "[INFO] Creation de base de donnees ..." << std::endl;
	
    EigenFacesDB egdb;
    egdb.buildBDFaces(paths,100);

    egdb.writeMeanImage("res/mean.png");
    
    std::string test_img_path = "../../../data/s25/1.pgm";
    vpImage<unsigned char> test_img;
    vpImageIo::read(test_img, test_img_path);
    
    vpColVector W = egdb.W(test_img);
    vpImage<unsigned char> reconstructed_img = egdb.Jp(W);
    
    vpImageIo::write(reconstructed_img, "res/reconstructed_img.png");
	
	return 0;
}
