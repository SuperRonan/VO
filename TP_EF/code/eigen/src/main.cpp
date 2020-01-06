#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <string>
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

template <class T>
double error(vpImage<T> const& I, vpImage<T> const& J)
{
    double res(0);
    for(size_t i=0; i<I.getSize(); ++i)
    {
        double pix_dif = J.bitmap[i] - I.bitmap[i];
        res += pix_dif * pix_dif;
    }

    return res / I.getSize();
}

// std::string operator+(std::string const& str, int i)
// {
//     char numstr[21]; // enough to hold all numbers up to 64-bits
//     sprintf(numstr, "%d", i);
//     return str + numstr;
// }

template <class T>
T pickOneElemenent(std::vector<T> & vec)
{
    assert(!vec.empty());
    int index = rand() % vec.size(); // it is a bit biased, but it will do
    T res = vec[index];
    vec[index] = vec.back();
    vec.resize(vec.size()-1);
    return res;
}




void test_k_error_sum()
{
    const std::vector<std::string> paths = buildPathImagesAttFaces();
    const int N = paths.size();
    const std::vector<int> K = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 
        10, 15, 20, 25, 
        30, 40, 50, 60, 70, 80, 90, 
        100, 125, 150, 175, 200, 225, 250, 275, 
        300, 330, 360,
        };
    std::vector<double> errors, ev_sums;
    errors.reserve(K.size());
    ev_sums.reserve(K.size());
    EigenFacesDB egdb;
    egdb.preBuild(paths);
    egdb.writeMeanImage(std::string("res/mean") + std::string(".png"));
    egdb.writeEigenFacesImage("res/", 30);

    for(int k : K)
    {
        std::cout<<"-----------------------"<<std::endl<<k<<std::endl<<"-----------------------"<<std::endl;
        
        double ev_sum = egdb.buildBDFaces(k);
        ev_sums.push_back(ev_sum);

        double avg_error=0;
        for(std::string test_img_path : paths)
        {
            vpImage<unsigned char> test_img;
            vpImageIo::read(test_img, test_img_path);
            
            vpColVector W = egdb.W(test_img);
            vpImage<unsigned char> reconstructed_img = egdb.Jp(W);
            
            double img_error = error(reconstructed_img, test_img);
            
            avg_error += img_error;
            //vpImageIo::write(reconstructed_img, std::string("res/reconstructed_img") + N + std::string(".png"));
            
        }

        avg_error /= N;
        errors.push_back(avg_error);
        
    }

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"k\te\t\tevs"<<std::endl;
    for(int i=0; i<K.size(); ++i)
    {
        std::cout<<K[i]<<"\t"<<errors[i]<<"\t"<<ev_sums[i]<<std::endl;
    }
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"k"<<std::endl;
    for(int i=0; i<K.size(); ++i)
    {
        std::cout<<K[i]<<std::endl;
    }
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"e"<<std::endl;
    for(int i=0; i<K.size(); ++i)
    {
        std::cout<<errors[i]<<std::endl;
    }
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"evs"<<std::endl;
    for(int i=0; i<K.size(); ++i)
    {
        std::cout<<ev_sums[i]<<std::endl;
    }
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"Done!"<<std::endl;
}

#define PRINT(var) std::cout << #var << ": " << var << std::endl;

void test_recognition(int iter=20, int k = 20)
{
    const std::vector<std::string> cpaths = buildPathImagesAttFaces();
    double mean=0, max = 0, min = 1e20;
    for(int i=0; i<iter; ++i)
    {
        std::vector<std::string> paths = cpaths;
        std::string test_img_path = pickOneElemenent(paths);
        EigenFacesDB egdb;
        egdb.preBuild(paths);
        egdb.buildBDFaces(k);
        vpImage<unsigned char> img;
        vpImageIo::read(img, test_img_path);
        double theta;
        egdb.closestImage(img, &theta);
        mean += theta;
        max = std::max(max, theta);
        min = std::min(min, theta);
        //PRINT(theta);
    }
    mean = mean / iter;
    PRINT(mean);
    PRINT(max);
    PRINT(min);
}





int main()
{
    srand(time(NULL));
    omp_set_num_threads(16);

    //test_k_error_sum();

    test_recognition(100, 50);

	return 0;
}
