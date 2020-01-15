#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <string>
#include <visp3/core/vpImage.h>
#include <algorithm>
#include <numeric>

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

namespace std
{
    template <class T>
    using dvector = std::vector<std::vector<T>>;
}

inline std::dvector<std::string> buildPathImagesAttFaces(int except=0)
{
	std::dvector<std::string> v;
	for(int nbDir=1; nbDir<=40; nbDir++)
    {
        std::vector<std::string> face_vec;
		for(int nbImage=1; nbImage<11-except;nbImage++)
		{
			std::ostringstream ss;
			ss << "../../../data/s" << nbDir << "/" << nbImage << ".pgm";
			face_vec.push_back(ss.str());
		}
        v.push_back(face_vec);
    }
	return v;
}

inline std::dvector<vpImage<unsigned char>> loadAllImages(std::dvector<std::string> const& paths)
{
    std::dvector<vpImage<unsigned char>> res(paths.size());
    OMP_PARALLEL_FOR
    for(int face=0; face < paths.size(); ++face)
    {
        res[face] = std::vector<vpImage<unsigned char>>(paths[face].size());
        for(int instance=0; instance<paths[face].size(); ++instance)
        {   
            vpImageIo::read(res[face][instance], paths[face][instance]);
        }
    }
    return res;
}


inline std::dvector<vpColVector> projectAll(std::dvector<vpImage<unsigned char>> const& imgs, EigenFacesDB const& egdb)
{
    std::dvector<vpColVector> res(imgs.size());
    for(int face=0; face<imgs.size(); ++face)
    {
        res[face] = std::vector<vpColVector>(imgs[face].size());
        for(int instance=0; instance<imgs[face].size(); ++instance)
        {
            res[face][instance] = egdb.W(imgs[face][instance]);
        }
    }
    return res;
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
T pickOne(std::vector<T> & vec)
{
    assert(!vec.empty());
    int index = rand() % vec.size(); // it is a bit biased, but it will do
    T res = vec[index];
    vec[index] = vec.back();
    vec.resize(vec.size()-1);
    return res;
}


template <class T>
size_t vecvecsize(std::vector<std::vector<T>> const& vec)
{
    return std::accumulate(vec.cbegin(), vec.cend(), 0,  [](size_t res, std::vector<T> const& v){return res + v.size();});
}



template <class Stream, class T>
Stream& matlabPrint(std::vector<T> const& vec, Stream & stream, bool semicolon=true)
{
    stream << "[";
    for(T const& elem : vec)
    {
        stream << elem << ", ";
    }
    stream << "]";
    if(semicolon)
        stream << ";";
    return stream;
}

template <class Stream, class T>
Stream& exelPrint(std::vector<T> const& vec, Stream & stream)
{
    for(T const& elem : vec)
    {
        stream << elem << std::endl;
    }
    return stream;
}

template <class Stream, class T>
Stream& printVector(std::vector<T> const& vec, Stream & stream, bool type=true)
{
    if(type)
        return matlabPrint(vec, stream);
    else
        return exelPrint(vec, stream);
}

void test_k_error_sum(bool type=true)
{
    std::cout<<"hi"<<std::endl;
    const auto paths = buildPathImagesAttFaces();
    const int N = vecvecsize(paths);
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
    std::cout<<"Eigen Values: "<<std::endl;
    printVector(egdb.getEigenValuesVector(), std::cout, type)<<std::endl;
    egdb.writeMeanImage(std::string("res/mean") + std::string(".png"));
    egdb.writeEigenFacesImage("res/", 30);
    vpImage<unsigned char> test_img, reconstructed_img;
    vpColVector W;
    for(int k : K)
    {
        //std::cout<<"-----------------------"<<std::endl;
        std::cout<<"\rk: "<<k<<std::flush;
        //std::cout<<"-----------------------"<<std::endl;
        
        double ev_sum = egdb.buildBDFaces(k);
        ev_sums.push_back(ev_sum);

        double avg_error=0;
        for(int face = 0; face < paths.size(); ++face)
        {
            for(int instance = 0; instance<paths[face].size(); ++instance)
            {
                std::string const& test_img_path = paths[face][instance];
                
                vpImageIo::read(test_img, test_img_path);
                
                W = egdb.W(test_img);
                reconstructed_img = egdb.Jp(W);
                
                double img_error = error(reconstructed_img, test_img);
                
                avg_error += img_error;
                //vpImageIo::write(reconstructed_img, std::string("res/reconstructed_img") + N + std::string(".png"));
                
            }
        }

        avg_error /= N;
        errors.push_back(avg_error);
        std::cout<<"    ";
    }
    std::cout<<std::endl;

    std::cout<<"k"<<std::endl;
    printVector(K, std::cout, type) << std::endl;

    std::cout<<"error"<<std::endl;
    printVector(errors, std::cout, type)<<std::endl;

    std::cout<<"evs"<<std::endl;
    printVector(ev_sums, std::cout, type)<<std::endl;

    std::cout<<"Done!"<<std::endl;
}

#define PRINT(var) std::cout << #var << ": " << var << std::endl;

void test_recognition(int iter=20, int k = 20, bool type=true)
{
    const auto cpaths = buildPathImagesAttFaces();
    double mean=0, max = 0, min = 1e20;
    for(int i=0; i<iter; ++i)
    {
        auto paths = cpaths;
        auto face = pickOne(paths);
        auto test_img_path = pickOne(face);
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

void stat_recognition(bool type)
{
    const auto paths = buildPathImagesAttFaces();
    const int N = vecvecsize(paths);
    const std::vector<int> K = {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 
        10, 15, 20, 25, 
        30, 40, 50, 60, 70, 80, 90, 
        100, 125, 150, 175, 200, 225, 250, 275, 
        300, 330, 360,
        };
    std::vector<double> img_errors, vec_errors;
    EigenFacesDB egdb;
    
    egdb.preBuild(paths);

    egdb.writeMeanImage(std::string("res/mean") + std::string(".png"));
    egdb.writeEigenFacesImage("res/", N);

    for(int k : K)
    {
        //std::cout<<"-----------------------"<<std::endl;
        std::cout<<"\rk: "<<k<<std::flush;
        //std::cout<<"-----------------------"<<std::endl;
        
        double ev_sum = egdb.buildBDFaces(k);
        double image_space_avg_error=0;
        double vector_space_avg_error=0;
        for(int face = 0; face<paths.size(); ++face)
            for(int instance = 0; instance < paths[face].size(); ++instance)
        {
            
            std::string test_img_path = paths[face][instance];
            vpImage<unsigned char> test_img;
            vpImageIo::read(test_img, test_img_path);
                
            vpColVector W = egdb.W(test_img);
            
            double theta;
            auto recognized_id = egdb.closestImage(W, &theta, {face, instance});

            vector_space_avg_error += theta / k;

            vpImage<unsigned char> recognized_img;
            
            vpImageIo::read(recognized_img, paths[recognized_id.first][recognized_id.second]);
            double img_error = error(recognized_img, test_img);
            image_space_avg_error += img_error;
        }

        image_space_avg_error /= N;
        vector_space_avg_error /= N;

        img_errors.push_back(image_space_avg_error);
        vec_errors.push_back(vector_space_avg_error);
        std::cout<<"        ";
    }
    std::cout<<std::endl;

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"k"<<std::endl;
    printVector(K, std::cout, type);
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"image space error"<<std::endl;
    printVector(img_errors, std::cout, type);
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"------------------------------"<<std::endl;
    std::cout<<"vector space error"<<std::endl;
    printVector(vec_errors, std::cout, type);
    std::cout<<"------------------------------"<<std::endl;

    std::cout<<"Done!"<<std::endl;
}

void test_matrix(const std::vector<int> & K = {20}, bool type=true)
{
    const auto paths = buildPathImagesAttFaces();
    const auto images = loadAllImages(paths);
    int N = vecvecsize(paths);

    vpMatrix matrix(N, N);
    matrix = 0;

    EigenFacesDB egdb;
    egdb.preBuild(paths);
    vpImage<unsigned char> img(egdb.m_h, egdb.m_w), img2(egdb.m_h, egdb.m_w);

    
    for(int k : K)
    {
        std::cout<<"k: "<<k<<std::endl;
        egdb.buildBDFaces(k);
        //vpColVector W(k), W2(k);
        const auto Ws = projectAll(images, egdb);
        int i=0;
        for(int face=0; face<paths.size(); ++face)
        {
            std::cout<<"\rface: "<<face<<std::flush;
            for(int instance=0; instance<paths[face].size(); ++instance)
            {
                //vpImageIo::read(img, paths[face][instance]);
                //W = egdb.W(images[face][instance]);
                int j=0;
                for(int face2=0; face2 < paths.size(); ++face2)
                {
                    for(int instance2=0; instance2<paths[face2].size(); ++instance2)
                    {
                        //vpImageIo::read(img2, paths[face2][instance2]);
                        //W2 = egdb.W(images[face][instance]);

                        double dist = (Ws[face][instance] - Ws[face2][instance2]).sumSquare();
                        matrix[i][j] = dist;
                        ++j;
                    }
                }
                ++i;
            }
        }
        std::cout<<std::endl<<"Writing matrix "<<k<<std::endl;
        std::string file_name = std::string("matrix") + k + std::string(".mat");
        std::ofstream file(file_name);
        file << "M"<<k<<" = ";
        matrix.matlabPrint(file);
        file << ";"<<std::endl;
        file.close();
    }
}

void computeCenteredImages(int N = 4)
{
    const auto paths = buildPathImagesAttFaces();
    EigenFacesDB egdb;
    egdb.preBuild(paths);
    vpImage<unsigned char> img(egdb.m_h, egdb.m_w);
    for(int i=0; i<N; ++i)
    {
        int face = rand() % paths.size();
        int instance = rand() % paths[face].size();
        vpImageIo::read(img, paths[face][instance]);

        egdb.centerImage(img);
        vpImageIo::write(img, std::string("res/centered_") + face + std::string("-") + instance + std::string(".png"));
    }
}

void savePng(int N = 4)
{
    const auto paths = buildPathImagesAttFaces();
    vpImage<unsigned char> img;
    for(int face=0; face < paths.size(); ++face)
    {
        for(int instance=0; instance < paths[face].size(); ++instance)
        {
            vpImageIo::read(img, paths[face][instance]);
            vpImageIo::write(img, std::string("data/face_") + face + std::string("-") + instance + std::string(".png"));
        }
    }
}

template <class numeric>
numeric evaluate_threshold(std::vector<numeric> & X)
{
    std::sort(X.begin(), X.end());
    numeric min = X.front(), max = X.back();
    numeric mean = std::accumulate(X.cbegin(), X.cend(), 0) / X.size();
    numeric var = std::accumulate(X.cbegin(), X.cend(), 0, [mean](numeric sum, numeric next){numeric dif = mean - next; return dif * dif;}) / X.size();
    numeric sigma = std::sqrt(var);

    return mean;
}

void evaluate_theta(std::vector<int> const& K = {20}, bool type=true)
{
    const auto paths = buildPathImagesAttFaces();
    const auto images = loadAllImages(paths);
    EigenFacesDB egdb;
    egdb.preBuild(paths);
    for(int k : K)
    {
        std::cout<<"k: "<<k<<std::endl;
        egdb.buildBDFaces(k);
        const auto Ws = projectAll(images, egdb);
        std::vector<double> thetas_same_face, thetas_is_face;
        for(int face=0; face<paths.size(); ++face)
        {
            for(int instance=0; instance<paths[face].size(); ++instance)
            {
                const vpColVector & W = Ws[face][instance];
                //find the best theta for the same face
                double theta_same_face = std::numeric_limits<double>::max(), theta_is_face=std::numeric_limits<double>::max();
                for(int f=0; f<paths.size(); ++f)
                {
                    for(int i=0; i<paths[f].size(); ++i)
                    {
                        const vpColVector & w = Ws[f][i];
                        if(f == face)
                        {
                            //same face
                            if(i == instance)   continue;
                            theta_same_face = std::min(theta_same_face, (W - w).sumSquare());
                        }
                        else
                        {
                            theta_is_face = std::min(theta_is_face, (W - w).sumSquare());
                        }
                    }
                }
                thetas_same_face.push_back(theta_same_face);
                thetas_is_face.push_back(theta_is_face);
            }
        }

        std::sort(thetas_same_face.begin(), thetas_same_face.end());
        std::sort(thetas_is_face.begin(), thetas_is_face.end());

        if(false)
        {
            std::cout<<"theta is face: "<<std::endl;
            printVector(thetas_is_face, std::cout, type)<<std::endl;

            std::cout<<"theta same face: "<<std::endl;
            printVector(thetas_same_face, std::cout, type)<<std::endl;
        }
        else
        {
            std::ofstream file(std::string("theta_") + k + std::string(".mat"));
            file << "TIF = ";
            printVector(thetas_is_face, file, type)<<std::endl<<std::endl;

            file<<"TSF = ";
            printVector(thetas_same_face, file, type)<<std::endl;
            file.close();
        }

        double theta_is_face = evaluate_threshold(thetas_is_face);

        double theta_same_face = evaluate_threshold(thetas_same_face);

        PRINT(theta_is_face);
        PRINT(theta_same_face);
        
        
    }
}

void evaluate_theta_alternative(std::vector<int> const& K = {20}, bool type=true)
{
    const auto paths = buildPathImagesAttFaces();
    const auto images = loadAllImages(paths);
    EigenFacesDB egdb;
    egdb.preBuild(paths);
    for(int k : K)
    {
        std::cout<<"k: "<<k<<std::endl;
        egdb.buildBDFaces(k);
        const auto Ws = projectAll(images, egdb);
        std::vector<double> thetas_same_face, thetas_is_face;
        for(int face=0; face<paths.size(); ++face)
        {
            for(int instance=0; instance<paths[face].size(); ++instance)
            {
                const vpColVector & W = Ws[face][instance];
                //find the best theta for the same face
                double theta_same_face = std::numeric_limits<double>::max(), theta_is_face=std::numeric_limits<double>::max();

                auto id = egdb.closestImage(W, &theta_same_face, {face, instance});
                if(id.first == face)//correct
                {
                    thetas_same_face.push_back(theta_same_face);
                }
                
                id = egdb.closestImage(W, &theta_is_face, {face, -1});
                thetas_is_face.push_back(theta_is_face);
            }
        }

        std::sort(thetas_same_face.begin(), thetas_same_face.end());
        std::sort(thetas_is_face.begin(), thetas_is_face.end());

        if(false)
        {
            std::cout<<"theta is face: "<<std::endl;
            printVector(thetas_is_face, std::cout, type)<<std::endl;

            std::cout<<"theta same face: "<<std::endl;
            printVector(thetas_same_face, std::cout, type)<<std::endl;
        }
        else
        {
            std::ofstream file(std::string("theta_alter_") + k + std::string(".mat"));
            file << "TIF = ";
            printVector(thetas_is_face, file, type)<<std::endl<<std::endl;

            file<<"TSF = ";
            printVector(thetas_same_face, file, type)<<std::endl;
            file.close();
        }

        double theta_is_face = evaluate_threshold(thetas_is_face);

        double theta_same_face = evaluate_threshold(thetas_same_face);

        PRINT(theta_is_face);
        PRINT(theta_same_face);
        
        
    }
}


int main()
{
    srand(time(NULL));
    omp_set_num_threads(16);

    bool type = 1; //matlab

    //test_k_error_sum(type);

    //test_recognition(100, 50, type);

    //stat_recognition(type);

    //test_matrix({1, 10, 20, 50, 100, 400}, type);

    evaluate_theta({1, 10, 20, 50, 100}, type);

    evaluate_theta_alternative({1, 10, 20, 50, 100}, type);

    //computeCenteredImages();

    //savePng();

    std::cout<<"Done!"<<std::endl;

	return 0;
}
