/*
 * EigenFacesDB.h
 *
 *  Created on: Dec 17, 2010
 *      Author: beltegeuse
 */

#ifndef EIGENFACESDB_H_
#define EIGENFACESDB_H_

#include <visp/vpColVector.h>
#include <visp/vpMatrix.h>
#include <vector>
#include <string>
#include <visp3/io/vpImageIo.h>
#include <omp.h>
#include <limits>
#include <visp3/core/vpImageConvert.h>
#include <algorithm>

#define OMP_PARALLEL_FOR _Pragma("omp parallel for")
//#define OMP_PARALLEL_FOR

std::string operator+(std::string const& str, int i)
{
    char numstr[21]; // enough to hold all numbers up to 64-bits
    sprintf(numstr, "%d", i);
    return str + numstr;
}

class EigenFacesDB
{
public:
	/**
	 * Attributs
	 */
	// Vecteur qui represente notre image moyenne
	vpColVector m_vMean;
	// Matrice contenant les vecteur U
	vpMatrix m_U;
	// Vecteur contenant les valeurs propres
	vpColVector m_eigenValue;
	// Les caracteristiques de l'image
	int m_h, m_w;
	int m_size;
	// Le nombre de vecteur U
	int m_maxEigenFace;

	vpMatrix m_A, full_U;
	vpColVector m_ev;
	

	std::vector<std::vector<std::string>> m_paths;

public:
	std::vector<vpColVector> m_projected_refs;

	double theta1, theta2;

	/**
	 *  Constructeurs et destructeurs
	 */
	EigenFacesDB(){}
	virtual ~EigenFacesDB(){}

	void preBuild(const std::vector<std::vector<std::string>>& paths)
	{
		m_paths = paths;
		// On calcul les attibuts de l'image
		vpImage<unsigned char> I;
		vpImageIo::read(I, paths.front().front());
		//std::cout<<"ici"<<std::endl;
		m_w = I.getWidth();
		m_h = I.getHeight();
		m_size = m_h*m_w;

		std::cout << " * Caracteristique de l'images : " << m_h << "x" << m_w << std::endl;
		std::cout << " * Nombre de visages de la base : " << paths.size()<< std::endl;
		

		// Creation du vpColVector pour le mean face
		std::cout << "[INFO] Creation Mean images ...." << std::endl;		
		buildMeanImage(paths);
		
		// Calcul de la matrice A
		std::cout << "[INFO] Calcul de A ... " << std::endl;
		m_A = vpMatrix(std::accumulate(paths.cbegin(), paths.cend(), 0,  [](size_t res, std::vector<std::string> const& v){return res + v.size();}), m_size);
		int k=0;
		for(int face=0; face<paths.size(); ++face)
		{
			for(int instance=0; instance<paths[face].size(); ++instance)
			{
				vpImageIo::read(I, paths[face][instance]);
				fillCenteredImage(I, m_A[k]);
				++k;
			}
		}
		m_A = m_A.t();

		full_U = m_A;
		vpMatrix V;
		//std::cout<<"[INFO] SVD ... "<<std::endl;
		full_U.svd(m_ev, V);

		m_ev = m_ev / m_ev.sum();


		//std::cout<<"Done!"<<std::endl;

		m_projected_refs = std::vector<vpColVector>(k);

	}



	// Construire notre base de donnee
	// returns the sum of the k first eigen values 
	double buildBDFaces(int maxEigenFace = 50)
	{
		m_maxEigenFace = std::min(maxEigenFace, (int)m_ev.size());
		//std::cout << " * Nombre de U : " << m_maxEigenFace << std::endl;
		m_U = full_U.extract(0, 0, full_U.getRows(), m_maxEigenFace);
		
		//std::cout << "[INFO] Fin calcul BD ... " << std::endl;

		double res = 0;
		for(int i=0; i<m_maxEigenFace; ++i)
		{
			res += m_ev[i];
		}

		vpImage<unsigned char> img(m_h, m_w);
		vpColVector img_W;
		int i=0;
		for(int face=0; face < m_paths.size(); ++face)
			for(int instance=0; instance<m_paths[face].size(); ++instance)
			{
				vpImageIo::read(img, m_paths[face][instance]);
				img_W = W(img);
				m_projected_refs[i] = img_W;
				++i;
			}
		

		return res;
	}

	void fillCenteredImage(vpImage<unsigned char> const& img, double * vec)const
	{
		OMP_PARALLEL_FOR
		for(int i=0; i<m_h; ++i)
			for(int j=0; j<m_w; ++j)
			{
				double mean = m_vMean[to1D(i, j)];
				double pixel = double(img[i][j]) / 255.0;
				vec[to1D(i, j)] = pixel - mean;
			}
	}

	std::vector<double> getEigenValuesVector()const
	{
		std::vector<double> res(m_ev.size());
		std::copy(m_ev.data, m_ev.data + m_ev.size(), res.begin());
		return res;
	}

	// Ecrire des images sur le systeme
	// * L'image moyenne
	void writeMeanImage(const std::string& path) const
	{
		vpImage<unsigned char> mean = getMeanImage();
		vpImageIo::write(mean, path);
	}

	vpImage<unsigned char> getMeanImage()const
	{
		vpImage<unsigned char> res(m_h, m_w);
		OMP_PARALLEL_FOR
		for(int i=0; i<m_h; ++i)
		{
			for(int j=0; j<m_w; ++j)
			{
				double pix = m_vMean[to1D(i, j)];
				res[i][j] = pix * 255;
			}
		}
		return res;
	}
	
	inline vpColVector W(vpImage<unsigned char> const& img)const
	{
		vpColVector vec_img(img.getSize());
		fillCenteredImage(img, vec_img.data);
		vpColVector & res = vec_img;
		res = m_U.t() * vec_img;
		return res;	
	}
	
	inline vpImage<unsigned char> Jp(vpColVector const& w)const
	{
		vpImage<unsigned char> res(m_h, m_w);
		vpColVector col = m_U * w;
		for(int i=0; i<col.getRows(); ++i)
		{
			double pix = std::max(0.0, std::min(255.0, (m_vMean[i] + col[i]) * 255));
			res.bitmap[i] = pix;
		}
		return res;
	}

	std::pair<int, int> closestImage(vpImage<unsigned char> const& img, double * theta = nullptr, std::pair<int, int> except={-1, -1})const
	{
		vpColVector img_W = W(img);
		return closestImage(img_W, theta, except);
	}

	std::pair<int, int> closestImage(vpColVector const& img_W, double * theta = nullptr, std::pair<int, int> except={-1, -1})const
	{
		std::pair<int, int> best_id;
		double best_score = std::numeric_limits<double>::max();
		int i=0;
		for(int face=0; face<m_paths.size(); ++face)
		{
			if(except.first == face && except.second == -1) //skip whole face
			{
				i += m_paths[face].size();
				continue;
			}
			for(int instance=0; instance<m_paths[face].size(); ++instance)
			{
				if(except.first == face && except.second == instance) //skip instance
				{
					++i;
					continue;
				}
				vpColVector dif = m_projected_refs[i] - img_W;
				double score = dif.sumSquare();
				if(score < best_score)
				{
					best_id = {face, instance};
					best_score = score;
				}
				++i;
			}
		}
		if(theta)	*theta = best_score;
		return best_id;
	}

	// * les images des eigenFaces
	void writeEigenFacesImage(const std::string& directory = "./", int max=0)
	{
		//vpColVector eigenFace;
		vpImage<double> ef_img(m_h, m_w);
		vpImage<unsigned char> uc_img(m_h, m_w);
		max = std::min(400, max);
		for(int i=0; i<max; ++i)
		{
			//eigenFace = full_U.getCol(i);
			OMP_PARALLEL_FOR
			for(int j=0; j<ef_img.getSize(); ++j)	ef_img.bitmap[j] = full_U[j][i];
			vpImageConvert::convert(ef_img, uc_img);
			vpImageIo::write(uc_img, directory + std::string("eigenFace") + i + std::string(".png"));
		}
	}

	void centerImage(vpImage<unsigned char> & img)const
	{
		vpColVector vec(m_size);
		fillCenteredImage(img, vec.data);
		vpImage<double> tmp(m_h, m_w);
		OMP_PARALLEL_FOR
		for(int i=0; i<tmp.getSize(); ++i)	tmp.bitmap[i] = vec[i];
		vpImageConvert::convert(tmp, img);
	}

	int recognize_face(vpImage<unsigned char> const& img, int except=-1)const
	{
		vpColVector w = W(img);
		return recognize_face(w, except);
	}
	// returns -2 if w is not a face image
	// returns -1 if w is a face, but not in the database
	int recognize_face(vpColVector const& w, int except=1)const
	{
		int naive_res = recognize_face_naive(w, except);
		if(naive_res < 0)
		    return naive_res;
		double alpha=0.9, beta=0.5, gamma=0.1;
		const auto score = [alpha, beta, gamma](std::vector<double> const& thetas, double sum, double min)
		{
			double mean = sum / thetas.size();
			double var = std::accumulate(thetas.cbegin(), thetas.cend(), 0, [mean](double x, double y){double dif=y - mean; return x + dif * dif;}) / thetas.size();
			return alpha * min + beta * mean + gamma * std::sqrt(var);
		};
		int face_id = -1;
		double best_score = std::numeric_limits<double>::max();
		int i=0; 
		for(int face=0; face<m_paths.size(); ++face)
		{
			std::vector<double> thetas;
			double sum = 0;	
			double min = std::numeric_limits<double>::max();
			for(int instance=0; instance<m_paths[face].size(); ++instance)
			{
				if(instance == except)	{++i; continue;}
				double error = (m_projected_refs[i] - w).sumSquare();
				min = std::min(min, error);
				sum += error;
				thetas.push_back(error);
				++i;
			}
			if(min < theta2)
			{
				double face_score = score(thetas, sum, min);
				if(face_score < best_score)
				{
					face_id = face;
					best_score = face_score;
				}
			}
		}
		

		return face_id;
	}

	int recognize_face_naive(vpImage<unsigned char> const& img, int except=-1)const
	{
		vpColVector w = W(img);
		return recognize_face_naive(w, except);
	}

	// returns -2 if w is not a face image
	// returns -1 if w is a face, but not in the database
	int recognize_face_naive(vpColVector const& w, int except=1)const
	{
		bool w_is_face = false;
		int face_id = -1;
		double best_score = std::numeric_limits<double>::max();
		int i=0; 
		for(int face=0; face<m_paths.size(); ++face)
		{
			for(int instance=0; instance<m_paths[face].size(); ++instance)
			{
				if(instance == except)	{++i; continue;}
				double error = (m_projected_refs[i] - w).sumSquare();
				w_is_face = w_is_face || error < theta2;
				
				if(error < best_score)
				{
					face_id = face;
					best_score = error;
				}

				++i;
			}
		}
		if(best_score < theta1)
			return face_id;
		else if(w_is_face)
			return -1;
		else
			return -2;
	}

	
private:
	void buildMeanImage(const std::vector<std::vector<std::string>>& paths)
	{
		m_vMean = vpColVector(m_size);
		m_vMean = 0;
		vpImage<unsigned char> tmp;
		int N = std::accumulate(paths.cbegin(), paths.cend(), 0,  [](size_t res, std::vector<std::string> const& v){return res + v.size();});
		double normalization_factor = 1.0 / (N * 255.0);
		for(int face=0; face<paths.size(); ++face)
			for(int instance=0; instance < paths[face].size(); ++instance)
		{
			
			vpImageIo::read(tmp, paths[face][instance]);

			OMP_PARALLEL_FOR
			for(int i=0; i<m_h; ++i)
			{
				for(int j=0; j<m_w; ++j)
				{
					double pix = tmp[i][j];
					pix *= normalization_factor;
					m_vMean[to1D(i, j)] += pix;
				}
			}
		}
	}

	int to1D(int i, int j)const
	{
		return i * m_w + j;
	}

};

#endif /* EIGENFACESDB_H_ */
