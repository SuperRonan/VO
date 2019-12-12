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

#define OMP_PARALLEL_FOR _Pragma("omp parallel for")
//#define OMP_PARALLEL_FOR

class EigenFacesDB
{
private:
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

public:
	/**
	 *  Constructeurs et destructeurs
	 */
	EigenFacesDB(){}
	virtual ~EigenFacesDB(){}

	// Construire notre base de donnee
	void buildBDFaces(const std::vector<std::string>& paths, int maxEigenFace = 50)
	{
		m_maxEigenFace = std::min(maxEigenFace, (int)paths.size());
    
		// On calcul les attibuts de l'image
		vpImage<unsigned char> I;
		vpImageIo::read(I,*paths.begin());
		m_w = I.getWidth();
		m_h = I.getHeight();
		m_size = m_h*m_w;
		
		std::cout << " * Caracteristique de l'images : " << m_h << "x" << m_w << std::endl;
		std::cout << " * Nombre d'image de la base : " << paths.size()<< std::endl;
		std::cout << " * Nombre de U : " << m_maxEigenFace << std::endl;
		
		// Creation du vpColVector pour le mean face
		std::cout << "[INFO] Creation Mean images ...." << std::endl;		
		buildMeanImage(paths);
		
		// Calcul de la matrice A
		std::cout << "[INFO] Calcul de A ... " << std::endl;
		vpMatrix A(paths.size(), m_size);
		for(int k=0; k<paths.size(); ++k)
		{
			vpImageIo::read(I, paths[k]);
			fillCenteredImage(I, A[k]);
		}
		A = A.t();

		vpMatrix U = A, V;
		vpColVector w;
		std::cout<<"[INFO] SVD ... "<<std::endl;
		U.svd(w, V);

		m_U = U.extract(0, 0, U.getRows(), m_maxEigenFace);
		
		
		std::cout << "[INFO] Fin calcul BD ... " << std::endl;
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
	
	vpColVector W(vpImage<unsigned char> const& img)const
	{
		vpColVector vec_img(img.getSize());
		fillCenteredImage(img, vec_img.data);
		vpColVector & res = vec_img;
		res = m_U.t() * vec_img;
		return res;	
	}
	
	vpImage<unsigned char> Jp(vpColVector const& w)const
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

	// * les images des eigenFaces
	void writeEigenFacesImage(const std::string& directory = "");
	// * l'image synthetisee
	vpColVector computeSynthesisImage(const std::string& image, int nbDim);
	void writeSynthesisImage(const std::string& image, const std::string& out, int nbDim);
	void writeSynthesisError(const std::string& image, const std::string& path);
	double computeSynthesisError(const std::string& image, int K);
	// Calculer le vecteur pour identifier une image
	vpColVector computeIdVector(const std::string& image, int nbDim);
	// Pour ecrire les valeurs propres
	void writeValeursPropres(const std::string& path);
	// Matrice de distance
	double computeDistanceImage(const std::string& i1, const std::string& i2);
	double distanceIdVectors(vpColVector& id1, vpColVector& id2);
	void writeMatriceDistance(std::vector<std::string> paths, const std::string& out, int nbDim);
	
private:
	void buildMeanImage(const std::vector<std::string>& paths)
	{
		m_vMean = vpColVector(m_size);
		m_vMean = 0;
		vpImage<unsigned char> tmp;
		double normalization_factor = 1.0 / (double(paths.size()) * 255.0);
		for(int k=0; k<paths.size(); ++k)
		{
			vpImageIo::read(tmp, paths[k]);

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

	vpMatrix buildMatrixA(const std::vector<std::string>& paths);
	void computeU(vpMatrix& A, vpMatrix& eigenVec);
	void writeEigenFaceInImage(const std::string& path, vpColVector v);

	int to1D(int i, int j)const
	{
		return i * m_w + j;
	}

};

#endif /* EIGENFACESDB_H_ */
