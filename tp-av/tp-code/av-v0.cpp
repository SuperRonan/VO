#define VP_TRACE

//! \example tutorial-viewer.cpp
//! [Include display]
#include <visp3/gui/vpDisplayD3D.h>
#include <visp3/gui/vpDisplayGDI.h>
#include <visp3/gui/vpDisplayGTK.h>
#include <visp3/gui/vpDisplayX.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/core/vpPoint.h>
#include <visp3/core/vpPlane.h>
#include <visp3/gui/vpPlot.h>
#include <visp3/core/vpMeterPixelConversion.h>
#include <visp3/core/vpCameraParameters.h>

#include <visp3/core/vpExponentialMap.h>
#include <visp3/core/vpVelocityTwistMatrix.h>



//! [Include display]
//! [Include io]
#include <visp3/io/vpImageIo.h>
//! [Include io]


#include <visp3/core/vpTime.h>

using namespace std ;


void
display(vpCameraParameters& cam, vpImage<unsigned char> &I, vpColVector &x, vpColVector &xd)
{
    for (int i = 0 ; i < x.getRows()/2 ; i++)
    {
        vpImagePoint u,ud ;
        vpMeterPixelConversion::convertPoint(cam,x[2*i],x[2*i+1],u) ;
        vpDisplay::displayPoint(I,u,vpColor::red, 2) ;

        vpMeterPixelConversion::convertPoint(cam,xd[2*i],xd[2*i+1],ud) ;
        vpDisplay::displayCircle(I,ud,10,vpColor::green) ;

    }

    vpDisplay::flush(I) ;
}


// Projection d'un point 3D sur le plane image  X(3), x(2)
void
project(vpColVector &X, vpColVector &x)
{
    x[0] = X[0] / X[2];
    x[1] = X[1] / X[2];
}

// Changement de repere bX(3), aX(3), aTb est une matrice homogène
void changeFrame(const vpColVector &bX, const vpHomogeneousMatrix &aTb,  vpColVector &aX)
{
    vpColVector bXHomo = vpColVector(bX.size() + 1);
    for (int i = 0; i < bX.size(); i++)
    {
        bXHomo[i] = bX[i];
    }
    bXHomo[bX.size()] = 1;

    vpColVector aXHomo = aTb * bXHomo;
    for (int i = 0; i < aX.size(); i++)
    {
        aX[i] = aXHomo[i];
    }
}



// Calcul de la matrice d'interaction d'un point 2D
void computeInteractionMatrix(const vpColVector &cX, double x, double y, vpMatrix &Lx)
{
    Lx = 0;
    Lx[0][0] = - 1 / cX[2];
    Lx[1][1] = - 1 / cX[2];
    Lx[0][2] = x / cX[2];
    Lx[1][2] = y / cX[2];
    Lx[0][3] = x * y;
    Lx[1][3] = 1 + y*y;
    Lx[0][4] = -(1 + x * x);
    Lx[1][4] = -x*y;
    Lx[0][5] = y;
    Lx[1][5] = -x;
}


void computeInteractionMatrixMultiPoints(const vpColVector * cX, vpColVector const& x, vpMatrix &Lx)
{
    Lx = 0;
    for(int i=0; i<x.size()/2; ++i)
    {
        Lx[2*i+0][0] = - 1 / cX[i][2];
        Lx[2*i+1][1] = - 1 / cX[i][2];
        Lx[2*i+0][2] = x[2*i+0] / cX[i][2];
        Lx[2*i+1][2] = x[2*i+1] / cX[i][2];
        Lx[2*i+0][3] = x[2*i+0] * x[2*i+1];
        Lx[2*i+1][3] = 1 + x[2*i+1] * x[2*i+1];
        Lx[2*i+0][4] = -(1 + x[2*i+0] * x[2*i+0]);
        Lx[2*i+1][4] = -x[2*i+0] * x[2*i+1];
        Lx[2*i+0][5] = x[2*i+1];
        Lx[2*i+1][5] = -x[2*i+0]; 
    }
}






void
tp2DVisualServoingOnePoint()
{

    //-------------------------------------------------------------
    // Mise en oeuvre des courbes


    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,2);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "Camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

   
    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);


    vpCameraParameters cam(400,400,300,200) ;


    //-------------------------------------------------------------

    //Definition de la scene

    //positions initiale (à tester)

    //Definition de la scene
    vpHomogeneousMatrix cTw (0,0,1,  0,0,0) ;

    //position of the point in the world frame
    vpColVector wX(3), cX(3);
    wX[0] = 0.5 ; wX[1] = 0.2 ; wX[2] = -0.5 ;


    // a chaque fois que vous verez size metter la bonne taille à la place
    int size = 2;

    vpColVector e(size) ; //
    e = 12;



    // position courante, position desiree
    vpColVector x(size), xd(size) ;

    xd = 0;
    //matrice d'interaction
    vpMatrix Lx(size,6), LxPlus(6, size);
    

    // position desirée  (au centre de l'image x=0, y=0)

    // vitesse de la camera
    vpColVector v(size) ;
    double lambda = 0.1 ;
    int iter = 0 ;
    while (fabs(e.sumSquare()) > 1e-6)
    {
        // calcul de la position des points dans l'image en fonction de cTw
       // instancier x
       
        changeFrame(wX, cTw, cX);
        project(cX, x);
        
        //calcul de l'erreur
        e = x - xd;

        // Calcul de la matrice d'interaction

        computeInteractionMatrix(cX, x[0], x[1], Lx);

        //calcul de la loi de commande v= ...

        LxPlus = Lx.pseudoInverse();
        v = -lambda * LxPlus * e;


        // Ne pas modifier la suite
        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;

        iter++ ;



        //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        display(cam,I,x,xd) ;


    }

    // sauvegarde des courbes

    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    int a ; cin >> a ;
    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"1pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}




void
tp2DVisualServoingFourPoint()
{
    //-------------------------------------------------------------
    // Mise en oeuvre des courbes

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,8);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------


    //positions initiale (à tester)
    //vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3, vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    //vpHomogeneousMatrix cTw (0.2,0.1,1.3,  0,0,vpMath::rad(5)) ;
    //vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180-10)) ;

    // position finale
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;





    // position des point dans le repere monde Rw
    vpColVector wX[4], cX[4], cdX[4];
    for (int i = 0 ; i < 4 ; i++)
    {
        wX[i].resize(3);
        cX[i].resize(3);
        cdX[i].resize(3);
    } 

    double M = 0.3 ;
    wX[0][0] = -M;      wX[0][1] = -M;      wX[0][2] = 0 ;
    wX[1][0] = M;       wX[1][1] = -M;      wX[1][2] = 0 ;
    wX[2][0] = M;       wX[2][1] =  M;      wX[2][2] = 0 ;
    wX[3][0] = -M;      wX[3][1] =  M;      wX[3][2] = 0 ;

    for(int i=0; i<4; ++i)
    {
        changeFrame(wX[i], cdTw, cdX[i]);
    }

   
    int size = 8;
    vpColVector e(size) ; //

    vpColVector x(size), xd(size) ;

    

    for(int i=0; i<4; ++i)
    {
        vpColVector cx(2), cX(3);
        changeFrame(wX[i], cdTw, cX);
        project(cX, cx);
        xd[2*i] = cx[0];
        xd[2*i+1] = cx[1];
        e[2*i] = e[2*i+1] = 1;
    }

    //initialisation de la position désire des points dans l'image en fonction de cdTw

    

    vpColVector v(6) ;
    double lambda = 0.1 ;
    int iter = 0 ;

    vpMatrix Lx(size,6), LxPlus(6, size);

    while (fabs(e.sumSquare()) > 1e-16)
    {
        vpColVector cx(2);
        // calcul de la position des points dans l'image en fonction de cTw
        for(int i=0; i<4; ++i)
        {
            changeFrame(wX[i], cTw, cX[i]);
            project(cX[i], cx);
            x[2*i] = cx[0];
            x[2*i+1] = cx[1];
        }
        

        // Calcul de la matrice d'interaction
        vpMatrix Lx(size,6) ;
        computeInteractionMatrixMultiPoints(cX, x, Lx);
        //computeInteractionMatrixMultiPoints(cdX, xd, Lx);

        //calcul de l'erreur
        e = x - xd;

        //calcul de la loi de commande
        LxPlus = Lx.pseudoInverse();
        v = -lambda * LxPlus * e;

        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;
        iter++ ;

       //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        display(cam,I,x,xd) ;

        vpTime::wait(30);
    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"4pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}



void
computeError3D(vpColVector const& cdtc, vpRotationVector const& cdrc, vpPoseVector & res)
{
    res = vpPoseVector(cdtc[0], cdtc[1], cdtc[2], cdrc[0], cdrc[1], cdrc[2]);
}

vpPoseVector computeError3D(vpHomogeneousMatrix const& cdTc)
{
    return vpPoseVector(cdTc.getTranslationVector(), cdTc.getThetaUVector());
}

void
computeInteractionMatrix3D(vpRotationMatrix * const& dcRc, int size, vpMatrix & res)
{
    res = vpMatrix(0, 0);
    for(int i=0; i<size; ++i)
    {
        vpMatrix Lx(6, 6);
        Lx = 0;
        for(int k=0; k<3; ++k) for(int l=0; l<3; ++l) Lx[k][l] = dcRc[i][k][l];
        Lx[3][3] = Lx[4][4] = Lx[5][5] = 1;
        res.stack(Lx);
    }
}


void
tp3DVisualServoing()
{
    vpTRACE("begin" ) ;

    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,6);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "Camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);



    //Definition de la scene

    //vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3, vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    //vpHomogeneousMatrix cTw (0.2,0.1,1.3,  0,0,vpMath::rad(5)) ;
    //vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180-10));


    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;

    vpHomogeneousMatrix cdTc = cdTw * cTw.inverse();

    vpPoseVector e = computeError3D(cdTc);
    

    vpColVector v(6) ;
    double lambda = 0.1;
    int iter = 0 ;
    while (fabs(vpColVector(e).sumSquare()) > 1e-6)
    {

        cdTc = cdTw * cTw.inverse();
        e = computeError3D(cdTc);
        //No need for a matrix
        v = vpColVector(e) * (-lambda);


        // Mis à jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;

        iter++ ;

    //mis a jour de courbes
        vpPoseVector crw(cTw) ;
        plot.plot(0,0,iter, vpColVector(e).sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter,crw) ;


    }

// sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    int a ; cin >> a ;

}


void
tp2DVisualServoingFourPointMvt()
{
    vpPlot plot(4, 700, 700, 100, 200, "Curves...");


    char title[40];
    strncpy( title, "||e||", 40 );
    plot.setTitle(0,title);
    plot.initGraph(0,1);

    strncpy( title, "x-xd", 40 );
    plot.setTitle(1, title);
    plot.initGraph(1,8);

    strncpy( title, "camera velocity", 40 );
    plot.setTitle(2, title);
    plot.initGraph(2,6);


    strncpy( title, "camera position", 40 );
    plot.setTitle(3, title);
    plot.initGraph(3,6);

    //-------------------------------------------------------------
    // Affichage des images
    vpImage<unsigned char> I(400,600) ;
    vpDisplayX d ;
    d.init(I) ;
    vpDisplay::display(I);
    vpCameraParameters cam(400,400,300,200) ;

    //-------------------------------------------------------------


    //positions initiale (à tester)
    vpHomogeneousMatrix cTw (-0.2, -0.1, 1.3, vpMath::rad(10), vpMath::rad(20), vpMath::rad(30) ) ;
    //vpHomogeneousMatrix cTw (0.2,0.1,1.3,  0,0,vpMath::rad(5)) ;
    //vpHomogeneousMatrix cTw (0,0,1,  0,0,vpMath::rad(45)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(90)) ;
    //vpHomogeneousMatrix cTw (0, 0, 1,  0, 0, vpMath::rad(180)) ;

    // position finale
    vpHomogeneousMatrix cdTw (0,0,1,  0,0,0) ;





    // position des point dans le repere monde Rw
    vpColVector wX[4], cX[4], cdX[4];
    for (int i = 0 ; i < 4 ; i++)
    {
        wX[i].resize(3);
        cX[i].resize(3);
        cdX[i].resize(3);
    } 

    double M = 0.3 ;
    wX[0][0] = -M;      wX[0][1] = -M;      wX[0][2] = 0 ;
    wX[1][0] = M;       wX[1][1] = -M;      wX[1][2] = 0 ;
    wX[2][0] = M;       wX[2][1] =  M;      wX[2][2] = 0 ;
    wX[3][0] = -M;      wX[3][1] =  M;      wX[3][2] = 0 ;

    for(int i=0; i<4; ++i)
    {
        changeFrame(wX[i], cdTw, cdX[i]);
    }

   
    int size = 8;
    vpColVector e(size) ; //

    vpColVector x(size), xd(size) ;

    

    for(int i=0; i<4; ++i)
    {
        vpColVector cx(2), cX(3);
        changeFrame(wX[i], cdTw, cX);
        project(cX, cx);
        xd[2*i] = cx[0];
        xd[2*i+1] = cx[1];
        e[2*i] = e[2*i+1] = 1;
    }

    //initialisation de la position désire des points dans l'image en fonction de cdTw

    

    vpColVector v(size) ;
    double lambda = 0.001 ;
    int iter = 0 ;

    vpMatrix Lx(size,6), LxPlus(6, size);

    double t = 0;
    while (fabs(e.sumSquare()) > 1e-16)
    {
        t += 1;
        // les points sont animés d'un mouvement de 1cm/s en x dans Rw
        for (int i = 0 ; i < 4 ; i++)
        {
            wX[i][0] += 0.001  * (t > 17000);
            wX[i][1] += 2 * sin(t / 120) * 0;
        }

        vpColVector cx(2);
        // calcul de la position des points dans l'image en fonction de cTw
        for(int i=0; i<4; ++i)
        {
            changeFrame(wX[i], cTw, cX[i]);
            project(cX[i], cx);
            x[2*i] = cx[0];
        }
        // Calcul de la matrice d'interaction
        vpMatrix Lx(size,6) ;
        computeInteractionMatrixMultiPoints(cX, x, Lx);
        //computeInteractionMatrixMultiPoints(cdX, xd, Lx);

        //calcul de l'erreur
        e = x - xd;

        //calcul de la loi de commande
        LxPlus = Lx.pseudoInverse();
        v = -lambda * LxPlus * e;

        //mise a jour de la position de la camera
        cTw = vpExponentialMap::direct(v).inverse()* cTw ;

        cout << "iter "<< iter <<" : "<< e.t() << endl ;
        iter++ ;

       //mise a jour des courbes
        vpPoseVector ctw(cTw) ;
        plot.plot(0,0,iter, e.sumSquare()) ;
        plot.plot(1,iter, e) ;
        plot.plot(2,iter, v) ;
        plot.plot(3,iter, ctw) ;
        //mise a jour de l'image
        display(cam,I,x,xd) ;

        vpTime::wait(0);
    }

    // sauvegarde des courbes
    plot.saveData(0,"e.txt","#");
    plot.saveData(1,"error.txt","#");
    plot.saveData(2,"v.txt","#");
    plot.saveData(3,"p.txt","#");

    // sauvegarde de l'image finale
    {
        vpImage<vpRGBa>  Irgb ;
        vpDisplay::getImage(I,Irgb) ;
        vpImageIo::write(Irgb,"4pt.jpg") ;
    }
    cout << "Clicker sur l'image pour terminer" << endl ;
    vpDisplay::getClick(I) ;
}




int main(int argc, char** argv)
{
    
    //tp2DVisualServoingOnePoint() ;
    //tp2DVisualServoingFourPoint() ;
    tp3DVisualServoing() ;
    //tp2DVisualServoingFourPointMvt() ;

}
