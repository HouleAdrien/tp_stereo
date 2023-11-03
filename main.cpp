/*
 *  FirstProgram.cpp
 *
 *
 *  Created by Olivier Strauss on 10/10/16.
 *  Copyright 2016 LIRMM. All rights reserved.
 *
 */

#include "CImg.h"
#include "string.h"
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <random>

using namespace cimg_library;

// Calcul de la multiplication C = A.B au sens des matrices.
// C est NlinxNcol, A est NlinxNcom, B est NcomxNcol
// Le pointeurs doivent etre alloues correctement avant l'appel
int MatMult(double *A, double *B, double *C, int Nlin, int Ncom, int Ncol)
{
    int lin, col, k;
    double *ptA, *ptB, *ptC;

    ptC = C;
    for (lin = 0; lin < Nlin; lin++)
    {
        for (col = 0; col < Ncol; col++, ptC++)
        {
            ptA = A + (lin * Ncom);
            ptB = B + col;
            (*ptC) = 0.0;
            for (k = 0; k < Ncom; k++)
            {
                (*ptC) += (*ptA) * (*ptB);
                ptA++;
                ptB += Ncol;
            }
        }
    }
    return 1;
}

// Donne des 2 points d'intersection entre la droite
// de parametre L et le carre de dimension Dx, Dy
int Intersection(double L[3], int Dx, int Dy, int x_inter[2], int y_inter[2])
{
    double a = L[0], b = L[1], c = L[2];
    double x[4], y[4];
    int nb_points_ok = 0, n;
    x[0] = 0;
    x[1] = Dx - 1;
    y[2] = 0;
    y[3] = Dy - 1;

    if (fabs(L[0]) > 1e-16)
    { // droite de la forme x = b'y + c' ;
        b = -L[1] / L[0];
        c = -L[2] / L[0];
        x[2] = b * y[2] + c;
        x[3] = b * y[3] + c;
    }
    else
    {
        x[2] = -Dx;
        x[3] = -Dx;
    }

    if (fabs(L[1]) > 1e-16)
    { // droite de la forme x = b'y + c' ;
        a = -L[0] / L[1];
        c = -L[2] / L[1];
        y[0] = a * x[0] + c;
        y[1] = a * x[1] + c;
    }
    else
    {
        y[0] = -Dy;
        y[1] = -Dy;
    }

    for (n = 0; n < 4; n++)
    {
        if ((x[n] >= 0.0) && (y[n] >= 0.0) && (x[n] <= Dx) && (y[n] <= Dy) && (nb_points_ok < 2))
        {
            x_inter[nb_points_ok] = (int)floor(x[n] + 0.5);
            y_inter[nb_points_ok] = (int)floor(y[n] + 0.5);
            nb_points_ok++;
        }
    }

    if (nb_points_ok == 2)
        return 1;
    else
        return 0;
}

std::vector<int> SelectRandomIndices(int n, int p)
{
    std::vector<int> indices;
    indices.reserve(p);

    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < p; ++i)
    {
        std::uniform_int_distribution<int> dist(0, n - 1);
        int index = dist(gen);

        // Vérifie si l'index est déjà dans le vecteur
        if (std::find(indices.begin(), indices.end(), index) == indices.end())
        {
            indices.push_back(index);
        }
        else
        {
            // Si l'index est déjà sélectionné, réessayer pour en obtenir un autre
            --i;
        }
    }

    return indices;
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    int n = values.size();
    if (n % 2 == 0)
    {
        return 0.5 * (values[n / 2 - 1] + values[n / 2]);
    }
    else
    {
        return values[n / 2];
    }
}

std::pair<double, double> Quartiles(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    int n = values.size();
    int q1 = n / 4;
    int q2 = n / 2;
    int q3 = 3 * n / 4;
    return std::make_pair(values[q1], values[q3]);
}

void ComputeEpipolarLine(double *F, int x, int y, int &x1, int &y1, int &x2, int &y2, int Dx, int Dy)
{
    // Convert the point to homogeneous coordinates
    double point[3] = {double(x), double(y), 1.0};

    // Compute the epipolar line
    double l[3];
    MatMult(F, point, l, 3, 3, 1);

    // Find intersection of the epipolar line with image boundaries
    int x_inter[2], y_inter[2];
    if (Intersection(l, Dx, Dy, x_inter, y_inter))
    {
        x1 = x_inter[0];
        y1 = y_inter[0];
        x2 = x_inter[1];
        y2 = y_inter[1];
    }
}

CImg<double> computeNormalizationTransform(double xPoints[], double yPoints[], int count)
{

    // Calculate centroid
    double meanX = 0, meanY = 0;
    for (int i = 0; i < count; i++)
    {
        meanX += xPoints[i];
        meanY += yPoints[i];
    }
    meanX /= count;
    meanY /= count;

    // Calculate mean distance from the centroid
    double meanDist = 0;
    for (int i = 0; i < count; i++)
    {
        double dx = xPoints[i] - meanX;
        double dy = yPoints[i] - meanY;
        double dist = sqrt(dx * dx + dy * dy);
        meanDist += dist;
    }
    meanDist /= count;

    // Calculate scaling factor
    double s = sqrt(2) / meanDist;

    // Create transformation matrix
    CImg<double> T(3, 3, 1, 1, 0);
    T(0, 0) = s;
    T(1, 1) = s;
    T(2, 2) = 1;
    T(0, 2) = -s * meanX;
    T(1, 2) = -s * meanY;

    return T;
}

int main(int argc, char *argv[])
{
    int nombre_de_points = 8, n = 0;
    double xd[nombre_de_points], yd[nombre_de_points], xg[nombre_de_points], yg[nombre_de_points];

    for (int i = 0; i < nombre_de_points; i++)
    {
        xd[i] = yd[i] = xg[i] = yg[i] = -1.0;
    }

    char droite_gauche = 'd';

    if (argc < 3)
    {
        printf("\nCe programme requiert deux images en argument: une pour la droite et une pour la gauche\n");
        printf("Utilisation: ./Stereo ImageDroite.tif ImageGauche.tif\n");
        return 0;
    }

    // Chargement des images
    CImg<unsigned char> imageD(argv[1]), imageG(argv[2]);
    const unsigned char red[] = {255, 0, 0}, green[] = {0, 255, 0}, blue[] = {0, 0, 255};
    CImgDisplay Droite_disp(imageD, "Image Droite"), Gauche_disp(imageG, "Image Gauche");

    while (!Droite_disp.is_closed() && !Gauche_disp.is_closed() && n < nombre_de_points)
    {
        switch (droite_gauche)
        {
        case 'd':
            Droite_disp.set_title("Cliquez ici");
            Gauche_disp.set_title("Image gauche");
            Droite_disp.wait();
            if (Droite_disp.button() && Droite_disp.mouse_y() >= 0)
            {
                xd[n] = Droite_disp.mouse_x();
                yd[n] = Droite_disp.mouse_y();
                imageD.draw_circle(xd[n], yd[n], 3, red).display(Droite_disp);
                droite_gauche = 'g';
            }
            break;
        case 'g':
            Gauche_disp.set_title("Cliquez ici");
            Droite_disp.set_title("Image Droite");
            Gauche_disp.wait();
            if (Gauche_disp.button() && Gauche_disp.mouse_y() >= 0)
            {
                xg[n] = Gauche_disp.mouse_x();
                yg[n] = Gauche_disp.mouse_y();
                imageG.draw_circle(xg[n], yg[n], 3, blue).display(Gauche_disp);
                droite_gauche = 'd';
                n++; // Increment n only after both points have been chosen
            }
            break;
        }
    }

    // Affichage des points en vert après sélection
    for (int i = 0; i < nombre_de_points; i++)
    {
        imageD.draw_circle(xd[i], yd[i], 3, green).display(Droite_disp);
        imageG.draw_circle(xg[i], yg[i], 3, green).display(Gauche_disp);
    }

    // Construisez la matrice A à partir des paires de points correspondants
    CImg<double> matrice_A(nombre_de_points, 9, 1, 1, 0);

    // Compute normalization matrices for both sets of points
    CImg<double> Tx = computeNormalizationTransform(xd, yd, nombre_de_points);
    CImg<double> Ty = computeNormalizationTransform(xg, yg, nombre_de_points);

    double xd_normalized[nombre_de_points];
    double yd_normalized[nombre_de_points];
    double xg_normalized[nombre_de_points];
    double yg_normalized[nombre_de_points];

    // Normalize points
    for (int i = 0; i < nombre_de_points; i++)
    {
        double xw = Tx(0, 0) * xd[i] + Tx(0, 1) * yd[i] + Tx(0, 2);
        double yw = Tx(1, 0) * xd[i] + Tx(1, 1) * yd[i] + Tx(1, 2);
        double w = Tx(2, 0) * xd[i] + Tx(2, 1) * yd[i] + Tx(2, 2);

        xd_normalized[i] = xw / w;
        yd_normalized[i] = yw / w;

        double xgw = Ty(0, 0) * xg[i] + Ty(0, 1) * yg[i] + Ty(0, 2);
        double ygw = Ty(1, 0) * xg[i] + Ty(1, 1) * yg[i] + Ty(1, 2);
        double wg = Ty(2, 0) * xg[i] + Ty(2, 1) * yg[i] + Ty(2, 2);

        xg_normalized[i] = xgw / wg;
        yg_normalized[i] = ygw / wg;
    }

    for (int i = 0; i < nombre_de_points; i++)
    {
        matrice_A(i, 0) = xg_normalized[i] * xd_normalized[i]; // x'x
        matrice_A(i, 1) = xg_normalized[i] * yd_normalized[i]; // x'y
        matrice_A(i, 2) = xg_normalized[i];                    // x'
        matrice_A(i, 3) = yg_normalized[i] * xd_normalized[i]; // y'x
        matrice_A(i, 4) = yg_normalized[i] * yd_normalized[i]; // y'y
        matrice_A(i, 5) = yg_normalized[i];                    // y'
        matrice_A(i, 6) = xd_normalized[i];                    // x
        matrice_A(i, 7) = yd_normalized[i];                    // y
        matrice_A(i, 8) = 1;
    }

    matrice_A = matrice_A.transpose();
    CImg<double> U, S, V;
    matrice_A.SVD(U, S, V);

    matrice_A = matrice_A.transpose();
    // Transpose V to get Vt
    CImg<double> Vt = V.get_transpose();

    CImg<double> f = Vt.get_shared_row(Vt.height() - 1);

    // La matrice fondamentale F est construite à partir de f
    CImg<double> F(3, 3, 1, 1, 0);
    F(0, 0) = f(0);
    F(0, 1) = f(1);
    F(0, 2) = f(2);
    F(1, 0) = f(3);
    F(1, 1) = f(4);
    F(1, 2) = f(5);
    F(2, 0) = f(6);
    F(2, 1) = f(7);
    F(2, 2) = f(8);
    F = F.transpose();
    CImgList<double> SVD_result = F.get_SVD();
    U = SVD_result[0];
    S = SVD_result[1];
    V = SVD_result[2];
    F = F.transpose();

    std::cout << "Singular values before denormalizaion:" << std::endl;
    for (int i = 0; i < S.height(); ++i)
    {
        if (std::abs(S(0, i)) > 1e-5)
            std::cout << S(0, i) << " ";
    }
    std::cout << std::endl;

    // Recompute F with the corrected singular values
    CImg<double> F_temp(3, 3, 1, 1, 0); // Temporary matrix to store the result
    MatMult(U.data(), S.data(), F_temp.data(), 3, 3, 3);

    CImg<double> F_temp3(3, 3, 1, 1, 0); // Temporary matrix to store the result
    MatMult(F_temp.data(), V.get_transpose().data(), F_temp3.data(), 3, 3, 3);

    CImg<double> F_corrected = F_temp3;

    // Assuming Tx and Ty are the normalization transformations for the two sets of points
    CImg<double> Tx_transpose = Tx.transpose();
    CImg<double> Ty_transpose = Ty.transpose();

    CImg<double> F_temp2(3, 3, 1, 1, 0); // Temporary matrix to store the result
    MatMult(Tx_transpose.data(), F_corrected.data(), F_temp2.data(), 3, 3, 3);

    // Perform the multiplication for the final denormalization step (F_temp * Ty)
    CImg<double> F_denormalized(3, 3, 1, 1, 0); // This will store the final denormalized fundamental matrix
    MatMult(F_temp.data(), Ty.data(), F_denormalized.data(), 3, 3, 3);

    // Display the denormalized F matrix
    CImg<double>::iterator it = F_denormalized.begin();
    int NlinF_denorm = F_denormalized.height();
    int NcolF_denorm = F_denormalized.width();

    std::cout << "F_denormalized:" << std::endl;
    for (int i = 0; i < NlinF_denorm; ++i)
    {
        for (int j = 0; j < NcolF_denorm; ++j)
        {
            std::cout << F_denormalized(i, j) << " ";
        }
        std::cout << std::endl;
    }
    // Compute SVD of the F matrix
    F_denormalized = F_denormalized.get_transpose();
    SVD_result = F_denormalized.get_SVD();
    U = SVD_result[0];
    S = SVD_result[1];
    V = SVD_result[2];
    F_denormalized = F_denormalized.get_transpose();

    // Print singular values before modification
    printf("\nSingular values before enforcing rank 2 condition:\n");
    for (int i = 0; i < S.height(); i++)
    {
        printf("%g ", S(i, i)); // Assuming S is a square diagonal matrix
    }
    printf("\n");

    // Find the two largest singular values and zero out the rest to enforce rank 2
    double first_max = 0.0;
    double second_max = 0.0;
    int first_index = -1;
    int second_index = -1;

    // Find indices of two largest singular values
    for (int i = 0; i < S.height(); i++)
    {
        if (S(i, i) > first_max)
        {
            second_max = first_max;
            second_index = first_index;
            first_max = S(i, i);
            first_index = i;
        }
        else if (S(i, i) > second_max)
        {
            second_max = S(i, i);
            second_index = i;
        }
    }

    // Set all singular values to zero except the two largest
    for (int i = 0; i < S.height(); i++)
    {
        if (i != first_index && i != second_index)
        {
            S(i, i) = 0.0;
        }
    }

    // Reconstruct the fundamental matrix with the modified singular values
    MatMult(U.data(), S.data(), F_temp.data(), 3, 3, 3);
    MatMult(F_temp.data(), V.get_transpose().data(), F_denormalized.data(), 3, 3, 3);

    // Print singular values after modification
    printf("\nSingular values after enforcing rank 2 condition:\n");
    for (int i = 0; i < S.height(); i++)
    {
        printf("%g ", S(i, i));
    }
    printf("\n");

    // Comptez les valeurs singulières non nulles
    int nonZeroSingularValues = 0;
    for (int i = 0; i < S.height(); i++)
    {
        if (std::abs(S(0, i)) > 1e-7)
        { // Vous pouvez ajuster la tolérance selon vos besoins
            nonZeroSingularValues++;
        }
    }
    if (nonZeroSingularValues == 2)
    {
        printf("La matrice fondamentale est correcte. Elle a exactement deux valeurs singulières non nulles.\n");
     
    }
    else
    {
        printf("La matrice fondamentale n'est pas correcte. Elle devrait avoir exactement deux valeurs singulières non nulles.\n");
    }

    while (!Droite_disp.is_closed() && !Gauche_disp.is_closed()) {

        Droite_disp.set_title("Cliquez ici");
        Gauche_disp.set_title("Image Gauche");
        Droite_disp.wait();
        if (Droite_disp.button() && Droite_disp.mouse_y() >= 0) {
            int x = Droite_disp.mouse_x();
            int y = Droite_disp.mouse_y();
            
            // Affiche les coordonnées du point cliqué
            std::cout << "Point cliqué sur l'image droite: (" << x << ", " << y << ")\n";

            imageD.draw_circle(x, y, 3, red).display(Droite_disp);

            CImg<double> pointDroite(3, 1);
            pointDroite(0, 0) = x;
            pointDroite(1, 0) = y;
            pointDroite(2, 0) = 1.0;
            CImg<double> ligneGauche;  

            int Nlin = F.height();    
            int Ncom = F.width();     
            int Ncol = 1;             

            ligneGauche.assign(Nlin, Ncol);

            MatMult(F.data(), pointDroite.data(), ligneGauche.data(), Nlin, Ncom, Ncol);
            
            double A = ligneGauche(0, 0);
            double B = ligneGauche(1, 0);
            double C = ligneGauche(2, 0);

            // Affiche les valeurs A, B et C pour la ligne
            std::cout << "Valeurs pour la ligne: A=" << A << ", B=" << B << ", C=" << C << "\n";
            
            int x_inter[2], y_inter[2];
            double L[3] = {A, B, C};
            int Dx = imageG.width();
            int Dy = imageG.height();

            if(Intersection(L, Dx, Dy, x_inter, y_inter)) {
                // Affiche les coordonnées d'intersection
                std::cout << "Points d'intersection sur l'image gauche: (" << x_inter[0] << ", " << y_inter[0] << ") et (" << x_inter[1] << ", " << y_inter[1] << ")\n";

                imageG.draw_line(x_inter[0], y_inter[0], x_inter[1], y_inter[1], red).display(Gauche_disp);
            } else {
                // Affiche un message si aucune intersection n'est trouvée
                std::cout << "Aucune intersection trouvée pour ces coordonnées.\n";
            }
        }
    }


    // Attente de la fermeture d'une des images pour arrÃªter le programme
    while (!Droite_disp.is_closed() && !Gauche_disp.is_closed())
        ;

    return 0;
}
