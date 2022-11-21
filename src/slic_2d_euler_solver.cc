/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   slic_2d_euler_solver.cc                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:44 by Zian Huang                               */
/*   Updated: 2022/11/20 19:44:58 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "slic_2d_euler_solver.hh"
#include "inline/primitive_tran.hh"
#include <math.h>

// local function used in updateMaxA
double calcMaxA_eachArr(std::array<double, 4> i_uVec)
{
    // (Euler equation system + ideal gas) only!
    double localGamma = 1.4;
    double maxA = sqrt(primitiveX_Vel(i_uVec) * primitiveX_Vel(i_uVec) + primitiveY_Vel(i_uVec) * primitiveY_Vel(i_uVec)) + sqrt(localGamma * primitivePressure(i_uVec) / i_uVec[0]);

    return maxA;
}

void writeToFileStream(std::ofstream &i_fstream, std::vector<std::vector<std::array<double, 4>>> const &i_inputVec, double i_x0, double i_dx, double i_y0, double i_dy, double i_t, int i_idx)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    // from 2 to nCell + 2 to discard the second-order boundary
    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            i_fstream << i_t << ' ' << i_y0 + (j - 2) * i_dy << ' ' << i_x0 + (i - 2) * i_dx << ' ' << i_inputVec[j][i][i_idx] << std::endl;
        }

        i_fstream << std::endl;
    }
    i_fstream << std::endl;
}

// Definitions #####################################################################################

SLIC_2D_EulerSolver::SLIC_2D_EulerSolver() {}

SLIC_2D_EulerSolver::SLIC_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> i_uVec, int i_nCell_x, int i_nCell_y)
{
    // define user input attributes
    nCell_x_ = i_nCell_x;
    nCell_y_ = i_nCell_y;

    // define the initial conditions, the initial full computational space
    uVec_ = i_uVec;
}

SLIC_2D_EulerSolver::SLIC_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> i_uVec, int i_nCell_x, int i_nCell_y, double i_x0, double i_x1, double i_y0, double i_y1, double i_c, double i_tStop)
{
    // define user input attributes
    nCell_x_ = i_nCell_x;
    nCell_y_ = i_nCell_y;
    x0_ = i_x0;
    x1_ = i_x1;
    y0_ = i_y0;
    y1_ = i_y1;
    c_ = i_c;
    tStop_ = i_tStop;

    // implied values
    dx_ = (i_x1 - i_x0) / i_nCell_x;
    dy_ = (i_y1 - i_y0) / i_nCell_y;

    // define the initial conditions, the initial full computational space
    uVec_ = i_uVec;
}

void SLIC_2D_EulerSolver::setName(std::string i_name)
{
    name_ = i_name;
}

void SLIC_2D_EulerSolver::setBound(double i_x0, double i_x1, double i_y0, double i_y1, double i_tStop)
{
    x0_ = i_x0;
    x1_ = i_x1;
    y0_ = i_y0;
    y1_ = i_y1;
    tStop_ = i_tStop;

    // implied values
    dx_ = (i_x1 - i_x0) / nCell_x_;
    dy_ = (i_y1 - i_y0) / nCell_y_;
}

void SLIC_2D_EulerSolver::setCFL(double i_c)
{
    c_ = i_c;
}

void SLIC_2D_EulerSolver::setRepoDir(std::string repoDir)
{
    repoDir_ = repoDir;
}

void SLIC_2D_EulerSolver::updateMaxA(int i_numIter)
{
    double localMaxA = 0.0;
    double dummieMaxA;

    // loop through the computational domain to find the maximum "aMax"
    for (int j = 2; j < nCell_y_ + 2; ++j)
    {
        for (int i = 2; i < nCell_x_ + 2; ++i)
        {
            dummieMaxA = calcMaxA_eachArr(uVec_[j][i]);

            if (dummieMaxA > localMaxA)
            {
                localMaxA = dummieMaxA;
            }
        }
    }

    // in the first 10 iterations of the numerical scheme, bigger "aMax"
    // this results in a smaller time step, "dt", to deal with sharp discontinuities
    if (i_numIter < 10)
    {
        aMax_ = localMaxA * 3;
    }
    else
    {
        aMax_ = localMaxA;
    }
}

void SLIC_2D_EulerSolver::updateDt()
{
    // see which is the smaller one, "dx" or "dy"
    double smallerD;
    if (dx_ < dy_)
    {
        smallerD = dx_;
    }
    else
    {
        smallerD = dy_;
    }

    // use the smaller delta space distance to determine the stable "dt"
    // dt = c*dx/aMax
    dt_ = c_ * smallerD / aMax_;
}

void SLIC_2D_EulerSolver::updateBoundary_transmissive()
{
    // transmissive boundary
    // updating boundaries in the x-axis
    for (int i = 2; i < nCell_y_ + 2; ++i)
    {
        uVec_[i][0] = uVec_[i][3];
        uVec_[i][1] = uVec_[i][2];
        uVec_[i][nCell_x_ + 3] = uVec_[i][nCell_x_];
        uVec_[i][nCell_x_ + 2] = uVec_[i][nCell_x_ + 1];
    }
    // updating boundaries in the y-axis
    for (int i = 2; i < nCell_x_ + 2; ++i)
    {
        uVec_[0][i] = uVec_[3][i];
        uVec_[1][i] = uVec_[2][i];
        uVec_[nCell_y_ + 3][i] = uVec_[nCell_y_][i];
        uVec_[nCell_y_ + 2][i] = uVec_[nCell_y_ + 1][i];
    }
}

void SLIC_2D_EulerSolver::updateBoundary_cylindrical()
{
    // cylindrical coordinate boundary update, reflective in vel_r
    // updating boundaries in the x-axis
    for (int i = 2; i < nCell_y_ + 2; ++i)
    {
        uVec_[i][0] = uVec_[i][3];
        uVec_[i][1] = uVec_[i][2];
        uVec_[i][nCell_x_ + 3] = uVec_[i][nCell_x_];
        uVec_[i][nCell_x_ + 2] = uVec_[i][nCell_x_ + 1];

        uVec_[i][0][1] = -uVec_[i][3][1];
        uVec_[i][1][1] = -uVec_[i][2][1];
        uVec_[i][nCell_x_ + 3][1] = -uVec_[i][nCell_x_][1];
        uVec_[i][nCell_x_ + 2][1] = -uVec_[i][nCell_x_ + 1][1];
    }
    // updating boundaries in the y-axis
    for (int i = 2; i < nCell_x_ + 2; ++i)
    {
        uVec_[0][i] = uVec_[3][i];
        uVec_[1][i] = uVec_[2][i];
        uVec_[nCell_y_ + 3][i] = uVec_[nCell_y_][i];
        uVec_[nCell_y_ + 2][i] = uVec_[nCell_y_ + 1][i];

        uVec_[i][0][1] = -uVec_[i][3][1];
        uVec_[i][1][1] = -uVec_[i][2][1];
        uVec_[i][nCell_x_ + 3][1] = -uVec_[i][nCell_x_][1];
        uVec_[i][nCell_x_ + 2][1] = -uVec_[i][nCell_x_ + 1][1];
    }
}

void SLIC_2D_EulerSolver::slicLeapX()
{
    uVec_ = vecTran.slicVecTran_x(uVec_, dx_, dt_);
}

void SLIC_2D_EulerSolver::slicLeapY()
{
    uVec_ = vecTran.slicVecTran_y(uVec_, dy_, dt_);
}

void SLIC_2D_EulerSolver::cylindricalSourceTermLeap()
{
    uVec_ = vecTran.cylindricalGeometricRK2(uVec_, dx_, dt_);
}

void SLIC_2D_EulerSolver::initiate()
{
    // set the output directory of the numerical results for further printing purpose
    // the directory is "{REPO_DIR}/data/{NAME}_rhoResults.dat" in the case of the quantity of interest "rho".
    rhoResults_.open(repoDir_ + (std::string) "data/" + name_ + (std::string) "_rhoResults.dat");
    momentumX_Results_.open(repoDir_ + (std::string) "data/" + name_ + (std::string) "_momentumX_Results.dat");
    momentumY_Results_.open(repoDir_ + (std::string) "data/" + name_ + (std::string) "_momentumY_Results.dat");
    energyResults_.open(repoDir_ + (std::string) "data/" + name_ + (std::string) "_energyResults.dat");

    // writing initial conditions to files, with time=0
    writeToFileStream(rhoResults_, uVec_, x0_, dx_, y0_, dy_, 0, 0);
    writeToFileStream(momentumX_Results_, uVec_, x0_, dx_, y0_, dy_, 0, 1);
    writeToFileStream(momentumY_Results_, uVec_, x0_, dx_, y0_, dy_, 0, 2);
    writeToFileStream(energyResults_, uVec_, x0_, dx_, y0_, dy_, 0, 3);
}

void SLIC_2D_EulerSolver::writeToFiles(double i_time)
{
    // using the private ofstream members, with input time, write to local files
    writeToFileStream(rhoResults_, uVec_, x0_, dx_, y0_, dy_, i_time, 0);
    writeToFileStream(momentumX_Results_, uVec_, x0_, dx_, y0_, dy_, i_time, 1);
    writeToFileStream(momentumY_Results_, uVec_, x0_, dx_, y0_, dy_, i_time, 2);
    writeToFileStream(energyResults_, uVec_, x0_, dx_, y0_, dy_, i_time, 3);
}

void SLIC_2D_EulerSolver::cleanUp()
{
    rhoResults_.close();
    momentumX_Results_.close();
    momentumY_Results_.close();
    energyResults_.close();
}

// simple one line definitions
int SLIC_2D_EulerSolver::nCell_x() { return nCell_x_; }
int SLIC_2D_EulerSolver::nCell_y() { return nCell_y_; }
double SLIC_2D_EulerSolver::x0() { return x0_; }
double SLIC_2D_EulerSolver::x1() { return x1_; }
double SLIC_2D_EulerSolver::dx() { return dx_; }
double SLIC_2D_EulerSolver::y0() { return y0_; }
double SLIC_2D_EulerSolver::y1() { return y1_; }
double SLIC_2D_EulerSolver::dy() { return dy_; }
double SLIC_2D_EulerSolver::tStart() { return tStart_; }
double SLIC_2D_EulerSolver::tStop() { return tStop_; }
double SLIC_2D_EulerSolver::c() { return c_; }
double SLIC_2D_EulerSolver::gamma() { return gamma_; }
double SLIC_2D_EulerSolver::aMax() { return aMax_; }
double SLIC_2D_EulerSolver::dt() { return dt_; }