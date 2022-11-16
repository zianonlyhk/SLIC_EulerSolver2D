/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   slic_2d_euler_solver.cc                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:44 by Zian Huang                               */
/*   Updated: 2022/11/16 10:00:29 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "slic_2d_euler_solver.hh"
#include "inline/primitive_tran.hh"
#include <math.h>

// local function used in updateMaxA
double calcMaxA_eachArr(std::array<double, 4> i_uVec)
{
    double localGamma = 1.4;
    double maxA = sqrt(primitiveX_Vel(i_uVec) * primitiveX_Vel(i_uVec) + primitiveY_Vel(i_uVec) * primitiveY_Vel(i_uVec)) + sqrt(localGamma * primitivePressure(i_uVec) / i_uVec[0]);

    return maxA;
}

// Definitions #####################################################################################

SLIC_2D_EulerSolver::SLIC_2D_EulerSolver() {}

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

    updateBoundary();
    updateMaxA();
    updateDt();
}

void SLIC_2D_EulerSolver::updateMaxA()
{
    double localMaxA = 0.0;
    double dummieMaxA;

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

    aMax_ = localMaxA;
}

void SLIC_2D_EulerSolver::updateDt()
{
    double smallerD;
    if (dx_ < dy_)
    {
        smallerD = dx_;
    }
    else
    {
        smallerD = dy_;
    }

    // dt = c*dx/aMax
    dt_ = c_ * smallerD / aMax_;
}

void SLIC_2D_EulerSolver::updateBoundary()
{
    // updating boundaries in the x-axis
    for (int i = 2; i < nCell_y_ + 2; ++i)
    {
        uVec_[i][0] = uVec_[i][2];
        uVec_[i][1] = uVec_[i][2];
        uVec_[i][nCell_x_ + 3] = uVec_[i][nCell_x_ + 1];
        uVec_[i][nCell_x_ + 2] = uVec_[i][nCell_x_ + 1];
    }
    // updating boundaries in the y-axis
    for (int i = 2; i < nCell_x_ + 2; ++i)
    {
        uVec_[0][i] = uVec_[2][i];
        uVec_[1][i] = uVec_[2][i];
        uVec_[nCell_y_ + 3][i] = uVec_[nCell_y_ + 1][i];
        uVec_[nCell_y_ + 2][i] = uVec_[nCell_y_ + 1][i];
    }

    // left upper corner
    uVec_[0][0] = uVec_[2][2];
    uVec_[0][1] = uVec_[2][2];
    uVec_[1][0] = uVec_[2][2];
    uVec_[1][1] = uVec_[2][2];
    // right upper corner
    uVec_[0][nCell_x_ + 3] = uVec_[2][nCell_x_ + 1];
    uVec_[0][nCell_x_ + 2] = uVec_[2][nCell_x_ + 1];
    uVec_[1][nCell_x_ + 2] = uVec_[2][nCell_x_ + 1];
    uVec_[1][nCell_x_ + 3] = uVec_[2][nCell_x_ + 1];
    // left lower corner
    uVec_[nCell_y_ + 2][0] = uVec_[nCell_y_ + 1][2];
    uVec_[nCell_y_ + 2][1] = uVec_[nCell_y_ + 1][2];
    uVec_[nCell_y_ + 3][0] = uVec_[nCell_y_ + 1][2];
    uVec_[nCell_y_ + 3][1] = uVec_[nCell_y_ + 1][2];
    // right lower corner
    uVec_[nCell_y_ + 3][nCell_x_ + 3] = uVec_[nCell_y_ + 1][nCell_x_ + 1];
    uVec_[nCell_y_ + 2][nCell_x_ + 2] = uVec_[nCell_y_ + 1][nCell_x_ + 1];
    uVec_[nCell_y_ + 2][nCell_x_ + 3] = uVec_[nCell_y_ + 1][nCell_x_ + 1];
    uVec_[nCell_y_ + 3][nCell_x_ + 2] = uVec_[nCell_y_ + 1][nCell_x_ + 1];
}

void SLIC_2D_EulerSolver::slicLeapX()
{
    uVec_ = vecTran.slicVecTran_x(uVec_, dx_, dt_);
}

void SLIC_2D_EulerSolver::slicLeapY()
{
    uVec_ = vecTran.slicVecTran_y(uVec_, dy_, dt_);
}
