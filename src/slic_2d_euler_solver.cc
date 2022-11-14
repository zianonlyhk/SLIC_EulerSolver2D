/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   slic_2d_euler_solver.cc                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:44 by Zian Huang                               */
/*   Updated: 2022/11/14 14:23:07 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "slic_2d_euler_solver.hh"
#include "inline/primitive_tran.hh"
#include <math.h>

SLIC_2D_EulerSolver::SLIC_2D_EulerSolver() {}

SLIC_2D_EulerSolver::SLIC_2D_EulerSolver(int i_nCell_x, int i_nCell_y, double i_x0, double i_x1, double i_y0, double i_y1, double i_c, double i_tStop)
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

    // transform the computational domain from "real" to "full"
    uVec_.resize(nCell_y_ + 4);
    for (int i = 0; i < nCell_y_ + 4; ++i)
    {
        uVec_[i].resize(nCell_x_ + 4);
    }

    // implied values
    dx_ = (i_x1 - i_x0) / i_nCell_x;
    dy_ = (i_y1 - i_y0) / i_nCell_y;
}

double calc_aMax_x(std::array<double, 4> i_uVec)
{
    double localGamma = 1.4;
    return fabs(primitiveX_Vel(i_uVec)) + sqrt(localGamma * primitivePressure(i_uVec) / i_uVec[0]);
}

double calc_aMax_y(std::array<double, 4> i_uVec)
{
    double localGamma = 1.4;
    return fabs(primitiveY_Vel(i_uVec)) + sqrt(localGamma * primitivePressure(i_uVec) / i_uVec[0]);
}

double SLIC_2D_EulerSolver::calcTimeStep()
{
    double local_aMax = 0.0;
    double temp_aMax_x;
    double temp_aMax_y;

    for (int j = 2; j < nCell_y_ + 2; ++j)
    {
        for (int i = 2; i < nCell_x_ + 2; ++i)
        {
            temp_aMax_x = calc_aMax_x(uVec_[j][i]);
            temp_aMax_y = calc_aMax_y(uVec_[j][i]);

            if (temp_aMax_x > local_aMax)
            {
                local_aMax = temp_aMax_x;
            }
            if (temp_aMax_y > local_aMax)
            {
                local_aMax = temp_aMax_y;
            }
        }
    }

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
    return c_ * smallerD / local_aMax;
}

void SLIC_2D_EulerSolver::updateBoundary()
{
    for (int i = 2; i < nCell_y_ + 2; ++i)
    {
        uVec_[i][0] = uVec_[i][3];
        uVec_[i][1] = uVec_[i][2];
        uVec_[i][nCell_x_ + 3] = uVec_[i][nCell_x_];
        uVec_[i][nCell_x_ + 2] = uVec_[i][nCell_x_ + 1];
    }
    for (int i = 2; i < nCell_x_ + 2; ++i)
    {
        uVec_[0][i] = uVec_[3][i];
        uVec_[1][i] = uVec_[2][i];
        uVec_[nCell_y_ + 3][i] = uVec_[nCell_y_][i];
        uVec_[nCell_y_ + 2][i] = uVec_[nCell_y_ + 1][i];
    }

    // left upper corner
    uVec_[0][0] = uVec_[3][3];
    uVec_[0][1] = uVec_[2][3];
    uVec_[1][0] = uVec_[3][2];
    uVec_[1][1] = uVec_[2][2];
    // right upper corner
    uVec_[0][nCell_x_ + 2] = uVec_[3][nCell_x_];
    uVec_[0][nCell_x_ + 3] = uVec_[2][nCell_x_ + 1];
    uVec_[1][nCell_x_ + 2] = uVec_[2][nCell_x_ + 1];
    uVec_[1][nCell_x_ + 3] = uVec_[3][nCell_x_ + 1];
    // left lower corner
    uVec_[nCell_y_ + 2][0] = uVec_[nCell_y_][2];
    uVec_[nCell_y_ + 2][1] = uVec_[nCell_y_ + 1][2];
    uVec_[nCell_y_ + 3][0] = uVec_[nCell_y_][3];
    uVec_[nCell_y_ + 3][1] = uVec_[nCell_y_ + 1][3];
    // right lower corner
    uVec_[nCell_y_ + 3][nCell_x_ + 3] = uVec_[nCell_y_][nCell_x_];
    uVec_[nCell_y_ + 2][nCell_x_ + 2] = uVec_[nCell_y_ + 1][nCell_x_ + 1];
    uVec_[nCell_y_ + 2][nCell_x_ + 3] = uVec_[nCell_y_][nCell_x_ + 1];
    uVec_[nCell_y_ + 3][nCell_x_ + 2] = uVec_[nCell_y_ + 1][nCell_x_];
}
