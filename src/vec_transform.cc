/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.cc                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:42 by Zian Huang                               */
/*   Updated: 2022/11/21 11:15:48 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "vec_transform.hh"
#include "inline/cell_operation.hh"
#include "inline/primitive_tran.hh"

// source term in cylicdrical symmetry
std::array<double, 4> cylindricalSource(std::array<double, 4> i_uArr, double i_dr, int i_idx)
{
    std::array<double, 4> dummieArr;

    // a negative sign
    dummieArr[0] = -i_uArr[0] * primitiveX_Vel(i_uArr) / (i_dr * (i_idx + 0.5));
    dummieArr[1] = -i_uArr[0] * primitiveX_Vel(i_uArr) * primitiveX_Vel(i_uArr) / (i_dr * (i_idx + 0.5));
    dummieArr[2] = -i_uArr[0] * primitiveY_Vel(i_uArr) * primitiveX_Vel(i_uArr) / (i_dr * (i_idx + 0.5));
    dummieArr[3] = -(i_uArr[3] + primitivePressure(i_uArr)) * primitiveX_Vel(i_uArr) / (i_dr * (i_idx + 0.5));

    return dummieArr;
}

// Definitions #####################################################################################

VecTran::VecTran() {}

std::vector<std::vector<std::array<double, 4>>> VecTran::slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dx, double i_dt)
{
    // aquire dimension information
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    // construct a return vector with the same dimension as the input
    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    // update the whole "real" computational domain along the x-axis
    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dx, diffCell(fluxFunc.slicFlux_x(i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_inputU_Vec[iter_y][iter_x + 2], i_dx, i_dt), fluxFunc.slicFlux_x(i_inputU_Vec[iter_y][iter_x - 2], i_inputU_Vec[iter_y][iter_x - 1], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y][iter_x + 1], i_dx, i_dt))));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &i_inputU_Vec, double i_dy, double i_dt)
{
    // aquire dimension information
    int xVecLen = i_inputU_Vec[0].size();
    int yVecLen = i_inputU_Vec.size();

    // construct a return vector with the same dimension as the input
    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }

    // update the whole "real" computational domain along the y-axis
    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            toBeReturnVec[iter_y][iter_x] = diffCell(i_inputU_Vec[iter_y][iter_x], scalingCell(i_dt / i_dy, diffCell(fluxFunc.slicFlux_y(i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_inputU_Vec[iter_y + 2][iter_x], i_dy, i_dt), fluxFunc.slicFlux_y(i_inputU_Vec[iter_y - 2][iter_x], i_inputU_Vec[iter_y - 1][iter_x], i_inputU_Vec[iter_y][iter_x], i_inputU_Vec[iter_y + 1][iter_x], i_dy, i_dt))));
        }
    }

    return toBeReturnVec;
}

std::vector<std::vector<std::array<double, 4>>> VecTran::cylindricalGeometricRK2(const std::vector<std::vector<std::array<double, 4>>> &i_uVec, double i_dx, double i_dt)
{
    // adapted from the above vectran methods:
    // ------------------------------
    // aquire dimension information
    int xVecLen = i_uVec[0].size();
    int yVecLen = i_uVec.size();

    // construct a return vector with the same dimension as the input
    std::vector<std::vector<std::array<double, 4>>> toBeReturnVec;
    toBeReturnVec.resize(yVecLen);
    for (int i = 0; i < yVecLen; ++i)
    {
        toBeReturnVec[i].resize(xVecLen);
    }
    // ------------------------------

    // 2nd order Runge-Kutta method
    for (int iter_y = 2; iter_y < yVecLen - 2; ++iter_y)
    {
        for (int iter_x = 2; iter_x < xVecLen - 2; ++iter_x)
        {
            std::array<double, 4> k1;
            std::array<double, 4> k2;

            // `iter_x - 2` was used here so `idx` starts from 0 and ends at xVecLen
            k1 = scalingCell(i_dt, cylindricalSource(i_uVec[iter_y][iter_x], i_dx, iter_x - 2));
            k2 = scalingCell(i_dt, cylindricalSource(sumCell(i_uVec[iter_y][iter_x], k1), i_dx, iter_x - 2));

            toBeReturnVec[iter_y][iter_x] = sumCell(i_uVec[iter_y][iter_x], scalingCell(0.5, sumCell(k1, k2)));
        }
    }

    return toBeReturnVec;
}
