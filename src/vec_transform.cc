/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.cc                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:42 by Zian Huang                               */
/*   Updated: 2022/11/14 23:00:01 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "vec_transform.hh"
#include "inline/cell_operation.hh"

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
