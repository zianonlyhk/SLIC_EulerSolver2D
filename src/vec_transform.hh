/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.hh                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:41 by Zian Huang                               */
/*   Updated: 2022/11/20 19:38:34 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef VEC_TRANSFORM_HH
#define VEC_TRANSFORM_HH

#include <vector>
#include <array>
#include "flux_func.hh"

// this class is considered more as a collection of functions/methods, than a C++ object
class VecTran
{
public:
    // constructor
    VecTran();
    FluxFVM fluxFunc;

    // force introduced by quantity gradient
    // pass by const reference to reduce memory pressure
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);

    // force as the source term
    std::vector<std::vector<std::array<double, 4>>> cylindricalGeometricRK2(const std::vector<std::vector<std::array<double, 4>>> &uVec, double dx, double dt);
};

#endif