/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   vec_transform.hh                                  Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:41 by Zian Huang                               */
/*   Updated: 2022/11/14 14:35:55 by Zian Huang                               */
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

    // pass by const reference to reduce memory pressure
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_x(const std::vector<std::vector<std::array<double, 4>>> &, double, double);
    std::vector<std::vector<std::array<double, 4>>> slicVecTran_y(const std::vector<std::vector<std::array<double, 4>>> &, double, double);
};

#endif