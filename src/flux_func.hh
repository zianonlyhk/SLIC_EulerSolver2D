/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   flux_func.hh                                      Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:31 by Zian Huang                               */
/*   Updated: 2022/11/13 17:37:04 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef FLUX_FUNC_HH
#define FLUX_FUNC_HH

#include <array>

// this class is considered more as a collection of functions/methods, than a C++ object
class FluxFVM
{
public:
    // constructor
    FluxFVM();

    // instantaneous flux based on the state of cell
    std::array<double, 4> conservationFlux_x(std::array<double, 4>);
    std::array<double, 4> conservationFlux_y(std::array<double, 4>);

    // centred scheme flux that requires 2 input, the 2 neighbours
    // since there are 2 neighbours involved, we also need the information of the distance between them
    // they are, "dx" and "dt"
    std::array<double, 4> forceFlux_x(std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> forceFlux_y(std::array<double, 4>, std::array<double, 4>, double, double);

    // SLIC scheme flux with 5 input cell arrays, very long!
    // same as the FORCE scheme, more than 1 cell instance so "dx" and "dt" required
    std::array<double, 4> slicFlux_x(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
    std::array<double, 4> slicFlux_y(std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, std::array<double, 4>, double, double);
};

#endif