/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   slic_2d_euler_solver.hh                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:43 by Zian Huang                               */
/*   Updated: 2022/11/14 17:16:22 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef SLIC_2D_EULER_SOLVER_HH
#define SLIC_2D_EULER_SOLVER_HH

#include "vec_transform.hh"

class SLIC_2D_EulerSolver
{
public:
    // constructors
    SLIC_2D_EulerSolver();
    SLIC_2D_EulerSolver(int, int, double, double, double, double, double, double);

    // we need to use external vector transform class methods
    VecTran vecTran;

    // class methods
    double calcTimeStep();
    void updateBoundary();

    // class members

    // user-defined attributes at construction
    int nCell_x_;
    int nCell_y_;
    double x0_;
    double x1_;
    double y0_;
    double y1_;
    double c_;
    double tStop_;
    // implied attributes as a consequence of user definition
    double tStart_ = 0.0;
    double dx_;
    double dy_;
    double gamma_ = 1.4;
    // being dynamically updated
    double dt_;

    // "full" computional domain
    std::vector<std::vector<std::array<double, 4>>> uVec_;
};

#endif