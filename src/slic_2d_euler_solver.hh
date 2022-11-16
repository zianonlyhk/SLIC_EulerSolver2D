/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   slic_2d_euler_solver.hh                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:43 by Zian Huang                               */
/*   Updated: 2022/11/14 21:25:17 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef SLIC_2D_EULER_SOLVER_HH
#define SLIC_2D_EULER_SOLVER_HH

#include "vec_transform.hh"

class SLIC_2D_EulerSolver
{
public:
    // CONSTRUCTORS
    SLIC_2D_EulerSolver();
    SLIC_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> i_uVec, int i_nCell_x, int i_nCell_y, double i_x0, double i_x1, double i_y0, double i_y1, double i_c, double i_tStop);

    // TOOLBOX
    // we need to use external a toolbox called "VecTran" - vector transform
    VecTran vecTran;
    // loosely packed tools
    void updateMaxA();
    void updateDt();
    void updateBoundary();
    void slicLeapX();
    void slicLeapY();

    // CONSTANT ATTRIBUTES
    // constant attributes after class construction
    int nCell_x_;
    int nCell_y_;
    double x0_;
    double x1_;
    double dx_;
    double y0_;
    double y1_;
    double dy_;
    double tStart_ = 0.0;
    double tStop_;
    double c_;
    double gamma_ = 1.4;

    // DYNAMIC ATTRIBUTES
    // dynamically updated attributes
    double aMax_;
    double dt_;
    // "full" computional domain
    std::vector<std::vector<std::array<double, 4>>> uVec_;
};

#endif