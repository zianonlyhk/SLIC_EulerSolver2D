/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   slic_2d_euler_solver.hh                           Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:43 by Zian Huang                               */
/*   Updated: 2022/11/20 19:42:40 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#ifndef SLIC_2D_EULER_SOLVER_HH
#define SLIC_2D_EULER_SOLVER_HH

#include "vec_transform.hh"
#include <string>
#include <fstream>

class SLIC_2D_EulerSolver
{
public:
    // CONSTRUCTORS
    // default constructor
    SLIC_2D_EulerSolver();
    // normal constructor, best for readability
    SLIC_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> uVec, int nCell_x, int nCell_y);
    // full constructor, all in one line
    SLIC_2D_EulerSolver(std::vector<std::vector<std::array<double, 4>>> uVec, int nCell_x, int nCell_y, double x0, double x1, double y0, double y1, double c, double tStop);

    // ADDITIONAL CONSTRUCTION
    // defining private members
    void setBound(double x0, double x1, double y0, double y1, double tStop);
    void setCFL(double c);
    void setName(std::string name);
    void setRepoDir(std::string repoDir);

    // NUMERICAL SCHEME TOOLBOX
    // we need to use external a toolbox called "VecTran" - vector transform
    VecTran vecTran;
    // loosely packed tools to update private members
    void updateMaxA(int numIter);
    void updateDt();
    void updateBoundary_transmissive();
    void updateBoundary_cylindrical();
    void slicLeapX();
    void slicLeapY();
    void cylindricalSourceTermLeap();
    // utilities for outputing results
    void initiate();
    void writeToFiles(double time);
    void cleanUp();

    // ACCESSING PRIVATE MEMBERS
    int nCell_x();
    int nCell_y();
    double x0();
    double x1();
    double dx();
    double y0();
    double y1();
    double dy();
    double tStart();
    double tStop();
    double c();
    double gamma();
    double aMax();
    double dt();

private:
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
    // attributes for outputing results to files
    std::string name_;
    std::string repoDir_;
    std::ofstream rhoResults_;
    std::ofstream momentumX_Results_;
    std::ofstream momentumY_Results_;
    std::ofstream energyResults_;

    // DYNAMIC ATTRIBUTES
    // dynamically updated attributes
    double aMax_;
    double dt_;
    // "full" computional domain
    std::vector<std::vector<std::array<double, 4>>> uVec_;
};

#endif