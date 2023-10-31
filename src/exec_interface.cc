/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   exec_interface.cc                                 Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:25 by Zian Huang                               */
/*   Updated: 2022/11/21 10:58:33 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "vec_transform.hh"
#include "slic_2d_euler_solver.hh"
#include <fstream>
#include <iostream>

// the following initial value set was copied from Steve's ACM1 practical 2
// these values are all in primitive form: density, xMomentum, yMomentum, pressure
void setInitRiemannValue(std::vector<std::vector<std::array<double, 4>>> &i_inputVec, double i_x0, double i_x1, double i_y0, double i_y1)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    double dx = (i_x1 - i_x0) / nCell_x;
    double dy = (i_x1 - i_x0) / nCell_y;

    double currX;
    double currY;
    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            currX = i_x0 + (i - 2) * dx;
            currY = i_y0 + (j - 2) * dy;

            if (currX < 0.5)
            {
                if (currY < 0.5)
                {
                    i_inputVec[j][i] = (std::array<double, 4>){1, -0.75, 0.5, 1};
                }
                else
                {
                    i_inputVec[j][i] = (std::array<double, 4>){2, 0.75, 0.5, 1};
                }
            }
            else
            {
                if (currY < 0.5)
                {
                    i_inputVec[j][i] = (std::array<double, 4>){3, -0.75, -0.5, 1};
                }
                else
                {
                    i_inputVec[j][i] = (std::array<double, 4>){1, 0.75, -0.5, 1};
                }
            }
        }
    }
}

// since we need to work in conservation form as we adapted finite volume method
// the following local function will transform our initial primitive conditions into conservative form
void conservativeFormTransform(std::vector<std::vector<std::array<double, 4>>> &i_inputVec)
{
    double localGamma = 1.4;
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    double velX_Copy;
    double velY_Copy;
    double pressureCopy;

    for (int j = 0; j < nCell_y + 4; ++j)
    {
        for (int i = 0; i < nCell_x + 4; ++i)
        {
            velX_Copy = i_inputVec[j][i][1];
            velY_Copy = i_inputVec[j][i][2];
            pressureCopy = i_inputVec[j][i][3];
            i_inputVec[j][i][1] = velX_Copy * i_inputVec[j][i][0];
            i_inputVec[j][i][2] = velY_Copy * i_inputVec[j][i][0];
            i_inputVec[j][i][3] = pressureCopy / (localGamma - 1) + 0.5 * i_inputVec[j][i][0] * (velX_Copy * velX_Copy + velY_Copy * velY_Copy);
        }
    }
}

int main()
{
    // some constant declarations and definitions here
    int nCells_x;
    int nCells_y;
    double x0;
    double x1;
    double y0;
    double y1;
    double c;
    double tStop;
    nCells_x = 100;
    nCells_y = 100;
    x0 = 0.0;
    x1 = 1.0;
    y0 = 0.0;
    y1 = 1.0;
    c = 0.8;
    tStop = 1.6;

    // prepare the vector of vectors of arrays
    std::vector<std::vector<std::array<double, 4>>> compDomain;
    compDomain.resize(nCells_y + 4);
    for (int i = 0; i < nCells_y + 4; ++i)
    {
        compDomain[i].resize(nCells_x + 4);
    }

    // define the arrays
    setInitRiemannValue(compDomain, x0, x1, y0, y1);
    // primitive input to conservative variables
    conservativeFormTransform(compDomain);

    // create a class instance
    SLIC_2D_EulerSolver testSolverClass(compDomain, nCells_x, nCells_y);
    // update its private members
    testSolverClass.setBound(x0, x1, y0, y1, tStop);
    testSolverClass.setCFL(c);
    testSolverClass.setName((std::string) "test");
    testSolverClass.setRepoDir((std::string) "/Users/zian/Room214N/github/SLIC_EulerSolver2D/");

    testSolverClass.updateBoundary_cylindrical();

    testSolverClass.initiate();

    // starting the numerical scheme #############################################################
    // preparing the while loop
    double t = 0.0;
    int numIter = 0;
    do
    {
        ++numIter;

        // preparation for time leap
        testSolverClass.updateMaxA(numIter);
        testSolverClass.updateDt();

        // time leaping forward
        t += testSolverClass.dt();

        // matrix transform
        testSolverClass.slicLeapX();
        testSolverClass.updateBoundary_cylindrical();
        testSolverClass.slicLeapY();
        testSolverClass.updateBoundary_cylindrical();

        testSolverClass.cylindricalSourceTermLeap();
        testSolverClass.updateBoundary_cylindrical();

        // writing current time frame to the output file
        if (numIter % 20 == 0)
        {
            testSolverClass.writeToFiles(t);
            std::cout << t << " / " << testSolverClass.tStop() << std::endl;
        }

    } while (t < testSolverClass.tStop());

    // closing the fstream
    testSolverClass.cleanUp();

    return 0;
}