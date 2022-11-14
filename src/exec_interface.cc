/* ************************************************************************** */
/*                                                                            */
/*                                                                            */
/*   exec_interface.cc                                 Personal Website       */
/*                                                    ##################      */
/*   By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||      */
/*                                                    ##################      */
/*   Created: 2022/11/09 19:14:25 by Zian Huang                               */
/*   Updated: 2022/11/14 15:16:15 by Zian Huang                               */
/*                                                                            */
/* ************************************************************************** */

#include "vec_transform.hh"
#include "slic_2d_euler_solver.hh"
#include <fstream>
#include <iostream>

void DEBUG_printOutVecRho(std::vector<std::vector<std::array<double, 4>>> const &i_inputVec)
{
    int nCellY = i_inputVec.size();
    int nCellX = i_inputVec[0].size();
    for (int j = 0; j < nCellY; ++j)
    {
        for (int i = 0; i < nCellX; ++i)
        {
            std::cout << i_inputVec[j][i][0] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

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
                    i_inputVec[j][i] = (std::array<double, 4>){0.8, 0, 0, 1};
                }
                else
                {
                    i_inputVec[j][i] = (std::array<double, 4>){1, 0.7276, 0, 1};
                }
            }
            else
            {
                if (currY < 0.5)
                {
                    i_inputVec[j][i] = (std::array<double, 4>){1, 0, 0.7276, 1};
                }
                else
                {
                    i_inputVec[j][i] = (std::array<double, 4>){0.5313, 0, 0, 0.4};
                }
            }
        }
    }
}

void conservativeFormTransform(std::vector<std::vector<std::array<double, 4>>> &i_inputVec)
{
    double localGamma = 1.4;
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    for (int j = 0; j < nCell_y + 4; ++j)
    {
        for (int i = 0; i < nCell_x + 4; ++i)
        {
            i_inputVec[j][i][1] = i_inputVec[j][i][1] * i_inputVec[j][i][0];
            i_inputVec[j][i][2] = i_inputVec[j][i][2] * i_inputVec[j][i][0];
            i_inputVec[j][i][3] = i_inputVec[j][i][3] / (localGamma - 1) + (i_inputVec[j][i][0] * i_inputVec[j][i][0] * i_inputVec[j][i][1] * i_inputVec[j][i][1] + i_inputVec[j][i][0] * i_inputVec[j][i][0] * i_inputVec[j][i][2] * i_inputVec[j][i][2]) / 2 / i_inputVec[j][i][0];
        }
    }
}

void writeToFileStream(std::ofstream &i_fstream, std::vector<std::vector<std::array<double, 4>>> const &i_inputVec, double i_x0, double i_dx, double i_y0, double i_dy, double i_t, int i_idx)
{
    int nCell_y = i_inputVec.size() - 4;
    int nCell_x = i_inputVec[0].size() - 4;

    for (int j = 2; j < nCell_y + 2; ++j)
    {
        for (int i = 2; i < nCell_x + 2; ++i)
        {
            i_fstream << i_t << ' ' << i_y0 + (j - 2) * i_dy << ' ' << i_x0 + (i - 2) * i_dx << ' ' << i_inputVec[j][i][i_idx] << std::endl;
        }

        i_fstream << std::endl;
    }
    i_fstream << std::endl;
    i_fstream << std::endl;
}

int main()
{
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
    tStop = 0.3;
    SLIC_2D_EulerSolver testSolverClass(nCells_x, nCells_y, x0, x1, y0, y1, c, tStop);

    std::vector<std::vector<std::array<double, 4>>> compDomain;
    compDomain.resize(nCells_y + 4);
    for (int i = 0; i < nCells_y + 4; ++i)
    {
        compDomain[i].resize(nCells_x + 4);
    }

    setInitRiemannValue(compDomain, x0, x1, y0, y1);
    conservativeFormTransform(compDomain);

    testSolverClass.uVec_ = compDomain;
    testSolverClass.updateBoundary();
    // DEBUG_printOutVecRho(testSolverClass.uVec_);

    // now the initial conditions are ready

    // setting fstream
    std::ofstream rhoResults;
    std::ofstream momentumX_Results;
    std::ofstream momentumY_Results;
    std::ofstream energyResults;
    std::string rhoResultsDir;
    std::string momentumX_ResultsDir;
    std::string momentumY_ResultsDir;
    std::string energyResultsDir;
    rhoResultsDir = "/Users/zianhuang/Room214N/dev/SLIC_EulerSolver2D/data/rhoResults.dat";
    momentumX_ResultsDir = "/Users/zianhuang/Room214N/dev/SLIC_EulerSolver2D/data/momentumX_Results.dat";
    momentumY_ResultsDir = "/Users/zianhuang/Room214N/dev/SLIC_EulerSolver2D/data/momentumY_Results.dat";
    energyResultsDir = "/Users/zianhuang/Room214N/dev/SLIC_EulerSolver2D/data/energyResults.dat";

    rhoResults.open(rhoResultsDir);
    momentumX_Results.open(momentumX_ResultsDir);
    momentumY_Results.open(momentumY_ResultsDir);
    energyResults.open(energyResultsDir);

    // writing initial conditions to files
    writeToFileStream(rhoResults, compDomain, x0, testSolverClass.dx_, y0, testSolverClass.dy_, 0, 0);
    writeToFileStream(momentumX_Results, compDomain, x0, testSolverClass.dx_, y0, testSolverClass.dy_, 0, 1);
    writeToFileStream(momentumY_Results, compDomain, x0, testSolverClass.dx_, y0, testSolverClass.dy_, 0, 2);
    writeToFileStream(energyResults, compDomain, x0, testSolverClass.dx_, y0, testSolverClass.dy_, 0, 3);

    // preparing the while loop
    double t = 0.0;
    double localDt;
    int numIter = 0;
    do
    {
        if (numIter < 10)
        {
            testSolverClass.c_ = 0.3;
            ++numIter;

            localDt = testSolverClass.calcTimeStep();
            t += localDt;
            std::cout << t << std::endl;

            testSolverClass.uVec_ = testSolverClass.vecTran.slicVecTran_x(testSolverClass.uVec_, testSolverClass.dx_, localDt);
            testSolverClass.updateBoundary();
            testSolverClass.uVec_ = testSolverClass.vecTran.slicVecTran_y(testSolverClass.uVec_, testSolverClass.dy_, localDt);
            testSolverClass.updateBoundary();

            writeToFileStream(rhoResults, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 0);
            writeToFileStream(momentumX_Results, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 1);
            writeToFileStream(momentumY_Results, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 2);
            writeToFileStream(energyResults, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 3);
        }
        else
        {
            testSolverClass.c_ = 0.8;

            localDt = testSolverClass.calcTimeStep();
            t += localDt;
            std::cout << t << std::endl;

            testSolverClass.uVec_ = testSolverClass.vecTran.slicVecTran_x(testSolverClass.uVec_, testSolverClass.dx_, localDt);
            testSolverClass.updateBoundary();
            testSolverClass.uVec_ = testSolverClass.vecTran.slicVecTran_y(testSolverClass.uVec_, testSolverClass.dy_, localDt);
            testSolverClass.updateBoundary();

            writeToFileStream(rhoResults, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 0);
            writeToFileStream(momentumX_Results, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 1);
            writeToFileStream(momentumY_Results, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 2);
            writeToFileStream(energyResults, testSolverClass.uVec_, x0, testSolverClass.dx_, y0, testSolverClass.dy_, t, 3);
        }

    } while (t < testSolverClass.tStop_);

    rhoResults.close();
    momentumX_Results.close();
    momentumY_Results.close();
    energyResults.close();

    return 0;
}