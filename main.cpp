#include <QCoreApplication>

#include "godunovsolver.h"
#include "exactsolver.h"
#include "musclsolver.h"
#include "tests.h"
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    cout << "Please enter which test you would like to run:" << endl;
    cout << "1 - Toro's First Test" << endl;
    cout << "2 - Toro's Second Test" << endl;
    cout << "3 - Toro's Third Test" << endl;
    cout << "4 - Toro's Fourth Test" << endl;
    cout << "5 - Toro's Fifth Test" << endl;
    cout << "6 - Chinnayya et al.'s Test" << endl;
    cout << "7 - Fedkiw et al.'s Test" << endl;
    cout << "8 - Wang et al.'s Test" << endl;
    cout << "9 - Mach 10 Test" << endl;
    int test;
    cin >> test;

    cout << "Please select which solution type you would like:" << endl;

    if (test == 1 || test == 7)
    {
        cout << "1 - Exact" << endl;
        cout << "2 - Original Ghost Fluid Method" << endl;
        cout << "3 - Real Ghost Fluid Method" << endl;
    }

    if (test == 2 || test == 3 || test == 6)
    {
        cout << "1 - Exact" << endl;
        cout << "2 - Real Ghost Fluid Method" << endl;
    }

    if (test == 8 || test == 9)
    {
        cout << "1 - Original Ghost Fluid Method" << endl;
        cout << "2 - Real Ghost Fluid Method" << endl;
    }
    int solutionType;
    cin >> solutionType;

    cout << "Please enter the desired resolution:" << endl;
    int resolution;
    cin >> resolution;

    int reinitialisation;
    int reinitialisationIteration;
    if (((test == 1 || test == 3 || test == 6) && solutionType != 1) || test == 9)
    {
        cout << "Would you like to run this test with reinitialisation?" << endl;
        cout << "1 - With reinitialisation" << endl;
        cout << "2 - Without reinitialisation" << endl;
        cin >> reinitialisation;

        if (reinitialisation == 1)
        {
            cout << "After how many iterations would you like to reinitialise?" << endl;
            cin >> reinitialisationIteration;
        }
    }

    cout << "Computing!" << endl;

    if (test == 1)
    {
        if (solutionType == 1)
        {
            Tests::solveToroTest1Exact(resolution);
        }
        if (solutionType == 2)
        {
            if (reinitialisation == 1)
            {
                Tests::solveToroTest1BasicGhostFluidReinitialisation(resolution, reinitialisationIteration);
            }
            else
            {
                Tests::solveToroTest1BasicGhostFluid(resolution);
            }
        }
        if (solutionType == 3)
        {
            if (reinitialisation == 1)
            {
                Tests::solveToroTest1RiemannGhostFluidReinitialisation(resolution, reinitialisationIteration);
            }
            else
            {
                Tests::solveToroTest1RiemannGhostFluid(resolution);
            }
        }
    }

    if (test == 2)
    {
        if (solutionType == 1)
        {
            Tests::solveToroTest2Exact(resolution);
        }

        if (solutionType == 2)
        {
            Tests::solveToroTest2RiemannGhostFluid(resolution);
        }
    }

    if (test == 3)
    {
        if (solutionType == 1)
        {
            Tests::solveToroTest3Exact(resolution);
        }

        if (solutionType == 2)
        {
            if (reinitialisation == 1)
            {
                Tests::solveToroTest3RiemannGhostFluidReinitialisation(resolution, reinitialisationIteration);
            }
            else
            {
                Tests::solveToroTest3RiemannGhostFluid(resolution);
            }
        }
    }

    if (test == 4)
    {
        Tests::solveToroTest4Exact(resolution);
    }

    if (test == 5)
    {
        Tests::solveToroTest5Exact(resolution);
    }

    if (test == 6)
    {
        if (solutionType == 1)
        {
            Tests::solveChinnayyaTestExact(resolution);
        }

        if (solutionType == 2)
        {
            if (reinitialisation == 1)
            {
                Tests::solveChinnayyaTestRiemannGhostFluidReinitialisation(resolution, reinitialisationIteration);
            }
            else
            {
                Tests::solveChinnayyaTestRiemannGhostFluid(resolution);
            }
        }
    }

    if (test == 7)
    {
        if (solutionType == 1)
        {
            Tests::solveFedkiwTestExact(resolution);
        }

        if (solutionType == 2)
        {
            Tests::solveFedkiwTestBasicGhostFluid(resolution);
        }

        if (solutionType == 3)
        {
            Tests::solveFedkiwTestRiemannGhostFluid(resolution);
        }
    }

    if (test == 8)
    {
        if (solutionType == 1)
        {
            Tests::solveWangTestBasicGhostFluid(resolution);
        }

        if (solutionType == 2)
        {
            Tests::solveWangTestRiemannGhostFluid(resolution);
        }
    }

    if (test == 9)
    {
        if (solutionType == 1)
        {
            if (reinitialisation == 1)
            {
                Tests::solveMach10TestBasicGhostFluidReinitialisation(resolution, reinitialisationIteration);
            }
            else
            {
                Tests::solveMach10TestBasicGhostFluid(resolution);
            }
        }

        if (solutionType == 2)
        {
            if (reinitialisation == 1)
            {
                Tests::solveMach10TestRiemannGhostFluidReinitialisation(resolution, reinitialisationIteration);
            }
            else
            {
                Tests::solveMach10TestRiemannGhostFluid(resolution);
            }
        }
    }

    cout << "Computation completed!" << endl;

    return a.exec();
}
