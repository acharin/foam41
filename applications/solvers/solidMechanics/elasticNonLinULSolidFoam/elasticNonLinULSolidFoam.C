/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    elasticNonLinULSolidFoam

Description
    Finite volume structural solver employing a incremental strain updated
    Lagrangian approach.

    Valid for small strains, finite displacements and finite rotations.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidInterface.H"
#include "volPointInterpolation.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "pointFields.H"
#include "twoDPointCorrector.H"
#include "leastSquaresVolPointInterpolation.H"
#include "processorFvPatchFields.H"
#include "transformGeometricField.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "readDivDSigmaExpMethod.H"
#   include "readDivDSigmaLargeStrainExpMethod.H"
#   include "readMoveMeshMethod.H"
#   include "createSolidInterfaceNonLin.H"
#   include "findGlobalFaceZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        int iCorr = 0;
        lduSolverPerformance solverPerf;
        scalar initialResidual = 1.0;
        scalar relativeResidual = 1.0;
        lduMatrix::debug = 0;

        do
        {
            DU.storePrevIter();

#           include "calculateDivDSigmaExp.H"
#           include "calculateDivDSigmaLargeStrainExp.H"

            //- Updated lagrangian momentum equation
            fvVectorMatrix DUEqn
            (
                fvm::d2dt2(rho,DU)
             ==
                fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
              + divDSigmaExp
              + divDSigmaLargeStrainExp
            );

            if (solidInterfaceCorr)
            {
                solidInterfacePtr->correct(DUEqn);
            }

            solverPerf = DUEqn.solve();

            if (iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            DU.relax();

            gradDU = fvc::grad(DU);

#           include "calculateDEpsilonDSigma.H"
#           include "calculateRelativeResidual.H"

            Info<< "\tTime " << runTime.value()
                << ", Corrector " << iCorr
                << ", Solving for " << DU.name()
                << " using " << solverPerf.solverName()
                << ", res = " << solverPerf.initialResidual()
                << ", rel res = " << relativeResidual
                << ", inner iters " << solverPerf.nIterations() << endl;
        }
        while
        (
            //solverPerf.initialResidual() > convergenceTolerance
            relativeResidual > convergenceTolerance
         && ++iCorr < nCorr
        );

      Info << nl << "Time " << runTime.value() << ", Solving for " << DU.name()
           << ", Initial residual = " << initialResidual
           << ", Final residual = " << solverPerf.initialResidual()
           << ", No outer iterations " << iCorr << endl;

#       include "moveMesh.H"
#       include "rotateFields.H"
#       include "writeFields.H"

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
