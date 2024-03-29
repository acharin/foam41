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
    deformedGeom

Description
    Deforms a polyMesh using a displacement field U and a scaling factor
    supplied as an argument.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "IStringStream.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"

    argList::validArgs.append("scaling factor");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedMesh.H"

    scalar scaleFactor(readScalar(IStringStream(args.additionalArgs()[0])()));

    instantList timeDirs = timeSelector::select0(runTime, args);

    pointField zeroPoints(mesh.points());

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        fvMesh::readUpdateState state = mesh.readUpdate();

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Create volPointInterpolation object after the mesh has been read
        volPointInterpolation pInterp(mesh);

        // Check U exists
        if (Uheader.headerOk())
        {
            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            pointField newPoints =
                zeroPoints
              + scaleFactor*pInterp.interpolate(U)().internalField();

            mesh.polyMesh::movePoints(newPoints);

            mesh.write();
        }
        else
        {
            Info<< "    No U" << endl;
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
