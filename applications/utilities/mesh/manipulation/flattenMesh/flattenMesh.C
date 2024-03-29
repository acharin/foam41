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
    flattenMesh

Description
    Flattens the front and back planes of a 2D cartesian mesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "emptyPolyPatch.H"
#include "twoDPointCorrector.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    boundBox bb(points);

    Info<< "bounding box: min = " << bb.min()
        << " max = " << bb.max() << " meters."
        << endl;


    point midPoint = gAverage(points);

    twoDPointCorrector twoDCorr(mesh);

    direction planeNormalCmpt = twoDCorr.normalDir();

    scalar midCmptVal = midPoint[planeNormalCmpt];
    scalar minCmptVal = bb.min()[planeNormalCmpt];
    scalar maxCmptVal = bb.max()[planeNormalCmpt];

    forAll(points, pointI)
    {
        if (points[pointI][planeNormalCmpt] < midCmptVal)
        {
            points[pointI][planeNormalCmpt] = minCmptVal;
        }
        else
        {
            points[pointI][planeNormalCmpt] = maxCmptVal;
        }
    }

    twoDCorr.correctPoints(points);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing points into directory " << points.path() << nl << endl;
    points.write();

    return 0;
}


// ************************************************************************* //
