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

Class
    Foam::cellSets

Description

\*---------------------------------------------------------------------------*/

#ifndef cellSets_H
#define cellSets_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellSets Declaration
\*---------------------------------------------------------------------------*/

class cellSets
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSets(const cellSets&);

        //- Disallow default bitwise assignment
        void operator=(const cellSets&);


public:

        label nHexesWedges;
        label nPrisms;
        label nPyrs;
        label nTets;
        label nPolys;

        labelList tets;
        labelList pyrs;
        labelList prisms;
        labelList wedges;
        labelList hexes;
        labelList polys;


    // Constructors

        //- Construct given the number ov cells
        cellSets(const label nCells)
        :
            nHexesWedges(0),
            nPrisms(0),
            nPyrs(0),
            nTets(0),
            nPolys(0),

            tets(nCells),
            pyrs(nCells),
            prisms(nCells),
            wedges(nCells),
            hexes(nCells),
            polys(nCells)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
