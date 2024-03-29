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
    smallStrainOrthotropicSolidInterface

Description
    solidInterface procedure for small strain.
    no non-orthogonal correction.
    Interface tangential gradients are extrapolted from adjacent
    cell centres

SourceFiles
    smallStrainOrthotropicSolidInterface.C

Author
    Philip Cardiff UCD
    Original procedure by Tukovic et al 2012
\*---------------------------------------------------------------------------*/

#ifndef smallStrainOrthotropicSolidInterface_H
#define smallStrainOrthotropicSolidInterface_H

#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class smallStrainOrthotropicSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class smallStrainOrthotropicSolidInterface
:
    public solidInterface
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        smallStrainOrthotropicSolidInterface
        (const smallStrainOrthotropicSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const smallStrainOrthotropicSolidInterface&);

public:

    //- Runtime type information
    TypeName("smallStrainOrthotropic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        smallStrainOrthotropicSolidInterface
        (
            const word& name,
            const fvMesh& mesh,
            const constitutiveModel& rheology
        );


    // Destructor

        virtual ~smallStrainOrthotropicSolidInterface();


    // Member Functions

        // correct UEqn
        virtual void correct(fvVectorMatrix& UEqn);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
