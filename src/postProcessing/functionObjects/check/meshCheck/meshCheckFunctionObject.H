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
    meshCheckFunctionObject

Description
    FunctionObject checking the mesh

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    meshCheckFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef meshCheckFunctionObject_H
#define meshCheckFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class meshCheckFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class meshCheckFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Check topology
        Switch checkTopology_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshCheckFunctionObject(const meshCheckFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const meshCheckFunctionObject&);


public:

    //- Runtime type information
    TypeName("meshCheck");


    // Constructors

        //- Construct from components
        meshCheckFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute(const bool forceWrite);

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet();

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
