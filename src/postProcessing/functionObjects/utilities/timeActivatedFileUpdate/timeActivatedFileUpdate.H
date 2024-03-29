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
    Foam::timeActivatedFileUpdate

Group
    grpUtilitiesFunctionObjects

Description
    Performs a file copy/replacement once a specified time has been reached.

    Example usage to update the fvSolution dictionary at various times
    throughout the calculation:

    \verbatim
    fileUpdate1
    {
        type              timeActivatedFileUpdate;
        functionObjectLibs ("libutilityFunctionObjects.so");
        outputControl     timeStep;
        outputInterval    1;
        fileToUpdate      "$FOAM_CASE/system/fvSolution";
        timeVsFile
        (
            (-1 "$FOAM_CASE/system/fvSolution.0")
            (0.10 "$FOAM_CASE/system/fvSolution.10")
            (0.20 "$FOAM_CASE/system/fvSolution.20")
            (0.35 "$FOAM_CASE/system/fvSolution.35")
        );
    }
    \endverbatim

SourceFiles
    timeActivatedFileUpdate.C
    IOtimeActivatedFileUpdate.H

\*---------------------------------------------------------------------------*/

#ifndef timeActivatedFileUpdate_H
#define timeActivatedFileUpdate_H

#include "Tuple2.H"
#include "pointFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                  Class timeActivatedFileUpdate Declaration
\*---------------------------------------------------------------------------*/

class timeActivatedFileUpdate
{
    // Private data

        //- Name of this set of timeActivatedFileUpdate objects
        word name_;

        //- Owner database
        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Name of file to update
        fileName fileToUpdate_;

        //- List of times vs filenames
        List<Tuple2<scalar, fileName> > timeVsFile_;

        //- Index of last file copied
        label lastIndex_;


    // Private Member Functions

        //- Update file
        void updateFile();

        //- Disallow default bitwise copy construct
        timeActivatedFileUpdate(const timeActivatedFileUpdate&);

        //- Disallow default bitwise assignment
        void operator=(const timeActivatedFileUpdate&);


public:

    //- Runtime type information
    TypeName("timeActivatedFileUpdate");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        timeActivatedFileUpdate
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~timeActivatedFileUpdate();


    // Member Functions

        //- Return name of the set of timeActivatedFileUpdate
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the timeActivatedFileUpdate data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Calculate the timeActivatedFileUpdate and write
        virtual void write();

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
