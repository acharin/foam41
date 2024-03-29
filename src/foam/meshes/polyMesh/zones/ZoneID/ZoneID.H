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
    Foam::ZoneID

Description
    A class that holds the data needed to identify a zone in a dynamic mesh.

    The zone is identified by name.
    Its index in the zoneMesh is updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef ZoneID_H
#define ZoneID_H

#include "label.H"
#include "word.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ZoneType, class MeshType> class ZoneMesh;

// Forward declaration of friend functions and operators

template<class ZoneType> class ZoneID;

template<class ZoneType>
Ostream& operator<<(Ostream&, const ZoneID<ZoneType>&);

/*---------------------------------------------------------------------------*\
                           Class ZoneID Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType>
class ZoneID
{
    // Private data

        //- Zone name
        word name_;

        //- Zone index
        label index_;


public:

    // Constructors

        //- Construct from name
        ZoneID(const word& name, const ZoneMesh<ZoneType, polyMesh>& zm)
        :
            name_(name),
            index_(zm.findZoneID(name))
        {}

        //- Construct from Istream
        ZoneID(Istream& is, const ZoneMesh<ZoneType, polyMesh>& zm)
        :
            name_(is),
            index_(zm.findZoneID(name_))
        {}


    // Destructor - default


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return index
            label index() const
            {
                return index_;
            }

            //- Has the zone been found
            bool active() const
            {
                return index_ > -1;
            }

        // Edit

            //- Update
            void update(const ZoneMesh<ZoneType, polyMesh>& zm)
            {
                index_ = zm.findZoneID(name_);
            }


    // IOstream Operators

        friend Ostream& operator<< <ZoneType>
        (
            Ostream& os, const ZoneID<ZoneType>& p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ZoneType>
Ostream& operator<<
(
    Ostream& os, const ZoneID<ZoneType>& p
)
{
    os  << token::BEGIN_LIST
        << p.name_ << token::SPACE
        << p.index_
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const ZoneID<ZoneType>&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
