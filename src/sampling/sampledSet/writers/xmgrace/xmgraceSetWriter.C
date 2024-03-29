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

\*---------------------------------------------------------------------------*/

#include "xmgraceSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::xmgraceSetWriter<Type>::xmgraceSetWriter()
:
    writer<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::xmgraceSetWriter<Type>::~xmgraceSetWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::xmgraceSetWriter<Type>::getFileName
(
    const coordSet& points,
    const wordList& valueSetNames
) const
{
    return this->getBaseName(points, valueSetNames) + ".agr";
}


template<class Type>
void Foam::xmgraceSetWriter<Type>::write
(
    const coordSet& points,
    const wordList& valueSetNames,
    const List<const Field<Type>*>& valueSets,
    Ostream& os
) const
{
    os  << "@g0 on" << nl
        << "@with g0" << nl
        << "@    title \"" << points.name() << '"' << nl
        << "@    xaxis label " << '"' << points.axis() << '"' << nl;

    forAll(valueSets, i)
    {
        os  << "@    s" << i << " legend " << '"'
            << valueSetNames[i] << '"' << nl
            << "@target G0.S" << i << nl;

        this->writeTable(points, *valueSets[i], os);

        os  << '&' << nl;
    }
}


template<class Type>
void Foam::xmgraceSetWriter<Type>::write
(
    const bool writeTracks,
    const PtrList<coordSet>& trackPoints,
    const wordList& valueSetNames,
    const List<List<Field<Type> > >& valueSets,
    Ostream& os
) const
{
    if (valueSets.size() != valueSetNames.size())
    {
        FatalErrorIn("gnuplotSetWriter<Type>::write(..)")
            << "Number of variables:" << valueSetNames.size() << endl
            << "Number of valueSets:" << valueSets.size()
            << exit(FatalError);
    }
    if (trackPoints.size() > 0)
    {
        os  << "@g0 on" << nl
            << "@with g0" << nl
            << "@    title \"" << trackPoints[0].name() << '"' << nl
            << "@    xaxis label " << '"' << trackPoints[0].axis() << '"' << nl;

        // Data index.
        label sI = 0;

        forAll(trackPoints, trackI)
        {
            forAll(valueSets, i)
            {
                os  << "@    s" << sI << " legend " << '"'
                    << valueSetNames[i] << "_track" << i << '"' << nl
                    << "@target G0.S" << sI << nl;
                this->writeTable(trackPoints[trackI], valueSets[i][trackI], os);
                os  << '&' << nl;

                sI++;
            }
        }
    }
}

// ************************************************************************* //
