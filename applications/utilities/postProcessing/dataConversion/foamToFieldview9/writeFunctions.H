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

InClass
    Foam::writeFunctions

Description

SourceFiles
    writeFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef writeFunctions_H
#define writeFunctions_H

#include "scalar.H"
#include "foamString.H"
#include "labelList.H"
#include <iostream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool writeStr80(std::ostream& os, const string& str);
bool writeInt(std::ostream& os, int val1);
bool writeFloat(std::ostream& os, scalar val1);
void writeBytes(char* start, int nBytes);
void writeWallFlags(Ostream& os, label cellI, const labelList& wallFlags);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
