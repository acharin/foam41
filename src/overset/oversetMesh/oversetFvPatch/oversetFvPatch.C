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

#include "fvPatch.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "oversetFvPatch.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(oversetFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, oversetFvPatch, polyPatch);

}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::oversetFvPatch::makeWeights(fvsPatchScalarField& w) const
{
    fvPatch::makeWeights(w);
}


void Foam::oversetFvPatch::makeDeltaCoeffs(fvsPatchScalarField& dc) const
{
    fvPatch::makeDeltaCoeffs(dc);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::oversetFvPatch::delta() const
{
    return fvPatch::delta();
}


const Foam::oversetMesh& Foam::oversetFvPatch::overset() const
{
    return oversetMesh::New(boundaryMesh().mesh());
}

Foam::tmp<Foam::labelField> Foam::oversetFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


void Foam::oversetFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    // Missing
    return labelField::null();
}


void Foam::oversetFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{}


Foam::tmp<Foam::labelField> Foam::oversetFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return labelField::null();
}


// ************************************************************************* //
