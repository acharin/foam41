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

#include "MGridGenGAMGAgglomeration.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MGridGenGAMGAgglomeration, 0);

    addToRunTimeSelectionTable
    (
        GAMGAgglomeration,
        MGridGenGAMGAgglomeration,
        lduMesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MGridGenGAMGAgglomeration::MGridGenGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& dict
)
:
    GAMGAgglomeration(mesh, dict),
    fvMesh_(refCast<const fvMesh>(mesh))
{
    // Min, max size of agglomerated cells
    label minSize(readLabel(dict.lookup("minSize")));
    label maxSize(readLabel(dict.lookup("maxSize")));


    // Get the finest-level interfaces from the mesh
    interfaceLevels_.set
    (
        0,
        new lduInterfacePtrsList(fvMesh_.boundary().interfaces())
    );

    // Start geometric agglomeration from the cell volumes and areas of the mesh
    scalarField* VPtr = const_cast<scalarField*>(&fvMesh_.cellVolumes());
    vectorField* SfPtr = const_cast<vectorField*>(&fvMesh_.faceAreas());

    // Create the boundary area cell field
    scalarField* SbPtr(new scalarField(fvMesh_.nCells(), 0));

    {
        scalarField& Sb = *SbPtr;

        const labelList& own = fvMesh_.faceOwner();
        const vectorField& Sf = fvMesh_.faceAreas();

        forAll(Sf, facei)
        {
            if (!fvMesh_.isInternalFace(facei))
            {
                Sb[own[facei]] += mag(Sf[facei]);
            }
        }
    }


    // Agglomerate until the required number of cells in the coarsest level
    // is reached

    label nCreatedLevels = 0;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        label nCoarseCells = -1;

        tmp<labelField> finalAgglomPtr = agglomerate
        (
            nCoarseCells,
            minSize,
            maxSize,
            meshLevel(nCreatedLevels).lduAddr(),
            *VPtr,
            *SfPtr,
            *SbPtr
        );

        if (continueAgglomerating(nCoarseCells))
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        agglomerateLduAddressing(nCreatedLevels);

        // Agglomerate the cell volumes field for the next level
        {
            scalarField* aggVPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            restrictField(*aggVPtr, *VPtr, nCreatedLevels);

            if (nCreatedLevels)
            {
                delete VPtr;
            }

            VPtr = aggVPtr;
        }

        // Agglomerate the face areas field for the next level
        {
            vectorField* aggSfPtr
            (
                new vectorField
                (
                    meshLevels_[nCreatedLevels].upperAddr().size(),
                    vector::zero
                )
            );

            restrictFaceField(*aggSfPtr, *SfPtr, nCreatedLevels);

            if (nCreatedLevels)
            {
                delete SfPtr;
            }

            SfPtr = aggSfPtr;
        }

        // Agglomerate the cell boundary areas field for the next level
        {
            scalarField* aggSbPtr
            (
                new scalarField(meshLevels_[nCreatedLevels].size())
            );

            restrictField(*aggSbPtr, *SbPtr, nCreatedLevels);

            delete SbPtr;
            SbPtr = aggSbPtr;
        }

        nCreatedLevels++;
    }

    // Shrink the storage of the levels to those created
    compactLevels(nCreatedLevels);

    // Delete temporary geometry storage
    if (nCreatedLevels)
    {
        delete VPtr;
        delete SfPtr;
    }
    delete SbPtr;
}


// ************************************************************************* //
