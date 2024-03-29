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

#include "mixedIbFvPatchField.H"
#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void mixedIbFvPatchField<Type>::updateIbValues()
{
    // Interpolate the values from tri surface using nearest triangle
    const labelList& nt = ibPatch_.ibPolyPatch().nearestTri();

    this->refValue() = Field<Type>(triValue_, nt);
    this->refGrad() = Field<Type>(triGrad_, nt);
    this->valueFraction() = scalarField(triValueFraction_, nt);
}


template<class Type>
void mixedIbFvPatchField<Type>::setDeadValues()
{
    // Fix the value in dead cells
    if (setDeadValue_)
    {
        const labelList& dc = ibPatch_.ibPolyPatch().deadCells();

        // Get non-const access to internal field
        Field<Type>& psiI = const_cast<Field<Type>&>(this->internalField());

        forAll (dc, dcI)
        {
            psiI[dc[dcI]] = deadValue_;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixedIbFvPatchField<Type>::mixedIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    triValue_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    triGrad_(ibPatch_.ibMesh().size(), pTraits<Type>::zero),
    triValueFraction_(false),
    setDeadValue_(false),
    deadValue_(pTraits<Type>::zero)
{}


template<class Type>
mixedIbFvPatchField<Type>::mixedIbFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),   // Do not read mixed data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    triValue_("triValue", dict, ibPatch_.ibMesh().size()),
    triGrad_("triGradient", dict, ibPatch_.ibMesh().size()),
    triValueFraction_("triValueFraction", dict, ibPatch_.ibMesh().size()),
    setDeadValue_(dict.lookup("setDeadValue")),
    deadValue_(pTraits<Type>(dict.lookup("deadValue")))
{
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "mixedIbFvPatchField<Type>::"
            "mixedIbFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<Type>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();

    mixedFvPatchField<Type>::evaluate();
}


template<class Type>
mixedIbFvPatchField<Type>::mixedIbFvPatchField
(
    const mixedIbFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(p, iF),  // Do not map mixed data
    ibPatch_(refCast<const immersedBoundaryFvPatch>(p)),
    triValue_(ptf.triValue()),
    triGrad_(ptf.triGrad()),
    triValueFraction_(ptf.triValueFraction()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{
    // Note: NO MAPPING.  Fields are created on the immersed boundary
    // HJ, 12/Apr/2012
    if (!isType<immersedBoundaryFvPatch>(p))
    {
        FatalErrorIn
        (
            "mixedIbFvPatchField<Type>::"
            "mixedIbFvPatchField\n"
            "(\n"
            "    const mixedIbFvPatchField<Type>&,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }

    // Re-interpolate the data related to immersed boundary
    this->updateIbValues();

    // On creation of the mapped field, the internal field is dummy and
    // cannot be used.  Initialise the value to avoid errors
    // HJ, 1/Dec/2017
    Field<Type>::operator=(pTraits<Type>::zero);
}


template<class Type>
mixedIbFvPatchField<Type>::mixedIbFvPatchField
(
    const mixedIbFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    ibPatch_(ptf.ibPatch()),
    triValue_(ptf.triValue()),
    triGrad_(ptf.triGrad()),
    triValueFraction_(ptf.triValueFraction()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{}


template<class Type>
mixedIbFvPatchField<Type>::mixedIbFvPatchField
(
    const mixedIbFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    ibPatch_(ptf.ibPatch()),
    triValue_(ptf.triValue()),
    triGrad_(ptf.triGrad()),
    triValueFraction_(ptf.triValueFraction()),
    setDeadValue_(ptf.setDeadValue_),
    deadValue_(ptf.deadValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixedIbFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    // Base fields do not map: re-interpolate them from tri data
    this->updateIbValues();
}


template<class Type>
void mixedIbFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList&
)
{
    // Base fields do not rmap: re-interpolate them from tri data

    const mixedIbFvPatchField<Type>& mptf =
        refCast<const mixedIbFvPatchField<Type> >(ptf);

    // Set rmap tri data
    triValue_ = mptf.triValue_;
    triGrad_ = mptf.triGrad_;
    triValueFraction_ = mptf.triValueFraction_;

    this->updateIbValues();
}


template<class Type>
void mixedIbFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    this->updateIbValues();

    // Set dead value
    this->setDeadValues();

    // Evaluate mixed condition
    mixedFvPatchField<Type>::evaluate();
}


template<class Type>
void mixedIbFvPatchField<Type>::write(Ostream& os) const
{
    // Resolve post-processing issues.  HJ, 1/Dec/2017
    fvPatchField<Type>::write(os);
    os.writeKeyword("patchType")
        << immersedBoundaryFvPatch::typeName << token::END_STATEMENT << nl;
    triValue_.writeEntry("triValue", os);
    triGrad_.writeEntry("triGradient", os);
    triValueFraction_.writeEntry("triValueFraction", os);
    os.writeKeyword("setDeadValue")
        << setDeadValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("deadValue")
        << deadValue_ << token::END_STATEMENT << nl;

    // The value entry needs to be written with zero size
    Field<Type>::null().writeEntry("value", os);
    // this->writeEntry("value", os);

    // Write VTK on master only
    if (Pstream::master())
    {
        // Add parallel reduction of all faces and data to proc 0
        // and write the whola patch together

        // Write immersed boundary data as a vtk file
        autoPtr<surfaceWriter> writerPtr = surfaceWriter::New("vtk");

        // Get the intersected patch
        const standAlonePatch& ts = ibPatch_.ibPolyPatch().ibPatch();

        writerPtr->write
        (
            this->dimensionedInternalField().path(),
            ibPatch_.name(),
            ts.points(),
            ts,
            this->dimensionedInternalField().name(),
            *this,
            false, // FACE_DATA
            false  // verbose
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
