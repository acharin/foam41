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

#include "nutURoughWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutURoughWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    // The flow velocity at the adjacent cell centre
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    const scalarField& yPlus = tyPlus();

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw();

    forAll(yPlus, faceI)
    {
        if (yPlus[faceI] > yPlusLam_)
        {
            const scalar Re = magUp[faceI]*y[faceI]/nuw[faceI] + ROOTVSMALL;
            nutw[faceI] = nuw[faceI]*(sqr(yPlus[faceI])/Re - 1);
        }
    }

    return tnutw;
}


tmp<scalarField> nutURoughWallFunctionFvPatchScalarField::calcYPlus
(
    const scalarField& magUp
) const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus();

    if (roughnessHeight_ > 0.0)
    {
        // Rough Walls
        const scalar c_1 = 1/(90 - 2.25) + roughnessConstant_;
        static const scalar c_2 = 2.25/(90 - 2.25);
        static const scalar c_3 = 2.0*atan(1.0)/log(90/2.25);
        static const scalar c_4 = c_3*log(2.25);

        //if (KsPlusBasedOnYPlus_)
        {
            // If KsPlus is based on YPlus the extra term added to the law
            // of the wall will depend on yPlus
            forAll(yPlus, faceI)
            {
                const scalar magUpara = magUp[faceI];
                const scalar Re = magUpara*y[faceI]/nuw[faceI];
                const scalar kappaRe = kappa_*Re;

                scalar yp = yPlusLam_;
                const scalar ryPlusLam = 1.0/yp;

                label iter = 0;
                scalar yPlusLast = 0.0;
                scalar dKsPlusdYPlus = roughnessHeight_/y[faceI];

                // Additional tuning parameter - nominally = 1
                dKsPlusdYPlus *= roughnessFactor_;

                do
                {
                    yPlusLast = yp;

                    // The non-dimensional roughness height
                    scalar KsPlus = yp*dKsPlusdYPlus;

                    // The extra term in the law-of-the-wall
                    scalar G = 0.0;

                    scalar yPlusGPrime = 0.0;

                    if (KsPlus >= 90)
                    {
                        const scalar t_1 = 1 + roughnessConstant_*KsPlus;
                        G = log(t_1);
                        yPlusGPrime = roughnessConstant_*KsPlus/t_1;
                    }
                    else if (KsPlus > 2.25)
                    {
                        const scalar t_1 = c_1*KsPlus - c_2;
                        const scalar t_2 = c_3*log(KsPlus) - c_4;
                        const scalar sint_2 = sin(t_2);
                        const scalar logt_1 = log(t_1);
                        G = logt_1*sint_2;
                        yPlusGPrime =
                            (c_1*sint_2*KsPlus/t_1) + (c_3*logt_1*cos(t_2));
                    }

                    scalar denom = 1.0 + log(E_*yp) - G - yPlusGPrime;
                    if (mag(denom) > VSMALL)
                    {
                        yp = (kappaRe + yp*(1 - yPlusGPrime))/denom;
                    }
                } while
                (
                    mag(ryPlusLam*(yp - yPlusLast)) > 0.0001
                 && ++iter < 10
                 && yp > VSMALL
                );

                yPlus[faceI] = max(0.0, yp);
            }
        }
    }
    else
    {
        // Smooth Walls
        forAll(yPlus, faceI)
        {
            const scalar magUpara = magUp[faceI];
            const scalar Re = magUpara*y[faceI]/nuw[faceI];
            const scalar kappaRe = kappa_*Re;

            scalar yp = yPlusLam_;
            const scalar ryPlusLam = 1.0/yp;

            label iter = 0;
            scalar yPlusLast = 0.0;

            do
            {
                yPlusLast = yp;
                yp = (kappaRe + yp)/(1.0 + log(E_*yp));

            } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 10);

            yPlus[faceI] = max(0.0, yp);
        }
    }

    return tyPlus;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    roughnessHeight_(0.0),
    roughnessConstant_(0.0),
    roughnessFactor_(0.0)
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const nutURoughWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    roughnessHeight_(ptf.roughnessHeight_),
    roughnessConstant_(ptf.roughnessConstant_),
    roughnessFactor_(ptf.roughnessFactor_)
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    roughnessHeight_(readScalar(dict.lookup("roughnessHeight"))),
    roughnessConstant_(readScalar(dict.lookup("roughnessConstant"))),
    roughnessFactor_(readScalar(dict.lookup("roughnessFactor")))
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const nutURoughWallFunctionFvPatchScalarField& rwfpsf
)
:
    nutWallFunctionFvPatchScalarField(rwfpsf),
    roughnessHeight_(rwfpsf.roughnessHeight_),
    roughnessConstant_(rwfpsf.roughnessConstant_),
    roughnessFactor_(rwfpsf.roughnessFactor_)
{}


nutURoughWallFunctionFvPatchScalarField::nutURoughWallFunctionFvPatchScalarField
(
    const nutURoughWallFunctionFvPatchScalarField& rwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(rwfpsf, iF),
    roughnessHeight_(rwfpsf.roughnessHeight_),
    roughnessConstant_(rwfpsf.roughnessConstant_),
    roughnessFactor_(rwfpsf.roughnessFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutURoughWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    tmp<scalarField> magUp = mag(Uw.patchInternalField() - Uw);

    return calcYPlus(magUp());
}


void nutURoughWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeKeyword("roughnessHeight")
        << roughnessHeight_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessConstant")
        << roughnessConstant_ << token::END_STATEMENT << nl;
    os.writeKeyword("roughnessFactor")
        << roughnessFactor_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutURoughWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
