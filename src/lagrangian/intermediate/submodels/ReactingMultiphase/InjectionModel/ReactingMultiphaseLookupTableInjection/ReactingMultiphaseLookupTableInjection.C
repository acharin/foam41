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

#include "ReactingMultiphaseLookupTableInjection.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label
Foam::ReactingMultiphaseLookupTableInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return round(injectorCells_.size()*(time1 - time0)*nParcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar
Foam::ReactingMultiphaseLookupTableInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    scalar volume = 0.0;
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        forAll(injectors_, i)
        {
            volume += injectors_[i].mDot()/injectors_[i].rho()*(time1 - time0);
        }
    }

    return volume;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseLookupTableInjection<CloudType>::
ReactingMultiphaseLookupTableInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    inputFileName_(this->coeffDict().lookup("inputFile")),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    nParcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    injectors_
    (
        IOobject
        (
            inputFileName_,
            owner.db().time().constant(),
            owner.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCells_(0)
{
    // Set/cache the injector cells
    injectorCells_.setSize(injectors_.size());
    forAll(injectors_, i)
    {
        this->findCellAtPosition(injectorCells_[i], injectors_[i].x());
    }

    // Determine volume of particles to inject
    this->volumeTotal_ = 0.0;
    forAll(injectors_, i)
    {
        this->volumeTotal_ += injectors_[i].mDot()/injectors_[i].rho();
    }
    this->volumeTotal_ *= duration_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingMultiphaseLookupTableInjection<CloudType>::
~ReactingMultiphaseLookupTableInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ReactingMultiphaseLookupTableInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar
Foam::ReactingMultiphaseLookupTableInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::ReactingMultiphaseLookupTableInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner
)
{
    label injectorI = parcelI*injectorCells_.size()/nParcels;

    position = injectors_[injectorI].x();
    cellOwner = injectorCells_[injectorI];
}


template<class CloudType>
void Foam::ReactingMultiphaseLookupTableInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    label injectorI = parcelI*injectorCells_.size()/nParcels;

    // set particle velocity
    parcel.U() = injectors_[injectorI].U();

    // set particle diameter
    parcel.d() = injectors_[injectorI].d();

    // set particle density
    parcel.rho() = injectors_[injectorI].rho();

    // set particle temperature
    parcel.T() = injectors_[injectorI].T();

    // set particle specific heat capacity
    parcel.cp() = injectors_[injectorI].cp();

    // set particle component mass fractions
    parcel.Y() = injectors_[injectorI].Y();

    // set particle gaseous component mass fractions
    parcel.YGas() = injectors_[injectorI].YGas();

    // set particle liquid component mass fractions
    parcel.YLiquid() = injectors_[injectorI].YLiquid();

    // set particle solid component mass fractions
    parcel.YSolid() = injectors_[injectorI].YSolid();
}


template<class CloudType>
bool
Foam::ReactingMultiphaseLookupTableInjection<CloudType>::fullyDescribed() const
{
    return true;
}


template<class CloudType>
bool Foam::ReactingMultiphaseLookupTableInjection<CloudType>::validInjection
(
    const label
)
{
    return true;
}


// ************************************************************************* //
