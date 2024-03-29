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

#include "transformList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class T>
Foam::List<T> Foam::transform
(
    const tensor& rotTensor,
    const UList<T>& field
)
{
    List<T> newField(field.size());

    forAll(field, i)
    {
        newField[i] = transform(rotTensor, field[i]);
    }

    return newField;
}


template <class T>
void Foam::transformList
(
    const tensorField& rotTensor,
    UList<T>& field
)
{
    if (rotTensor.size() == 1)
    {
        forAll(field, i)
        {
            field[i] = transform(rotTensor[0], field[i]);
        }
    }
    else if (rotTensor.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] = transform(rotTensor[i], field[i]);
        }
    }
    else
    {
        FatalErrorIn
        (
            "transformList(const tensorField&, UList<T>&)"
        )   << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


template <class T>
void Foam::transformList
(
    const tensorField& rotTensor,
    Map<T>& field
)
{
    if (rotTensor.size() == 1)
    {
        forAllIter(typename Map<T>, field, iter)
        {
            iter() = transform(rotTensor[0], iter());
        }
    }
    else
    {
        FatalErrorIn
        (
            "transformList(const tensorField&, Map<T>&)"
        )   << "Multiple transformation tensors not supported. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


template <class T>
void Foam::transformList
(
    const tensorField& rotTensor,
    EdgeMap<T>& field
)
{
    if (rotTensor.size() == 1)
    {
        forAllIter(typename EdgeMap<T>, field, iter)
        {
            iter() = transform(rotTensor[0], iter());
        }
    }
    else
    {
        FatalErrorIn
        (
            "transformList(const tensorField&, EdgeMap<T>&)"
        )   << "Multiple transformation tensors not supported. field:"
            << field.size() << " transformation:" << rotTensor.size()
            << abort(FatalError);
    }
}


// ************************************************************************* //
