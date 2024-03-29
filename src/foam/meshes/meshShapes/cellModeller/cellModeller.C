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

Description
    Constructor of cellModeller: just sets the cellModeller's params.

\*---------------------------------------------------------------------------*/

#include "cellModeller.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellModeller::cellModeller()
{
    if (modelPtrs_.size())
    {
        FatalErrorIn("cellModeller::cellModeller(const fileName&)")
            << "attempt to re-construct cellModeller when it already exists"
            << exit(FatalError);
    }

    label maxIndex = 0;
    forAll(models_, i)
    {
        if (models_[i].index() > maxIndex) maxIndex = models_[i].index();
    }

    modelPtrs_.setSize(maxIndex + 1);
    modelPtrs_ = NULL;

    // For all the words in the wordlist, set the details of the model
    // to those specified by the word name and the other parameters
    // given. This should result in an automatic 'read' of the model
    // from its File (see cellModel class).
    forAll(models_, i)
    {
        if (modelPtrs_[models_[i].index()])
        {
            FatalErrorIn("cellModeller::cellModeller(const fileName&)")
                << "more than one model share the index "
                << models_[i].index()
                << exit(FatalError);
        }

        modelPtrs_[models_[i].index()] = &models_[i];

        if (modelDictionary_.found(models_[i].name()))
        {
            FatalErrorIn("cellModeller::cellModeller(const fileName&)")
                << "more than one model share the name "
                << models_[i].name()
                << exit(FatalError);
        }

        modelDictionary_.insert(models_[i].name(), &models_[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellModeller::~cellModeller()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Returns a pointer to a model which matches the string symbol
// supplied. A null pointer is returned if there is no suitable match.

const Foam::cellModel* Foam::cellModeller::lookup(const word& name)
{
    HashTable<const cellModel*>::iterator iter = modelDictionary_.find(name);

    if (iter != modelDictionary_.end())
    {
        return iter();
    }
    else
    {
        return NULL;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
