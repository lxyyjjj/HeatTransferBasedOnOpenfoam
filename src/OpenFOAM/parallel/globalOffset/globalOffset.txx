/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Mark Olesen
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "globalOffset.H"
#include <array>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IntType>
Foam::GlobalOffset<IntType>::GlobalOffset(Istream& is)
:
    OffsetRange<IntType>()
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IntType>
template<class IntType2>
Foam::List<IntType2> Foam::GlobalOffset<IntType>::toGlobal
(
    const UList<IntType2>& labels
) const
{
    // Or using std::transform

    //std::transform(labels.begin(), labels.end(), result.begin(),
    //    [=](auto id) { return id += start_ });

    List<IntType2> result(labels);
    inplaceToGlobal(result);
    return result;
}


template<class IntType>
template<class IntType2>
void Foam::GlobalOffset<IntType>::inplaceToGlobal
(
    UList<IntType2>& labels
) const
{
    if (const auto beg = this->start(); beg)
    {
        for (auto& val : labels)
        {
            val += beg;
        }
    }
}


template<class IntType>
template<class IntType2>
IntType2 Foam::GlobalOffset<IntType>::toLocal(const IntType2 i) const
{
    // !this->contains(i)
    if (i < this->begin_value() || this->end_value() <= i)
    {
        FatalErrorInFunction
            << "Global id:" << i << " not contained in the interval ["
            << this->begin_value() << "," << this->end_value() << "]\n"
            << abort(FatalError);
    }

    return (i - this->start());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class IntType>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::GlobalOffset<IntType>& rhs
)
{
    is.readBegin("GlobalOffset");
    is >> rhs.start() >> rhs.size() >> rhs.total();
    is.readEnd("GlobalOffset");

    is.check(FUNCTION_NAME);
    return is;
}


template<class IntType>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::GlobalOffset<IntType>& rhs
)
{
    os  << token::BEGIN_LIST
        << rhs.start() << token::SPACE
        << rhs.size() << token::SPACE
        << rhs.total()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
