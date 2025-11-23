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

Class
    Foam::GlobalOffset

Description
    A Foam::OffsetRange tuple of start/size/total with some additional
    parallel functionality.

    In general use OffsetRange for parameters requiring addressing
    information, and use GlobalOffset for creating and managing
    that addressing locally.

SourceFiles
    globalOffset.txx
    globalOffsetI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_globalOffset_H
#define Foam_globalOffset_H

#include "UPstream.H"
#include "OffsetRange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class IntType> class GlobalOffset;

template<class IntType>
Istream& operator>>(Istream& os, GlobalOffset<IntType>&);

template<class IntType>
Ostream& operator<<(Ostream& os, const GlobalOffset<IntType>&);

// Common types

//- A GlobalOffset with label for the addressing
typedef GlobalOffset<label> globalOffset;


/*---------------------------------------------------------------------------*\
                        Class GlobalOffset Declaration
\*---------------------------------------------------------------------------*/

template<class IntType>
class GlobalOffset
:
    public OffsetRange<IntType>
{
public:

    // Generated Methods: copy/move construct, copy/move assignment

    // Constructors

        //- Default construct as (0,0,0)
        inline constexpr GlobalOffset() noexcept;

        //- Construct with specified length, starting at zero (0,len,len).
        //- No communication.
        inline constexpr GlobalOffset(IntType len) noexcept;

        //- Read construct from start/size/total tuple
        explicit GlobalOffset(Istream& is);


    // Member Functions

    // Edit

        //- Reset to zero
        inline void clear() noexcept;

        //- Reset from local size with zero offset.
        inline void reset(IntType len = 0) noexcept;


    // Some functionality as per globalIndex

    // Local queries and renumbering

        //- From local to global index
        template<class IntType2>
        IntType2 toGlobal(IntType2 i) const noexcept
        {
            return (this->start() + i);
        }

        //- Convert a list from local to global index
        template<class IntType2>
        List<IntType2> toGlobal(const UList<IntType2>& labels) const;

        //- Convert a list (inplace) from local to global index
        template<class IntType2>
        void inplaceToGlobal(UList<IntType2>& labels) const;

        //- From global to local
        //  FatalError for an out-of-range value
        template<class IntType2>
        IntType2 toLocal(const IntType2 i) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalOffsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalOffset.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
