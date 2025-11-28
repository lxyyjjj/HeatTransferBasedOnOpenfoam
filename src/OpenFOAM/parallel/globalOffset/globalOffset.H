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

#include "OffsetRange.H"
#include "OffsetRangeIO.H"
#include "UPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class IntType> class GlobalOffset;

// Common Types
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

    // Static Functions

        //- Return a null GlobalOffset (reference to a nullObject).
        //- Behaves like an empty range.
        static const GlobalOffset<IntType>& null() noexcept
        {
            return NullObjectRef<GlobalOffset<IntType>>();
        }


    // Generated Methods: copy/move construct, copy/move assignment

    // Constructors

        //- Default construct as (0,0,0)
        inline constexpr GlobalOffset() noexcept;

        //- Construct with specified length, starting at zero (0,len,len).
        //- No communication.
        inline constexpr GlobalOffset(IntType len) noexcept;

        //- Construct from local size, with parallel coordination of
        //- the offsets and total size
        GlobalOffset
        (
            IntType len,
            int communicator,                   //!< communicator
            bool parallel = UPstream::parRun()  //!< use parallel comms
        );

        //- Copy construct from OffsetRange.
        //  Should be left as implicit to allow assignment from
        //  the results of the calculate() method
        inline GlobalOffset(const OffsetRange<IntType>& range) noexcept;

        //- Read construct from start/size/total tuple
        explicit GlobalOffset(Istream& is);


    // Member Functions

    // Access

        //- Same as total() - method name as per globalIndex
        IntType totalSize() const noexcept { return this->total(); }


    // Edit

        //- Reset from local size with zero offset.
        inline void reset(IntType len) noexcept;

        //- Reset from local size, with parallel coordination of
        //- the offsets and total size
        inline void reset
        (
            IntType len,
            int communicator,                   //!< communicator
            bool parallel = UPstream::parRun()  //!< use parallel comms
        );


    // Parallel Functions

        //- Based on the local input size,
        //- determine the globally-consistent start offset and total size
        void reduce(int communicator = UPstream::worldComm);


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


    // Helper Functions

        //- Calculate a globally-consistent OffsetRange based on the
        //- the local input size.
        static OffsetRange<IntType> calcOffsetRange
        (
            IntType localSize,
            const int communicator = UPstream::worldComm  //!< communicator
        );
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalOffsetI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalOffset.txx"
#endif

// * * * * * * * * * * * * * * * Global Reduction  * * * * * * * * * * * * * //

namespace Foam
{

//- Parallel reduction of OffsetRange (eg, GlobalOffset) on its size to
//- yield the globally-consistent offset and the total size across all ranks
template<class IntType>
void reduceOffset
(
    OffsetRange<IntType>& range,
    const int communicator = UPstream::worldComm
)
{
    // Single-item reduction
    Foam::PstreamDetail::reduce_offsetRange<IntType>
    (
        range,
        communicator
    );
}


//- Parallel reduction of multiple OffsetRange (eg, GlobalOffset) items.
//  To avoid MPI multiple times it packs/unpacks values into a local work
//  array.
//
//  \note OffsetRange is already restricted to integral types,
//  but this wrapper also requires that the parameters have the same
//  representation (ie, cannot mix int32/int64, signed/unsigned etc)
template
<
    class First,
    class... Rest,
    class = std::enable_if_t
    <
    (
        std::is_integral_v<typename First::value_type>
     && std::is_base_of_v
        <Foam::OffsetRange<typename First::value_type>, First>
     && (std::is_base_of_v
        <Foam::OffsetRange<typename First::value_type>, Rest> && ...)
    )>
>
void reduceOffsets
(
    const int communicator,
    First& first,
    Rest&... rest
)
{
    if constexpr (sizeof...(rest) == 0)
    {
        // Single-item reduction
        Foam::PstreamDetail::reduce_offsetRange<typename First::value_type>
        (
            first,
            communicator
        );
    }
    else
    {
        // Reduce multiple items, use pack/unpack
        Foam::PstreamDetail::reduce_offsetRanges<typename First::value_type>
        (
            communicator,
            std::make_index_sequence<1 + sizeof...(rest)>{},
            first,
            rest...
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
