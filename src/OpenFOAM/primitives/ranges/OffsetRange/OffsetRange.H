/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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
    Foam::OffsetRange

Description
    A tuple of integers comprising start, size, total.
    This frequently corresponds to processor addressing
    (eg, Foam::globalIndex)
    but representing the local rank addressing only.

Note
    Only a minimum of IO operators are defined, to avoid incurring
    too many dependencies or cyclic dependencies.

SourceFiles
    OffsetRangeI.H
    OffsetRanges.cxx

\*---------------------------------------------------------------------------*/

#ifndef Foam_OffsetRange_H
#define Foam_OffsetRange_H

#include "labelFwd.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                         Class OffsetRange Declaration
\*---------------------------------------------------------------------------*/

template<class IntType>
class OffsetRange
{
    static_assert(std::is_integral_v<IntType>, "Integral required");

    // Private Data

        //- The start (offset)
        IntType start_;

        //- The size
        IntType size_;

        //- The total size
        IntType total_;


public:

    // STL type definitions

        //- Type of values contained
        typedef IntType value_type;

        //- The type that can represent the size
        typedef IntType size_type;


    // Generated Methods: copy/move construct, copy/move assignment

    // Constructors

        //- Default construct as (0,0,0)
        inline constexpr OffsetRange() noexcept;

        //- Construct with specified length, starting at zero.
        //- Sets total = size
        inline constexpr OffsetRange(IntType len) noexcept;

        //- Construct from start/length, no checks
        //- Sets total = size
        inline constexpr OffsetRange
        (
            IntType beg,
            IntType len
        ) noexcept;

        //- Construct from all components, no checks
        inline constexpr OffsetRange
        (
            IntType beg,
            IntType len,
            IntType tot
        ) noexcept;


    // Member Functions

    // Access

        //- True if zero-sized
        bool empty() const noexcept { return !size_; }

        //- The size
        IntType size() const noexcept { return size_; }

        //- Non-const access to the size
        IntType& size() noexcept { return size_; }

        //- The lower value of the range
        IntType start() const noexcept { return start_; }

        //- Non-const access to start of the range
        IntType& start() noexcept { return start_; }

        //- The total size
        IntType total() const noexcept { return total_; }

        //- Non-const access to the total size
        IntType& total() noexcept { return total_; }

        //- The (inclusive) lower value of the range,
        //- same as start(), begin_value()
        IntType min() const noexcept { return start_; }

        //- The (inclusive) upper value of the range,
        //- same as rbegin_value(). Ill-defined for an empty range
        IntType max() const noexcept { return (start_ + (size_ - 1)); }


    // Edit

        //- Reset to zero
        inline void clear() noexcept;

        //- Reset from size, start=0, total = size
        inline void reset(IntType len = 0) noexcept;

        //- Reset from start/size. Assigns total = size
        inline void reset(IntType beg, IntType len) noexcept;

        //- Reset all components
        inline void reset(IntType beg, IntType len, IntType tot) noexcept;


    // Comparison

        //- Test equality of start/size/total.
        inline bool equals(const OffsetRange& b) const noexcept;

        //- Compare start/size/total in that order.
        //  \return -1,0,+1
        inline int compare(const OffsetRange& b) const noexcept;


    // Search

        //- True if the (global) value is within the start/size interval.
        inline bool contains(IntType value) const noexcept;


    // Iterator ranges (as per IntRange)

        //- The value at the beginning of the range - same as start()
        inline IntType begin_value() const noexcept;

        //- The value 1 beyond the end of the range.
        inline IntType end_value() const noexcept;

        //- The max value of the range.
        inline IntType rbegin_value() const noexcept;

        //- The value 1 before the begin of range
        inline IntType rend_value() const noexcept;


    // Method name compatibility with globalIndex
    // (may change in the future)

        //- Same as start() - method name as per globalIndex
        IntType localStart() const noexcept { return start_; }

        //- Same as size() - method name as per globalIndex
        IntType localSize() const noexcept { return size_; }

        //- Same as total() - method name as per globalIndex
        IntType totalSize() const noexcept { return total_; }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Input operators for common integer types

//- Read IntRange from Istream as bracketed (start size total) tuple
Istream& operator>>(Istream& os, OffsetRange<int32_t>& range);

//- Read OffsetRange from Istream as bracketed (start size total) tuple
Istream& operator>>(Istream& os, OffsetRange<int64_t>& range);


// Output operators for common integer types

//- Write OffsetRange to Ostream as bracketed (start size total) tuple
Ostream& operator<<(Ostream& os, const OffsetRange<int32_t>& range);

//- Write OffsetRange to Ostream as bracketed (start size total) tuple
Ostream& operator<<(Ostream& os, const OffsetRange<int64_t>& range);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Test for equality of begin/size/total values
template<class IntType>
inline bool operator==
(
    const OffsetRange<IntType>& a,
    const OffsetRange<IntType>& b
) noexcept
{
    return a.equals(b);
}


//- Comparison function for sorting, compares the start.
//  If the start values are equal, also compares the size, etc
template<class IntType>
inline bool operator<
(
    const OffsetRange<IntType>& a,
    const OffsetRange<IntType>& b
) noexcept
{
    return (a.compare(b) < 0);
}


// Derived comparisons

template<class IntType>
inline bool operator!=
(
    const OffsetRange<IntType>& a,
    const OffsetRange<IntType>& b
) noexcept
{
    return !a.equals(b);
}

template<class IntType>
inline bool operator<=
(
    const OffsetRange<IntType>& a,
    const OffsetRange<IntType>& b
) noexcept
{
    return (a.compare(b) <= 0);
}

template<class IntType>
inline bool operator>
(
    const OffsetRange<IntType>& a,
    const OffsetRange<IntType>& b
) noexcept
{
    return (a.compare(b) > 0);
}

template<class IntType>
inline bool operator>=
(
    const OffsetRange<IntType>& a,
    const OffsetRange<IntType>& b
) noexcept
{
    return (a.compare(b) >= 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OffsetRangeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
