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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{
namespace PstreamDetail
{

// Reduction of OffsetRange is similar to globalIndex::calcOffsetTotal
// but retaining all of the values as member data
template<class IntType>
void reduce_offsetRange
(
    Foam::OffsetRange<IntType>& range,
    const int communicator  // The parallel communicator
)
{
    if (UPstream::is_parallel(communicator))
    {
        // Exscan (sum) yields the offsets, assigns 0 for rank=0
        IntType work = range.size();
        UPstream::mpiExscan_sum(&work, 1, communicator);

        // For the truly paranoid:
        // if (UPstream::master(communicator)) work = 0;

        range.start() = work;

        // The rank=(nProcs-1) knows the total - broadcast to others
        const auto root = (UPstream::nProcs(communicator)-1);
        if (root == UPstream::myProcNo(communicator))
        {
            // Update work as total == (start + size)
            work += range.size();
        }
        UPstream::broadcast(&work, 1, communicator, root);

        range.total() = work;
    }
}


// Implementation for reduceOffsets
//
// Equivalent to GlobalOffset::reduce (like globalIndex::calcOffsetTotal)
// but bundles values and performs operations on multiple values,
// which avoids calling MPI repeatedly
template
<
    class IntType,      // Should match GlobalOffset::value_type
    std::size_t... Is,
    class... OffsetRanges
>
void reduce_offsetRanges
(
    const int communicator,       // The parallel communicator
    std::index_sequence<Is...>,   // Indices into items
    OffsetRanges&... items
)
{
    if (UPstream::is_parallel(communicator))
    {
        using value_type = IntType;

        // Like globalIndex::calcOffsetTotal
        // but handling multiple items at once to reduce communication

        // Pack all sizes into the work buffer
        std::array<value_type, sizeof...(items)> work{ (items.size())... };

        // Exscan (sum) yields the offsets, assigns 0 for rank=0
        UPstream::mpiExscan_sum(work.data(), work.size(), communicator);

        // For the truly paranoid:
        // if (UPstream::master(communicator)) work.fill(0);

        // The work buffer now contains the offsets, copy back to starts
        ((items.start() = work[Is]), ...);

        // The rank=(nProcs-1) knows the total - broadcast to others
        const auto root = (UPstream::nProcs(communicator)-1);
        if (root == UPstream::myProcNo(communicator))
        {
            // Update work buffer as total == (start + size)
            ((work[Is] += items.size()), ...);
        }
        UPstream::broadcast(work.data(), work.size(), communicator, root);

        // The work buffer now contains the totals, copy back to total_
        ((items.total() = work[Is]), ...);
    }
}

} // End namespace PstreamDetail
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class IntType>
Foam::OffsetRange<IntType>
Foam::GlobalOffset<IntType>::calcOffsetRange
(
    IntType localSize,
    const int communicator
)
{
    OffsetRange<IntType> result(localSize);

    // Single-item reduction
    Foam::PstreamDetail::reduce_offsetRange(result, communicator);

    return result;
}


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
void Foam::GlobalOffset<IntType>::reduce
(
    const int communicator
)
{
    // Single-item reduction
    Foam::PstreamDetail::reduce_offsetRange(*this, communicator);
}


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


// ************************************************************************* //
