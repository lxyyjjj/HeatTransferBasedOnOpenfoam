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

Application
    Test-parallel-scan

Description
    Simple tests for MPI_Scan/MPI_Exscan

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "argList.H"
#include "Time.H"
#include "vector.H"
#include "Pstream.H"
#include "globalIndex.H"

using namespace Foam;


template<class Type>
void testGlobalIndex(const Type& localSize)
{
    const auto comm = UPstream::worldComm;

    // With calcOffset
    {
        Type start = globalIndex::calcOffset(localSize, comm);

        Type total(start + localSize);

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, UPstream::nProcs(comm)-1);
        }

        Pout<< "size: " << localSize
            << " offset: " << start
            << " total: " << total
            << " (globalIndex calcOffset)" << endl;
    }

    // With calcRange
    {
        IntRange<Type> range = globalIndex::calcRange(localSize, comm);

        Type total = range.end_value();

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, UPstream::nProcs(comm)-1);
        }

        Pout<< "range: " << range
            << " total: " << total
            << " (globalIndex calcRange)" << endl;
    }
}


template<class Type>
void testScan(const Type& localValue, const bool exclusive)
{
    const auto comm = UPstream::worldComm;

    // sum
    {
        Type result = UPstream::mpiScan_sum(localValue, comm, exclusive);

        Type total(result);
        if (exclusive)
        {
            total += localValue;
        }

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, (UPstream::nProcs(comm)-1));
        }

        if (exclusive)
        {
            Pout<< "input: " << localValue
                << " beg-sum: " << result
                << " total: " << total << endl;
        }
        else
        {
            Pout<< "input: " << localValue
                << " end-sum: " << result
                << " total: " << total << endl;
        }
    }

    // min
    {
        Type result = UPstream::mpiScan_min(localValue, comm, exclusive);

        if (exclusive)
        {
            Pout<< "input: " << localValue
                << " (exclusive)min : " << result
                << endl;
        }
        else
        {
            Pout<< "input: " << localValue
                << " (inclusive)min : " << result
                << endl;
        }
    }

    // max (this type with templated version)
    {
        constexpr UPstream::opCodes opCode = UPstream::opCodes::op_max;

        Type result = UPstream::mpiScan<opCode>(localValue, comm, exclusive);

        if (exclusive)
        {
            Pout<< "input: " << localValue
                << " (exclusive)max : " << result
                << endl;
        }
        else
        {
            Pout<< "input: " << localValue
                << " (inclusive)max : " << result
                << endl;
        }
    }


    // With calcOffsetRange
    if constexpr (std::is_integral_v<Type>)
    {
        const auto& size = localValue;

        // Exscan ony!
        Type start = UPstream::mpiExscan_sum(size, comm);

        Type total(start + size);

        if (UPstream::is_parallel(comm))
        {
            // Broadcast total from rank=N-1
            UPstream::broadcast(&total, 1, comm, (UPstream::nProcs(comm)-1));
        }

        OffsetRange<Type> range(start, size, total);

        Pout<< "input: " << size
            << " range: " << range << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    argList::addBoolOption("exclusive", "Perform exclusive scan");
    argList::addOption("factor", "Multiplier for rank values");

    #include "setRootCase.H"

    if (!UPstream::parRun())
    {
        Info<< "parallel only!" << nl;
        return 1;
    }

    const auto myProci = UPstream::myProcNo();

    const bool exclusive = args.found("exclusive");
    const scalar scaleFactor = args.getOrDefault<scalar>("factor", 1);

    Info<< nl;
    {
        typedef label Type;

        Type localValue = (myProci+1)*scaleFactor;
        List<Type> allSizes = UPstream::listGatherValues(localValue);

        Info<< "all sizes: " << flatOutput(allSizes) << nl;

        testScan(localValue, exclusive);

        testGlobalIndex(localValue);
    }

    #if 0
    Info<< nl;
    {
        typedef int64_t Type;

        Type localValue = (myProci+1)*scaleFactor;
        List<Type> allSizes = UPstream::listGatherValues(localValue);

        Info<< "all sizes: " << flatOutput(allSizes) << nl;

        testScan(localValue, exclusive);

        testGlobalIndex(localValue);
    }
    #endif

    Info<< nl;
    {
        typedef vector Type;

        Type localValue = vector::uniform(myProci+1);
        localValue *= scaleFactor;

        testScan(localValue, exclusive);
    }

    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
