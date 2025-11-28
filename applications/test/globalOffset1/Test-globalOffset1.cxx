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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-globalOffset1

Description
    Basic tests for OffsetRange and GlobalOffset
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelPair.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "IOstreams.H"

#include "globalIndex.H"
#include "globalOffset.H"

using namespace Foam;

template<class T>
void printInfo(const OffsetRange<T>& range)
{
    Pout<< "min:" << range.begin_value()
        << " max:" << range.rbegin_value()
        << " size:" << range.size()
        << " (total:" << range.total() << ')' << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noFunctionObjects();
    argList::noCheckProcessorDirectories();
    #include "setRootCase.H"
    #include "createTime.H"

    // Fails static_assert
    // Info<< "Default construct float: " << OffsetRange<float>().size() << nl;

    // Does not really make sense, but...
    // Info<< "Default construct bool: " << OffsetRange<bool>().size() << nl;

    typedef OffsetRange<int> intRange;
    typedef GlobalOffset<int> intOffset;

    Info<< "Default construct int32_t: " << OffsetRange<int32_t>() << nl
        << "Default construct int64_t: " << OffsetRange<int64_t>() << nl;

    Info<< "  one: " << intRange(10) << nl
        << "three: " << intRange(5, 10, 25) << nl;

    Pout<< nl;
    {
        intRange range(10);

        Pout<< "  input: " << range << nl;

        Foam::reduceOffset(range);

        Pout<< "  reduced: " << range << nl;

        printInfo(range);
    }

    // Now the same thing in one shot
    Pout<< nl;
    {
        intRange range0(4);
        intRange range1(10);
        intRange range2(12);
        intRange range3(8);

        Pout<< "  input: "
            << range0 << ", "
            << range1 << ", "
            << range2 << ", "
            << range3 << nl;

        Foam::reduceOffsets
        (
            UPstream::worldComm,
            range0, range1, range2, range3
        );

        Pout<< "  reduced: "
            << range0 << ", "
            << range1 << ", "
            << range2 << ", "
            << range3 << nl;
    }

    // Test mixing types (all derived from OffsetRange)
    Pout<< nl;
    {
        intOffset range0(4);
        intRange range1(10);

        [[maybe_unused]] OffsetRange<int32_t> range_32(4);
        [[maybe_unused]] OffsetRange<int64_t> range_64(4);
        [[maybe_unused]] label value(10);

        Pout<< "  input: "
            << range0 << ", "
            << range1 << nl;

        Foam::reduceOffsets
        (
            UPstream::worldComm,
            // Cannot mix representations - will not compile!
            // value,
            // range_32,
            // range_64,
            range0,
            range1
        );

        Pout<< "  reduced: "
            << range0 << ", "
            << range1 << nl;
    }

    // Other tests
    Pout<< nl;
    {
        label len1 = 10;
        label len2 = 15;

        globalOffset go1(len1);

        Pout<< "  value=" << len1 << " unreduced: " << go1 << endl;

        globalOffset go2(len2, UPstream::worldComm);

        Pout<< "  value=" << len2 << " reduced: " << go2 << endl;
        printInfo(go2);

        // Don't expect the generated move operators to do anything
        Pout<< "  move..." << endl;
        Pout<< "  old: " << go1 << endl;

        go1 = std::move(go2);
        Pout<< "  dst: " << go1 << endl;
        Pout<< "  src: " << go2 << endl;
    }

    Pout<< nl;
    {
        auto range4 = globalOffset::calcOffsetRange(25);
        Pout<< "  calcOffsetRange(25) = " << range4 << endl;
    }

    // Output tests for (mostly) arbitrary types
    {
        Info<< "range uint8: "
            << OffsetRange<direction>{sizeof(direction)} << nl
            << "range int16: " << OffsetRange<int16_t>{sizeof(int16_t)} << nl
            << "range uint32: " << OffsetRange<uint32_t>{sizeof(int32_t)} << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
