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

Application
    Test-parallel-file-write1

Description
    Simple test of writing with MPI/IO

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Switch.H"
#include "UPstreamFile.H"
#include "SpanStream.H"

using namespace Foam;

template<class IntType>
void zeropadded(std::ostream& os, IntType val, char fillch = '0')
{
    // set fill char and width
    os.setf(std::ios_base::right, std::ios_base::adjustfield);
    fillch = os.fill(fillch);
    os.width(std::numeric_limits<IntType>::digits10+1);
    os  << val;
    // restore fill char
    os.fill(fillch);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption();
    argList::addBoolOption("master-footer", "Write footer from master");

    #include "setRootCase.H"

    if (!UPstream::parRun())
    {
        Info<< "###############" << nl
            << "Not running in parallel. Stopping now" << nl
            << "###############" << endl;
        return 1;
    }

    const bool optMasterFooter = args.found("master-footer");

    Info<< nl << "Write master-footer: " << Switch::name(optMasterFooter)
        << nl << nl;

    Info<< "Create time (without controlDict)\n" << endl;

    auto runTimePtr = Time::New();
    auto& runTime = runTimePtr();

    const auto myProc = UPstream::myProcNo();
    const auto nProcs = UPstream::nProcs();

    // Some content
    OCharStream charset;
    for (int i = 0; i < 10; ++i)
    {
        charset<< char('A' + i);
    }

    // Header/footer buffers - these can be separate or bundled into
    // the first/last blocks

    OCharStream header;
    OCharStream footer;

    // Content buffer
    OCharStream os(IOstream::BINARY);

    {
        const auto v = charset.view();

        os << nl;
        os.beginBlock(word("rank" + Foam::name(myProc)));

        for (int repeat = 0;  repeat <= myProc; ++repeat)
        {
            os << indent << word("entry" + Foam::name(repeat))
                << ' ' << word("List<char>");
            // os << nl;
            os << ' ';
            os << label(v.size());
            os.write(v.data(), v.size());
            // os << nl;
            os.endEntry();
        }

        os.endBlock();
    }

    // Bundle the footer into the last block
    if (!optMasterFooter && (myProc == nProcs-1))
    {
        IOobject::writeEndDivider(os);
    }


    // All content now exists - commit to disk
    const std::string_view blockData(os.view());
    const int64_t blockSize(blockData.size());

    // Collect sizes
    const List<int64_t> sizes
    (
        UPstream::allGatherValues(blockSize, UPstream::worldComm)
    );


    // Format header with size information
    if (UPstream::master())
    {
        header
            << "Simple MPI/IO test with " << nProcs << " ranks" << nl << nl;

        ocharstream labelbuf;
        labelbuf.reserve_exact(32);

        // Position before writing a label
        auto labelBegin = header.tellp();

        header.beginBlock("meta");
        {
            header << indent << word("data.start") << ' ';

            labelBegin = header.tellp();

            // Add the start value (placeholder)
            {
                labelbuf.rewind();
                zeropadded(labelbuf, label(0));
                header.append(labelbuf.view());
            }

            header.endEntry();

            header << indent << word("data.sizes") << nl;
            sizes.writeList(header);  // flatOutput
            header.endEntry();
        }
        header.endBlock();

        header << nl;
        IOobject::writeDivider(header);

        // Now update with the correct size
        {
            labelbuf.rewind();
            zeropadded(labelbuf, label(header.view().size()));
            header.overwrite(labelBegin, labelbuf.view());
        }

        // Bundled the footer into the last block or from master?
        if (optMasterFooter)
        {
            IOobject::writeEndDivider(footer);
        }
    }

    // With additional header/footer
    int64_t headerSize(header.view().size());
    int64_t footerSize(footer.view().size());

    Pstream::broadcast(headerSize);
    if (optMasterFooter)
    {
        Pstream::broadcast(footerSize);
    }


    int64_t totalSize(headerSize);
    for (int i = 0; i < myProc; ++i)
    {
        totalSize += sizes[i];
    }

    const int64_t blockOffset(totalSize);

    for (int i = myProc; i < nProcs; ++i)
    {
        totalSize += sizes[i];
    }
    const int64_t footerOffset(totalSize);
    totalSize += footerSize;


    Pout<< "write size=" << label(blockSize)
        << " at=" << label(blockOffset) << " total=" << label(totalSize) << nl;

    {
        UPstream::File file;

        bool ok = file.open_write
        (
            UPstream::worldComm,
            runTime.globalPath()/"mpiio-test1.txt",
            IOstreamOption::ATOMIC
        );

        if (ok)
        {
            Info<< "writing: " << file.name() << nl;

            if (UPstream::master())
            {
                // A no-op for empty buffer
                ok = file.write_at(0, header.view());
            }

            ok = file.write_at_all(blockOffset, blockData);

            if (UPstream::master())
            {
                // A no-op for empty buffer
                ok = file.write_at(footerOffset, footer.view());
            }
        }

        file.set_size(totalSize);
        file.close();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
