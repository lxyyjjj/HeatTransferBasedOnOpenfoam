/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2025 OpenCFD Ltd.
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
    foamListRegions

Group
    grpPostProcessingUtilities

Description
    List regions from constant/regionProperties.

Usage
    \b foamListRegions [OPTION]

Note
    The OpenFOAM banner information is suppressed so that the output can be
    piped into another command.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "regionProperties.H"

// Same as faMesh::prefix() but without additional linkage
constexpr const char* const faMeshPrefix = "finite-area";

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "List volume regions from constant/regionProperties,\n"
        "or area regions from constant/finite-area/regionProperties"
    );

    argList::noBanner();
    argList::noParallel();
    argList::noJobInfo();
    argList::noFunctionObjects();  // Never use function objects
    // No profiling since there is no time loop

    argList::addBoolOption
    (
        "finite-area",
        "List constant/finite-area/regionProperties (if available)"
    );

    argList::addDryRunOption
    (
        "Make reading optional and add verbosity"
    );
    argList::addVerboseOption("Additional verbosity");

    // Arguments are optional (non-mandatory)
    argList::noMandatoryArgs();
    argList::addArgument("regionType ... regionType");

    #include "setRootCase.H"

    const bool dryRun = args.dryRun();
    int optVerbose = args.verbose();

    if (dryRun && !optVerbose)
    {
        ++optVerbose;
    }

    // Use finite-area regions
    const bool doFiniteArea = args.found("finite-area");

    IOobjectOption::readOption readOpt(IOobjectOption::MUST_READ);

    if (dryRun || doFiniteArea)
    {
        // Always treat finite-area regionProperties as optional
        readOpt = IOobjectOption::READ_IF_PRESENT;
    }

    // Silent version of "createTime.H", without libraries
    Time runTime
    (
        Time::controlDictName,
        args,
        false,  // no enableFunctionObjects
        false   // no enableLibs
    );

    regionProperties regionProps;
    if (doFiniteArea)
    {
        regionProps = regionProperties(runTime, faMeshPrefix, readOpt);

        if (regionProps.empty())
        {
            InfoErr<< "No finite-area region types" << nl;
        }
        else if (optVerbose)
        {
            InfoErr
                << "Have " << regionProps.size()
                << " finite-area region types, "
                << regionProps.count() << " regions" << nl << nl;
        }
    }
    else
    {
        regionProps = regionProperties(runTime, readOpt);

        if (regionProps.empty())
        {
            // Probably only occurs with -dry-run option
            InfoErr<< "No region types" << nl;
        }
        else if (optVerbose)
        {
            InfoErr
                << "Have " << regionProps.size() << " region types, "
                << regionProps.count() << " regions" << nl << nl;
        }
    }

    // We now handle checking args and general sanity etc.
    wordList regionTypes;

    if (args.size() > 1)
    {
        regionTypes.resize(args.size()-1);

        // No duplicates
        wordHashSet uniq;

        label nTypes = 0;
        for (label argi = 1; argi < args.size(); ++argi)
        {
            regionTypes[nTypes] = args[argi];

            const word& regType = regionTypes[nTypes];

            if (uniq.insert(regType))
            {
                if (regionProps.contains(regType))
                {
                    ++nTypes;
                }
                else
                {
                    InfoErr<< "No region-type: " << regType << nl;
                }
            }
        }

        regionTypes.resize(nTypes);
    }
    else
    {
        regionTypes = regionProps.sortedToc();
    }


    for (const word& regionType : regionTypes)
    {
        const wordList& regionNames = regionProps[regionType];

        for (const word& regionName : regionNames)
        {
            Info<< regionName << nl;
        }
    }

    return 0;
}


// ************************************************************************* //
