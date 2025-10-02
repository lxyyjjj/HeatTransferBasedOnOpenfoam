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
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Gets single or multiple area region names based on command-line options:
    (-all-area-regions | -area-regions | -area-region)

Priority
    1. -all-area-regions
    2. -area-regions = specify multiple area regions to select,
        or a single area region
    3. -area-region = specify a single area region

Note
    There is no semantical difference between "-area-regions name"
    and "-area-region name".
    If regionProperties is missing (or empty) and the "-area-regions"
    option is supplied with literal names, these will be used directly
    instead of being used to filter the regionProperties.

Required Classes
    - Foam::polyMesh
    - Foam::regionProperties
    - Foam::IOobjectOption

Required Variables
    - args [argList]
    - runTime [Time]

Provides Variables
    - areaRegionNames [wordList]

See Also
    addAllFaRegionOptions.H

\*---------------------------------------------------------------------------*/

wordList areaRegionNames;
{
    // Local alias
    auto& names = areaRegionNames;
    // Exit or fallback if nothing matches
    constexpr bool exitOnNoMatches = false;

    if (args.found("all-area-regions"))
    {
        names =
        (
            Foam::regionProperties
            (
                runTime,
                "finite-area",  // ie, faMesh::prefix()
                IOobjectOption::READ_IF_PRESENT
            ).names()
        );

        if (names.empty())
        {
            InfoErr
                << "Warning: No finite-area regionProperties, "
                   "assume default region" << nl << endl;
        }
        else
        {
            Info<< "Using all area-regions: "
                << flatOutput(names) << nl;
        }
    }
    else if
    (
        wordRes select;
        args.readListIfPresent<wordRe>("area-regions", select)
    )
    {
        select.uniq();  // Normally a no-op

        if (select.size() == 1 && select.front().isLiteral())
        {
            // Identical to -area-region NAME
            names.resize(1);
            names.front() = select.front();
        }
        else if (select.size())
        {
            names =
            (
                Foam::regionProperties
                (
                    runTime,
                    "finite-area",  // ie, faMesh::prefix()
                    IOobjectOption::READ_IF_PRESENT
                ).names()
            );

            if (!names.empty())
            {
                // Filter the region names
                const auto indices = select.matching(names);

                if (indices.empty())
                {
                    InfoErr
                        << "No match in regions: "
                        << flatOutput(names) << nl;

                    if constexpr (exitOnNoMatches)
                    {
                        InfoErr<< "... stopping" << nl << endl;
                        return 1;
                    }
                    else
                    {
                        InfoErr<< "... ignoring and using default region"
                            << nl << endl;

                        // Arbitrary fallback, but can handle this further on
                        names.resize(1);
                        names.front() = Foam::polyMesh::defaultRegion;
                    }
                }
                else
                {
                    names = wordList(names, indices);
                    Info<< "Using area-regions: "
                        << flatOutput(names) << nl;
                }
            }
            else
            {
                // No information from regionProperties:
                // - use literal names from select

                names.resize(select.size());

                label nGood(0);
                for (const auto& val : select)
                {
                    if (val.isLiteral())
                    {
                        names[nGood++] = val;
                    }
                }
                if (nGood)
                {
                    names.resize(nGood);
                    InfoErr
                        << "Warning: No finite-area regionProperties, "
                           "using specified regions" << nl;

                    Info<< "Using area-regions: "
                        << flatOutput(names) << nl;
                }
                else
                {
                    names.resize(1);
                    names.front() = Foam::polyMesh::defaultRegion;

                    InfoErr
                        << "Warning: No finite-area regionProperties, "
                           "assume default region" << nl << endl;
                }
            }
        }
    }
    else if (args.found("area-region"))
    {
        names.resize(1);
        names.front() = args.get<word>("area-region");
    }

    // Fallback to defaultRegion
    if (names.empty())
    {
        names.resize(1);
        names.front() = Foam::polyMesh::defaultRegion;
    }
}


// ************************************************************************* //
