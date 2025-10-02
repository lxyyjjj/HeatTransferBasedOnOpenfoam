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
    Add multi-region command-line options:
    -all-area-regions, -area-regions, -area-region

Required Classes
    - Foam::argList

See Also
    getAllFaRegionOptions.H

\*---------------------------------------------------------------------------*/

{
    Foam::argList::addBoolOption
    (
        "all-area-regions",
        "Use all regions in finite-area regionProperties"
    );

    Foam::argList::addOption
    (
        "area-regions",
        "wordRes",
        "Use specified area region. Eg, -area-regions film\n"
        "Or from regionProperties.  Eg, -area-regions '(film \"solid.*\")'"
    );

    Foam::argList::addOption
    (
        "area-region",
        "name",
        "Specify area-mesh region. Eg, -region shell"
    );
}


// ************************************************************************* //
