    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    UEqn =
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      - fvm::laplacian(thermo.mu(), U)
      - fvc::div(thermo.mu()*dev2(T(fvc::grad(U))))
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (momentumPredictor)
    {
        solve
        (
            UEqn
        ==
            fvc::reconstruct
            (
                (
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
