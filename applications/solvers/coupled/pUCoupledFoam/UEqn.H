{
    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    rAU = 1.0/UEqn.A();

    if (gMin(rAU) < 0)
    {
        WarningIn(args.executable())
            << "Negative rAU detected.  Check mesh and discretisation"
            << endl;
    }

    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    // Old term reinserted here // GGSF 21/08/18
    UEqn += fvm::SuSp(-divPhi,U) + divPhi*U;

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    UpEqn.insertEquation(0, UEqn);

#   include "addBlockCoupledBC.H"
}
