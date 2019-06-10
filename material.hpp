#pragma once

class MaterialParameters 
{
public:
    MaterialParameters() {}
    MaterialParameters(double youngsModulus, double poissonsRatio, 
        double hardening, double criticalCompression, double criticalStretch, double density);
    double E, nu, xsi, thetaC, thetaS, lambda, mu, rho;
};