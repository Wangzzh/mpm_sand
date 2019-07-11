#include "material.hpp"

MaterialParameters::MaterialParameters(double youngsModulus, double poissonsRatio, 
        double hardening, double criticalCompression, double criticalStretch, double density, bool fluid) {
    this->E = youngsModulus;
    this->nu = poissonsRatio;
    this->xsi = hardening;
    this->thetaC = criticalCompression;
    this->thetaS = criticalStretch;
    this->lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    this->mu = E / 2 / (1 + nu);
    this->rho = density;
    this->fluid = fluid;
}