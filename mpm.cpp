#include "mpm.hpp"

MaterialParameters::MaterialParameters(double youngsModulus, double poissonsRatio, 
        double hardening, double criticalCompression, double criticalStretch) {
    this->E = youngsModulus;
    this->nu = poissonsRatio;
    this->xsi = hardening;
    this->thetaC = criticalCompression;
    this->thetaS = criticalStretch;
    this->lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    this->mu = E / 2 / (1 + nu);
}

MPM::MPM(int nGrid, double timeStep, MaterialParameters material) {
    this->nGrid = nGrid;
    this->time = 0;
    this->timeStep = timeStep;
    this->material = material;

    // Particle* p1 = new Particle();
    // p1 -> mass = 1;
    // p1 -> position << 0.47, 0.63;
    // p1 -> velocity << 1.0, 0.0;
    // particles.push_back(p1);
    
    // Particle* p2 = new Particle();
    // p2 -> mass = 2;
    // p2 -> position << 0.5, 0.52;
    // p2 -> velocity << 0.5, 0.5;
    // particles.push_back(p2);

    // Particle* p3 = new Particle();
    // p3 -> mass = 1;
    // p3 -> position << 0.81, 0.87;
    // p3 -> velocity << -0.2, 0.5;
    // particles.push_back(p3);

    int divs = 100;
    for (int i = 0; i < divs; i++) {
        Particle* p = new Particle();
        p -> mass = 0.04;
        double x = 0.5 + (rand() % 10000) / 100000.;
        double y = 0.8 + (rand() % 10000) / 100000.;
        p -> position << x, y;
        p -> velocity << 0, 0;
        particles.push_back(p);
    }

    grids = std::vector<std::vector<Grid*>>(nGrid, std::vector<Grid*>(nGrid));
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j] = new Grid();
            grids[i][j]->position << (((double)i + 0.5) / nGrid), (((double)j + 0.5) / nGrid);
            grids[i][j]->linearMomentum << i * 0.1, j * 0.1;
        }
    }
}

MPM::~MPM() {
    for (auto& particle : particles) {
        delete particle;
    }
    for (auto& gridVec : grids) {
        for (auto& grid : gridVec) {
            delete grid;
        }
    }
}

void MPM::step() {
    particleToGrid();
    // if (time == 0.) computeParticleDensity();
    computeGridForce();
    computeGridVelocity();
    updateDeformation();
    updateParticleVelocity();
    handleParticleCollision();
    updateParticlePosition();

    time += timeStep;
}

void MPM::render() {
    // Rendering particles
    glColor3f(1, 0, 0);
    for (auto& particle : particles) {
        particle->render();
    }
    
    // Rendering grids
    glColor3f(1, 1, 0);
    for (auto& gridVec : grids) {
        for (auto& grid : gridVec) {
            grid -> render();
        }
    }
}

void MPM::particleToGrid() {
    for (auto& gridVec : grids) {
        for (auto& grid : gridVec) {
            grid -> mass = 0;
            grid -> linearMomentum = Eigen::Vector2d::Zero();
        }
    }
    for (auto& particle : particles) {
        particle->xLeft = (int) (particle->position[0] * nGrid - 0.5);
        particle->yLeft = (int) (particle->position[1] * nGrid - 0.5);
        particle->xDiff = (particle->position[0] * nGrid) - particle->xLeft - 0.0;
        particle->yDiff = (particle->position[1] * nGrid) - particle->yLeft - 0.0;
        // particle->xDiff = (particle->position[0] * nGrid) - particle->xLeft - 0.5;
        // particle->yDiff = (particle->position[1] * nGrid) - particle->yLeft - 0.5;

        particle->calculateWeights();
        Eigen::Matrix2d Dinv = Eigen::Matrix2d::Identity() * 3 * nGrid * nGrid;
        
        for (int dx = -1; dx <= 2; dx++) {
            for (int dy = -1; dy <= 2; dy++) {
                int xGrid = particle->xLeft + dx;
                int yGrid = particle->yLeft + dy;
                Eigen::Vector2d displacement;
                displacement << (particle->xDiff - dx) / nGrid, (particle->yDiff - dy) / nGrid;
                if (xGrid >= 0 && xGrid < nGrid && yGrid >= 0 && yGrid < nGrid) {
                    grids[xGrid][yGrid] -> mass += particle->xWeight[dx+1] * particle->yWeight[dy+1] * particle->mass;
                    grids[xGrid][yGrid] -> linearMomentum += particle->xWeight[dx+1] * particle->yWeight[dy+1] * 
                        particle->mass * (particle-> velocity + particle->B * Dinv * displacement);
                }
                // std::cout << "displacement:" << displacement << std::endl;
            }
        }
    }
}

void MPM::computeParticleDensity() {
    // for (auto& particle : particles) {        
    //     particle->density = 0;
    //     for (int dx = -1; dx <= 2; dx++) {
    //         for (int dy = -1; dy <= 2; dy++) {
    //             int xGrid = particle->xLeft + dx;
    //             int yGrid = particle->yLeft + dy;
    //             if (xGrid >= 0 && xGrid < nGrid && yGrid >= 0 && yGrid < nGrid) {
    //                 particle->density += grids[xGrid][yGrid] -> mass * particle->xWeight[dx+1] * particle->yWeight[dy+1] * nGrid * nGrid;
    //             }
    //         }
    //     }
    //     particle->volume = particle->mass / particle->density;

    //     std::cout << "r " << particle->density << std::endl;
    //     std::cout << "v " << particle->volume << std::endl;
    // }
}

void MPM::computeGridForce() {
    Eigen::Vector2d gravity;
    gravity << 0., -10;

    for (auto& gridVec : grids) {
        for (auto& grid : gridVec) {
            grid -> force = gravity * grid->mass;
        }
    }

    for (auto& particle : particles) {   
        double JE = particle->FE.determinant();
        double JP = particle->FP.determinant();
        double lambda = material.lambda * exp(material.xsi * (1 - JP));
        double mu = material.mu * exp(material.xsi * (1 - JP));
        // std::cout << "JE: " << JE << std::endl;
        // std::cout << "JP: " << JP << std::endl;
        // std::cout << "lambda: " << lambda << std::endl;
        // std::cout << "mu: " << mu << std::endl;

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(particle->FE, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d RE = U * V.transpose();
        Eigen::Matrix2d Sig = U.inverse() * particle->FE * V.transpose().inverse();
        Eigen::Matrix2d S = V * Sig * V.transpose();
        // std::cout << "U" << std::endl << U << std::endl;
        // std::cout << "V" << std::endl << V << std::endl;
        // std::cout << "S" << std::endl << S << std::endl;
        // std::cout << "R" << std::endl << RE << std::endl;
        // std::cout << "SR" << std::endl << RE * S << std::endl << particle->FE << std::endl;

        Eigen::Matrix2d PF = 2 * mu * (particle->FE - RE) + lambda * (JE - 1) * JE * particle->FE.transpose().inverse();
        // std::cout << "F-R" << std::endl << particle->FE - RE << std::endl;
        // std::cout << "JE*:" << (JE-1) * JE << std::endl; 
        // std::cout << "PF:" << std::endl << PF << std::endl; 

        for (int dx = -1; dx <= 2; dx++) {
            for (int dy = -1; dy <= 2; dy++) {
                int xGrid = particle->xLeft + dx;
                int yGrid = particle->yLeft + dy;
                if (xGrid >= 0 && xGrid < nGrid && yGrid >= 0 && yGrid < nGrid) {
                    Eigen::Vector2d weightGradient;
                    weightGradient << particle->xWeightGradient[dx+1] * particle->yWeight[dy+1], 
                                      particle->yWeightGradient[dy+1] * particle->xWeight[dx+1];
                    // std::cout << "gradW:" << std::endl << weightGradient << std::endl; 
                    grids[xGrid][yGrid] -> force -= particle->volume * PF * particle->FE.transpose() * weightGradient;
                }
            }
        }
    }
}

void MPM::computeGridVelocity() {
    for (auto& gridVec : grids) {
        for (auto& grid : gridVec) {
            // std::cout << grid->force << std::endl << std::endl;
            grid -> newLinearMomentum = grid -> linearMomentum + timeStep * grid->force;
        }
    }
}

void MPM::updateDeformation() {
    for (auto& particle : particles) {   
        Eigen::Matrix2d velocityGradient = Eigen::Matrix2d::Zero();

        for (int dx = -1; dx <= 2; dx++) {
            for (int dy = -1; dy <= 2; dy++) {
                int xGrid = particle->xLeft + dx;
                int yGrid = particle->yLeft + dy;
                if (xGrid >= 0 && xGrid < nGrid && yGrid >= 0 && yGrid < nGrid) {
                    Eigen::Vector2d weightGradient;
                    weightGradient << particle->xWeightGradient[dx+1] * particle->yWeight[dy+1], 
                                      particle->yWeightGradient[dy+1] * particle->xWeight[dx+1];
                    velocityGradient += grids[xGrid][yGrid]->newLinearMomentum / grids[xGrid][yGrid]->mass * weightGradient.transpose();
                }
            }
        }


        // std::cout << "vgrad: " << std::endl << velocityGradient << std::endl;
        Eigen::Matrix2d Fnew = (Eigen::Matrix2d::Identity() + timeStep * velocityGradient) * particle->FE;
        particle->FE = Fnew;
        
        // std::cout << "Fnew: " << std::endl << Fnew << std::endl;

        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(Fnew, Eigen::ComputeThinU | Eigen::ComputeThinV);
        // Eigen::Matrix2d U = svd.matrixU();
        // Eigen::Matrix2d V = svd.matrixV();
        // Eigen::Matrix2d S = U.inverse() * Fnew * V.transpose().inverse();
        
        
        // std::cout << "S: " << std::endl << S << std::endl;
        // for (int i = 0; i < 2; i++) {
        //     if (S(i, i) > 1 + material.thetaS) S(i, i) = 1 + material.thetaS;
        //     if (S(i, i) < 1 - material.thetaC) S(i, i) = 1 - material.thetaC;
        // }

        // particle->FE = U * S * V.transpose();
        // particle->FP = (V * S.inverse() * U.transpose() * Fnew) * particle->FP;
        
        // std::cout << "FE: " << std::endl << particle->FE << std::endl;
        // std::cout << "FP: " << std::endl << particle->FP << std::endl;

    }

}

void MPM::updateParticleVelocity() {
    double alpha = 0.;
    for (auto& particle : particles) {   
        Eigen::Vector2d vPIC = Eigen::Vector2d::Zero();
        // Eigen::Vector2d vFLIP = particle->velocity;
        particle->B = Eigen::Matrix2d::Zero();

        for (int dx = -1; dx <= 2; dx++) {
            for (int dy = -1; dy <= 2; dy++) {
                int xGrid = particle->xLeft + dx;
                int yGrid = particle->yLeft + dy;
                Eigen::Vector2d displacement;
                displacement << (particle->xDiff - dx) / nGrid, (particle->yDiff - dy) / nGrid;
                if (xGrid >= 0 && xGrid < nGrid && yGrid >= 0 && yGrid < nGrid) {
                    Eigen::Vector2d gridV = grids[xGrid][yGrid]->linearMomentum / grids[xGrid][yGrid] -> mass;
                    Eigen::Vector2d gridNewV = grids[xGrid][yGrid]->newLinearMomentum / grids[xGrid][yGrid] -> mass;
                    vPIC += particle->xWeight[dx+1] * particle->yWeight[dy+1] * gridNewV;
                    // vFLIP += particle->xWeight[dx+1] * particle->yWeight[dy+1] * (gridNewV - gridV);
                    particle->B += particle->xWeight[dx+1] * particle->yWeight[dy+1] * gridNewV * displacement.transpose(); 
                }
            }
        }

        // particle->velocity = (1-alpha) * vPIC + alpha * vFLIP;
        particle->velocity = vPIC;
    }
}

void MPM::handleParticleCollision() {
    for (auto& particle : particles) {
        if (particle -> position[1] <= 0.1 && particle -> velocity[1] <= 0.) {
            particle -> velocity[1] = 0;
        }
    }
}

void MPM::updateParticlePosition() {
    for (auto& particle : particles) {
        particle -> position += timeStep * particle -> velocity;
    }
}
