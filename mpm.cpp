#include "mpm.hpp"


MPM::MPM(int nGrid, double timeStep) {
    this->nGrid = nGrid;
    this->time = 0;
    this->timeStep = timeStep;
    
    MaterialParameters* material = new MaterialParameters(5000, 0.2, 1, 0.005, 0.0015, 1000); 
    materials.push_back(material);

    addCube(Eigen::Vector2d(0.5, 0.25), Eigen::Vector2d(0.1, 0.3), 0, 17, 1, materials[0], Eigen::Vector3d(1, 1, 1));
    // addCube(Eigen::Vector2d(0.52, 0.2), Eigen::Vector2d(0.1, 0.1), -0.1, 10, 1, materials[0], Eigen::Vector3d(0.5, 0.5, 1));
    // addCube(Eigen::Vector2d(0.3, 0.8), Eigen::Vector2d(0.15, 0.15), 0.2, 12, 1, materials[0], Eigen::Vector3d(0.5, 1, 0.5));

    grids = std::vector<std::vector<Grid*>>(nGrid, std::vector<Grid*>(nGrid));
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j] = new Grid();
            grids[i][j]->position << (((double)i + 0.5) / nGrid), (((double)j + 0.5) / nGrid);
            grids[i][j]->linearMomentum << 0, 0;
            grids[i][j]->color = Eigen::Vector3d(0.5, 0, 0);
        }
    }
}


void MPM::addCube(Eigen::Vector2d position, Eigen::Vector2d size, double angle, 
                    int div, int random, 
                    MaterialParameters* material, Eigen::Vector3d color)
{
    for (int i = 0; i < div; i++) {
        for (int j = 0; j < div; j++) {
            Particle* p = new Particle();
            p -> mass = material->rho * size(0) * size(1) / (double)div / (double)div;
            double x, y;
            if (random == 1) {
                double rx = rand() % 10000 / 10000. - 0.5;
                double ry = rand() % 10000 / 10000. - 0.5;
                x = position(0) + cos(angle) * size(0) * rx + sin(angle) * size(1) * ry;
                y = position(1) + cos(angle) * size(1) * ry - sin(angle) * size(0) * rx;
            } else {
                x = position(0) + cos(angle) * size(0) * (1./((double)div-1.)*(double)i-0.5) + sin(angle) * size(1) * (1./((double)div-1.)*(double)j-0.5);
                y = position(1) + cos(angle) * size(1) * (1./((double)div-1.)*(double)j-0.5) - sin(angle) * size(0) * (1./((double)div-1.)*(double)i-0.5);
            }
            p -> position << x, y;
            p -> velocity << 0., 0;
            p -> material = material;
            p -> color = color;
            particles.push_back(p);
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
    for (auto& material : materials) {
        delete material;
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
    glColor3f(1, 1, 0);
    for (auto& particle : particles) {
        particle->render();
    }
    
    // Rendering grids
    glColor3f(0.5, 0, 0);
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
        // particle->xLeft = (int) (particle->position[0] * nGrid - 0.5);
        // particle->yLeft = (int) (particle->position[1] * nGrid - 0.5);
        particle->xLeft = floor(particle->position[0] * nGrid);
        particle->yLeft = floor(particle->position[1] * nGrid);
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
        // double lambda = particle->material->lambda * exp(particle->material->xsi * (1 - JP));
        // double mu = particle->material->mu * exp(particle->material->xsi * (1 - JP));
        
        // Hardening removed 
        double lambda = particle->material->lambda;
        double mu = particle->material->mu;

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(particle->FE, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d RE = U * V.transpose();
        Eigen::Matrix2d Sig = U.inverse() * particle->FE * V.transpose().inverse();
        Eigen::Matrix2d S = V * Sig * V.transpose();

        // Eigen::Matrix2d PF = 2 * mu * (particle->FE - RE) + lambda * (JE - 1) * JE * particle->FE.transpose().inverse();
        
        // Eigen::Matrix2d PF = (0.2 * 1000 * (JP - 1)) * particle->FE.transpose().inverse();
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
    double friction_mu = 3.;
    for (auto& gridVec : grids) {
        for (auto& grid : gridVec) {
            // std::cout << grid->force << std::endl << std::endl;
            grid -> newLinearMomentum = grid -> linearMomentum + timeStep * grid->force;
            if (grid -> position(1) <= 0.1) {
                if (abs(grid -> newLinearMomentum(0)) < friction_mu * abs(grid -> newLinearMomentum(1))) {
                    grid -> newLinearMomentum(1) = 0;
                    grid -> newLinearMomentum(0) = 0;
                    grid -> color = Eigen::Vector3d(0.2, 0.8, 0.8);
                } else {
                    grid -> newLinearMomentum(0) = grid -> newLinearMomentum(0) / abs(grid -> newLinearMomentum(0)) * 
                        (abs(grid -> newLinearMomentum(0)) - friction_mu * abs(grid -> newLinearMomentum(1)));
                    grid -> newLinearMomentum(1) = 0;
                    grid -> color = Eigen::Vector3d(0.8, 0.2, 0.8);
                }
            } else {
                grid -> color = Eigen::Vector3d(0.2, 0.2, 0.2);
            }
        }
    }
}

void MPM::updateDeformation() {
    double maxS = -10000;
    double minS = 10000;

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
        // particle->FE = Fnew;
        
        // std::cout << "Fnew: " << std::endl << Fnew << std::endl;

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Fnew, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Matrix2d U = svd.matrixU();
        Eigen::Matrix2d V = svd.matrixV();
        Eigen::Matrix2d S = U.inverse() * Fnew * V.transpose().inverse();

        double k = 5;
        Eigen::Matrix2d e;
        e << log(S(0, 0)), 0, 0, log(S(1, 1));
        
        // if (e(0, 0) + e(1, 1) >= 0) {
        //     e(0, 0) = 0;
        //     e(1, 1) = 0;
        //     particle->color = Eigen::Vector3d(1, 0, 0);
        // }
        // else if (e(1, 1) < k * e(0, 0)) {
        //     e(0, 0) = (e(0, 0) + e(1, 1)) * 1. / (k + 1.);
        //     e(1, 1) = (e(0, 0) + e(1, 1)) * k / (k + 1.);
        //     particle->color = Eigen::Vector3d(0, 0, 1);
        // }
        // else if (e(0, 0) < k * e(1, 1)) {
        //     e(0, 0) = (e(0, 0) + e(1, 1)) * k / (k + 1.);
        //     e(1, 1) = (e(0, 0) + e(1, 1)) * 1. / (k + 1.);
        //     particle->color = Eigen::Vector3d(0, 0, 1);
        // } else {
        //     particle->color = Eigen::Vector3d(0, 1, 0);
        // }

        // S << exp(e(0, 0)), 0, 0, exp(e(1, 1));

        // Sand plasticity *UPDATE*
        Eigen::Matrix2d ehat = e - e.trace() / 2. * Eigen::Matrix2d::Identity();
        double normF_ehat = sqrt((ehat * ehat).trace());
        double dgamma = normF_ehat + 
            (particle->material->lambda + particle->material->mu) / particle->material->mu * 
            e.trace() * particle->alpha;
        if (dgamma < 0) {
            particle->color = Eigen::Vector3d(0, 1, 0);
        } 
        else if (normF_ehat == 0 || e.trace()>0) {
            S = Eigen::Matrix2d::Identity();
            particle->q += sqrt((e * e).trace());
            particle->color = Eigen::Vector3d(1, 0, 0);
        } 
        else {
            Eigen::Matrix2d H = e - dgamma * ehat / normF_ehat;
            S << exp(H(0, 0)), 0, 0, exp(H(1, 1));
            particle->q += dgamma;
            particle->color = Eigen::Vector3d(0, 0, 1);
        }
        particle->phi = particle->h0 + (particle->h1 * particle->q - particle->h3) * 
            exp(-particle->h2 * particle->q);
        particle->alpha = sqrt(2./3.) * 2 * sin(particle->phi*3.1416/180) / (3 - sin(particle->phi*3.1416/180));

        particle->FE = U * S * V.transpose();
        particle->FP = (V * S.inverse() * U.transpose() * Fnew) * particle->FP;
        
        // particle->FE = Eigen::Matrix2d::Identity();
        // particle->FP = Fnew * particle->FP;

        // std::cout << "FE: " << std::endl << particle->FE << std::endl;
        // std::cout << "FP: " << std::endl << particle->FP << std::endl;

    }
    
    // std::cout << minS << " ~ " << maxS << std::endl;

}

void MPM::updateParticleVelocity() {
    double alpha = 0.;
    for (auto& particle : particles) {   
        Eigen::Vector2d vPIC = Eigen::Vector2d::Zero();
        Eigen::Vector2d vFLIP = particle->velocity;
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
                    vFLIP += particle->xWeight[dx+1] * particle->yWeight[dy+1] * (gridNewV - gridV);
                    particle->B += particle->xWeight[dx+1] * particle->yWeight[dy+1] * gridNewV * displacement.transpose(); 
                }
            }
        }

        particle->velocity = (1-alpha) * vPIC + alpha * vFLIP;
        // particle->velocity = vPIC;
    }
}

void MPM::handleParticleCollision() {
    // for (auto& particle : particles) {
    //     if (particle -> position[1] <= 0.1 && particle -> velocity[1] <= 0.) {
    //         particle -> velocity[1] = 0;
    //         particle -> velocity[0] *= 0.99;
    //     }
    // }
}

void MPM::updateParticlePosition() {
    for (auto& particle : particles) {
        particle -> position += timeStep * particle -> velocity;
    }
}
