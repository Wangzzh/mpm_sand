#include "mpm.hpp"


MPM::MPM(int nGrid, double timeStep) {
    this->nGrid = nGrid;
    this->time = 0;
    this->timeStep = timeStep;
    
    // MaterialParameters material = MaterialParameters(1000, 0.2, 10, 0.025, 0.0075, 1000); // elastic
    MaterialParameters* material = new MaterialParameters(10000, 0.2, 1, 0.005, 0.0015, 1000); // collapsing
    materials.push_back(material);

    addCube(Eigen::Vector2d(0.5, 0.5), Eigen::Vector2d(0.1, 0.1), 0.4, 12, 1, materials[0], Eigen::Vector3d(1, 1, 1));
    addCube(Eigen::Vector2d(0.52, 0.2), Eigen::Vector2d(0.1, 0.1), -0.1, 12, 1, materials[0], Eigen::Vector3d(1, 1, 0.5));

    grids = std::vector<std::vector<Grid*>>(nGrid, std::vector<Grid*>(nGrid));
    for (int i = 0; i < nGrid; i++) {
        for (int j = 0; j < nGrid; j++) {
            grids[i][j] = new Grid();
            grids[i][j]->position << (((double)i + 0.5) / nGrid), (((double)j + 0.5) / nGrid);
            grids[i][j]->linearMomentum << i * 0.1, j * 0.1;
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
    glColor3f(1, 0, 0);
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
        double lambda = particle->material->lambda * exp(particle->material->xsi * (1 - JP));
        double mu = particle->material->mu * exp(particle->material->xsi * (1 - JP));
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
        //2 * mu * (particle->FE - RE) + lambda * (JE - 1) * JE * particle->FE.transpose().inverse();

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
            if (grid -> position(1) < 0.1 && grid -> newLinearMomentum(1) < 0.) {
                grid -> newLinearMomentum(1) = 0;
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
        
        
        // std::cout << "S: " << std::endl << S << std::endl;
        for (int i = 0; i < 2; i++) {
            if (S(i, i) > maxS) maxS = S(i, i);
            if (S(i, i) < minS) minS = S(i, i);
            if (S(i, i) > 1 + particle->material->thetaS) S(i, i) = 1 + particle->material->thetaS;
            if (S(i, i) < 1 - particle->material->thetaC) S(i, i) = 1 - particle->material->thetaC;
        }


        particle->FE = U * S * V.transpose();
        particle->FP = (V * S.inverse() * U.transpose() * Fnew) * particle->FP;
        
        // std::cout << "FE: " << std::endl << particle->FE << std::endl;
        // std::cout << "FP: " << std::endl << particle->FP << std::endl;

    }
    
    // std::cout << minS << " ~ " << maxS << std::endl;

}

void MPM::updateParticleVelocity() {
    double alpha = 0.05;
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
    //         particle -> velocity[0] *= 0.999;
    //     }
    // }
}

void MPM::updateParticlePosition() {
    for (auto& particle : particles) {
        particle -> position += timeStep * particle -> velocity;
    }
}
