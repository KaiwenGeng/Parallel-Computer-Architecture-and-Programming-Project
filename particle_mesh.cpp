#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <chrono>
#include <random>
#include <omp.h>
class Particle {
public:
    double mass;
    double charge;
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> force;

    Particle(double m, double q, std::vector<double> pos, std::vector<double> vel)
        : mass(m), charge(q), position(pos), velocity(vel), force({0.0, 0.0, 0.0}) {}
};

class P3MSimulator {
private:
    std::vector<Particle> particles;
    double timeStep;
    double kCoulomb = 8.9875517873681764e9; // Coulomb's constant in N*m^2/C^2
    double G = 6.67430e-11; // Gravitational constant in m^3/kg/s^2

public:
    P3MSimulator(double dt) : timeStep(dt) {}

    void addParticle(const Particle& particle) {
        particles.push_back(particle);
    }

    void computeForces() {
        // Reset forces
        #pragma omp parallel for 
        for (auto& particle : particles) {
            std::fill(particle.force.begin(), particle.force.end(), 0.0);
        }
        #pragma omp parallel for 

        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                std::vector<double> diff(3);
                double rSquared = 0.0;
                for (int k = 0; k < 3; ++k) {
                    diff[k] = particles[i].position[k] - particles[j].position[k];
                    rSquared += diff[k] * diff[k];
                }
                double r = sqrt(rSquared);
                double forceMagnitudeCoulomb  = kCoulomb * particles[i].charge * particles[j].charge / rSquared;
                double forceMagnitudeGravitational = G * particles[i].mass * particles[j].mass / rSquared;
                for (int k = 0; k < 3; ++k) {
                    double fCoulomb = forceMagnitudeCoulomb * diff[k] / r;
                    double fGravitational = forceMagnitudeGravitational * diff[k] / r;
                    particles[i].force[k] += (fCoulomb - fGravitational); // Coulomb force is repulsive or attractive, gravitational force is always attractive
                    particles[j].force[k] -= (fCoulomb - fGravitational); // Newton's third law
                }
            }
        }
    }

    void integrate() {
        #pragma omp parallel for 
        for (auto& particle : particles) {
            for (int k = 0; k < 3; ++k) {
                particle.velocity[k] += timeStep * particle.force[k] / particle.mass;
                particle.position[k] += timeStep * particle.velocity[k];
            }
        }
    }
    void generateRandomParticles(int n) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> massDist(0.5, 2);  
        std::uniform_real_distribution<> chargeDist(-0.05, 0.05);  
        std::uniform_real_distribution<> posDist(-100, 100);  
        std::uniform_real_distribution<> velDist(-0.05, 0.05);  
        #pragma omp parallel for 
        for (int i = 0; i < n; ++i) {
            double mass = massDist(gen);
            double charge = chargeDist(gen);
            std::vector<double> position = {posDist(gen), posDist(gen), posDist(gen)};
            std::vector<double> velocity = {velDist(gen), velDist(gen), velDist(gen)};
            addParticle(Particle(mass, charge, position, velocity));
        }
    }

    std::vector<std::tuple<double, double, double>> get_particle_position() {
        std::vector<std::tuple<double, double, double>> coordinates;
        #pragma omp parallel for 
        for (const auto& p : particles) {
            // std::cout << "Position: (" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << ") ";
            // std::cout << "Velocity: (" << p.velocity[0] << ", " << p.velocity[1] << ", " << p.velocity[2] << ")\n";
            coordinates.push_back({p.position[0], p.position[1], p.position[2]});
        }
        return coordinates;
    }
    std::vector<std::tuple<double, double>> get_particle_weight_charge(){
        std::vector<std::tuple<double, double>> weight_charge;
        #pragma omp parallel for 
        for (const auto& p : particles) {
            weight_charge.push_back({p.mass, p.charge});
        }
        return weight_charge;
    }
};

int main() {
    P3MSimulator sim(0.0025);
    int iterations = 200;
    sim.generateRandomParticles(20);
    sim.addParticle(Particle(5, 100000, {500, 500, 500}, {0, 0, 0}));


    std::vector<std::tuple<double, double, double>> temp_coordinates;
    std::map<std::string, std::vector<std::tuple<double, double, double>>> position_data;
    std::vector<std::tuple<double, double>> temp_weight_charge;
    std::map<std::string, std::vector<std::tuple<double, double>>> weight_charge_data;
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for 
    for (int i = 0; i < iterations; ++i) {
        sim.computeForces();
        sim.integrate();
        temp_coordinates = sim.get_particle_position();
        std::string key = "iteration_" + std::to_string(i);
        position_data[key] = temp_coordinates;

        std::ofstream file("position_result.txt");
        if (file.is_open()) {
            for (const auto& pair : position_data) {
                file << pair.first << " : "; // Write the key
                for (const auto& tuple : pair.second) {
                    // Write each tuple in the vector
                    file << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ", " << std::get<2>(tuple) << ") ";
                }
                file << "\n";
            }
            file.close();
        } else {
            std::cout << "Unable to open file.";
            return 1;
        }
    }
    temp_weight_charge = sim.get_particle_weight_charge();
    std::string key2 = "weight_charge";
    weight_charge_data[key2] = temp_weight_charge;
    std::ofstream file2("weight_charge_result.txt");
    if (file2.is_open()) {
        for (const auto& pair : weight_charge_data) {
            file2 << pair.first << " : "; // Write the key
            for (const auto& tuple : pair.second) {
                // Write each tuple in the vector
                file2 << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ") ";
            }
            file2 << "\n";
        }
        file2.close();
    } else {
        std::cout << "Unable to open file.";
        return 1;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << "s\n";
}
