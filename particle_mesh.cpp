#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <map>
#include <chrono>

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

public:
    P3MSimulator(double dt) : timeStep(dt) {}

    void addParticle(const Particle& particle) {
        particles.push_back(particle);
    }

    void computeForces() {
        // Reset forces
        for (auto& particle : particles) {
            std::fill(particle.force.begin(), particle.force.end(), 0.0);
        }

        // Calculate forces (simple direct method for illustration)
        for (size_t i = 0; i < particles.size(); ++i) {
            for (size_t j = i + 1; j < particles.size(); ++j) {
                std::vector<double> diff(3);
                double rSquared = 0.0;
                for (int k = 0; k < 3; ++k) {
                    diff[k] = particles[i].position[k] - particles[j].position[k];
                    rSquared += diff[k] * diff[k];
                }
                double r = sqrt(rSquared);
                double forceMagnitude = kCoulomb * particles[i].charge * particles[j].charge / rSquared;

                for (int k = 0; k < 3; ++k) {
                    double f = forceMagnitude * diff[k] / r;
                    particles[i].force[k] += f;
                    particles[j].force[k] -= f; // Newton's third law
                }
            }
        }
    }

    void integrate() {
        for (auto& particle : particles) {
            for (int k = 0; k < 3; ++k) {
                particle.velocity[k] += timeStep * particle.force[k] / particle.mass;
                particle.position[k] += timeStep * particle.velocity[k];
            }
        }
    }

    std::vector<std::tuple<double, double, double>> printParticles() {
        std::vector<std::tuple<double, double, double>> coordinates;
        for (const auto& p : particles) {
            std::cout << "Position: (" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << ") ";
            std::cout << "Velocity: (" << p.velocity[0] << ", " << p.velocity[1] << ", " << p.velocity[2] << ")\n";
            coordinates.push_back({p.position[0], p.position[1], p.position[2]});
        }
        return coordinates;
    }
};

int main() {
    P3MSimulator sim(0.001); // timestep of 0.01s
    int iterations = 10;
    sim.addParticle(Particle(1.0, 0.01, {0.0, 0.0, 0-.0}, {0.0, 0.0, 0.5}));
    sim.addParticle(Particle(0.5, -0.01, {1.0, 0.0, 0.5}, {0.0, 0.0, 0.0}));
    sim.addParticle(Particle(2.0, 0.01, {0.5, 0.5, 0.0}, {0.1, 0.0, 0.0}));
    sim.addParticle(Particle(0.75, -0.01, {0.5, -0.5, 0.0}, {0.0, 0.1, 0.0}));
    sim.addParticle(Particle(0.25, 0.02, {-0.5, 0.5, 0.0}, {-0.1, -0.1, 0.0}));
    sim.addParticle(Particle(0.5, 0.02, {-1.0, -1.0, 0.0}, {0.05, 0.05, 0.0}));
    sim.addParticle(Particle(1.2, 0.015, {0.3, 0.2, 0.1}, {0.05, 0.02, 0.03}));
    sim.addParticle(Particle(0.8, -0.020, {-0.4, 0.1, 0.2}, {0.0, -0.03, 0.01}));
    sim.addParticle(Particle(0.9, 0.010, {0.7, -0.6, 0.3}, {-0.02, 0.04, -0.05}));
    sim.addParticle(Particle(0.7, -0.005, {0.2, 0.8, -0.1}, {0.01, -0.01, 0.02}));
    sim.addParticle(Particle(1.1, 0.025, {-0.3, -0.2, 0.4}, {0.03, 0.05, 0.0}));
    sim.addParticle(Particle(0.6, 0.018, {0.5, -0.4, 0.2}, {0.04, 0.01, -0.03}));
    sim.addParticle(Particle(0.5, -0.012, {-0.5, 0.5, -0.2}, {-0.03, 0.0, 0.02}));
    sim.addParticle(Particle(1.0, 0.020, {0.1, -0.3, 0.5}, {0.02, -0.02, -0.01}));
    sim.addParticle(Particle(0.4, -0.008, {-0.2, 0.4, 0.1}, {0.01, 0.03, -0.04}));


    std::vector<std::tuple<double, double, double>> temp_coordinates;
    std::map<std::string, std::vector<std::tuple<double, double, double>>> data;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < iterations; ++i) {
        sim.computeForces();
        sim.integrate();
        temp_coordinates = sim.printParticles();
        std::string key = "iteration_" + std::to_string(i);
        data[key] = temp_coordinates;
        std::ofstream file("output.txt");
        if (file.is_open()) {
            for (const auto& pair : data) {
                file << pair.first << ": "; // Write the key
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
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << "s\n";
}
