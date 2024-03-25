#include <stdio.h>
#include <iostream>
#include <SFML/Graphics.hpp>
#include <tuple>
#include <vector>
#include <cmath>
#include <random>

using namespace std;



class Colors {
public:
    Colors(); // Constructor

    // Public member variables
    std::tuple<int, int, int> white;
    std::tuple<int, int, int> red;
    std::tuple<int, int, int> green;
    std::tuple<int, int, int> blue;
    std::tuple<int, int, int> black;
};
class Particle
{
	public:
	Particle();
    Particle(int t_idx, std::array<double, 2> t_location);
	void detectCorners();
	void detectNeighbours(std::vector<std::vector<std::vector<int> > > grid);
	double KernalFunction(double x);
    double KernalDerivative(double x);
    void calculateProperties(std::vector<Particle>& particleArray);
    std::array<double,2> calculatePressureForce(double r, Particle P);
	std::array<double,2> calculateViscousForce(double r, Particle P);
	std::array<double, 2> calculateBodyForces(double r);
	void calculateNetForce(std::vector<Particle>& particleArray);
	void move();
	
    int idx;
    double attractionRadius;
    double AttractionCoefficient;
    double dt;
    double damping;
    double radius;
    std::vector<int> neighbours;
    int domainSize;
    Colors C;
    sf::Color color;

    // Physics Properties
    std::array<double, 2> location;
    double mass;
    double density;
    double pressure;
    std::array<double, 2> velocity;
    std::array<double, 2> force;
    double rhoZero;
    double gamma;
    double p0;
    double K;
    double inertiaFactor;
};




class World 
{
private:
    Colors color;
    int width, height;
    std::vector<std::vector<std::vector<int> > > domainGrid;
    std::vector<Particle> particleArray;
    sf::Clock clock;
    std::vector<double> meandensity;
    std::vector<double> stdevDensity;
    std::vector<double> time;
    double dt;
    int particleCount;

    void generateParticleArray(int count);
    void mapParticleToArray(Particle& p);
    void clearGrid();
    void runStep(sf::RenderWindow& window);
    void addDietoFluid();
    void addParticles(double x, double y);

public:
    World();
    void drawParticle(sf::RenderWindow& window, Particle& p);
    void run();
		
};

