#include "SPH.h"

using namespace std;

Colors::Colors() :
    white(std::make_tuple(255, 255, 255)),
    red(std::make_tuple(255, 0, 0)),
    green(std::make_tuple(0, 255, 0)),
    blue(std::make_tuple(0, 0, 255)),
    black(std::make_tuple(0, 0, 0)) {
}


//Functions for Particle Class
Particle::Particle()
{
	
}
Particle::Particle(int t_idx, std::array<double, 2> t_location)
{
    idx = t_idx;
    attractionRadius = 100;
    
    dt = 0.001;
    damping = 1;
    radius = 3;
    domainSize = 1000;
    color = sf::Color::Blue;
    location = t_location;
    mass = 2;
    density = 1;
    pressure = 0;
    velocity = {0, 0};
    force = {0, 0};
    rhoZero = 8 * 10E-7;
    gamma = 3;
    AttractionCoefficient = 0.00007;
    p0 = 0;
    K = 0.006;
    inertiaFactor = 0.5;
}

void Particle::detectCorners() {
    if (location[0] < 0 + radius) {
        location[0] = 0 + radius;
        velocity[0] = -damping * velocity[0];
    }

    if (location[0] > domainSize - radius) {
        location[0] = domainSize - radius;
        velocity[0] = -damping * velocity[0];
    }

    if (location[1] < 0 + radius) {
        location[1] = 0 + radius;
        velocity[1] = -damping * velocity[1];
    }

    if (location[1] > domainSize - radius) {
        location[1] = domainSize - radius;
        velocity[1] = -damping * velocity[1];
    }
}


void Particle::detectNeighbours(std::vector<std::vector<std::vector<int>  > > grid) {
	int x, y;
	try {
		x = static_cast<int>(location[0] / attractionRadius);
		y = static_cast<int>(location[1] / attractionRadius);

		int ulim = grid.size();

		if (x >= ulim) x = ulim - 1;
		if (y >= ulim) y = ulim - 1;
		if (x <= 0) x = 0;
		if (y <= 0) y = 0;

		neighbours = grid[x][y]; // my grid

		if (x < ulim - 1) 
		{
			neighbours.insert(neighbours.end(), grid[x + 1][y].begin(), grid[x + 1][y].end()); // right
			if (y < ulim - 1) 
			{
				neighbours.insert(neighbours.end(), grid[x + 1][y + 1].begin(), grid[x + 1][y + 1].end()); // top right
			}
			if (y > 0) 
			{
				neighbours.insert(neighbours.end(), grid[x + 1][y - 1].begin(), grid[x + 1][y - 1].end()); // bottom right
			}
		}

		if (x > 0) 
		{
			neighbours.insert(neighbours.end(), grid[x - 1][y].begin(), grid[x - 1][y].end()); // left
			if (y < ulim - 1) 
			{
				neighbours.insert(neighbours.end(), grid[x - 1][y + 1].begin(), grid[x - 1][y + 1].end()); // top left
			}
			if (y > 0) 
			{
				neighbours.insert(neighbours.end(), grid[x - 1][y - 1].begin(), grid[x - 1][y - 1].end()); // bottom left
			}
		}

		if (y > 0) 
		{
			neighbours.insert(neighbours.end(), grid[x][y - 1].begin(), grid[x][y - 1].end()); // bottom
		}
		if (y < ulim - 1) 
		{
			neighbours.insert(neighbours.end(), grid[x][y + 1].begin(), grid[x][y + 1].end()); // top
		}
		/*for(int kk = 0; kk < neighbours.size(); kk++)
		{
			cout << neighbours[kk] << ",";
		}
		cout << endl << endl;*/
	} 
	catch (...) 
	{
		std::cout << "Error in finding neighbours" << std::endl;
		std::cout << x << " " << y << " " << location[0] << " " << location[1] << std::endl;
	}
}
	
	
double Particle::KernalFunction(double x) 
{
	double h = attractionRadius;
	return 6 * (h - x) * (h - x) / (M_PI * h * h * h * h);
}

double Particle::KernalDerivative(double x) 
{
	double h = attractionRadius;
	return 12 * (x - h) / (M_PI* h * h * h * h);
}

void Particle::calculateProperties(std::vector<Particle>& particleArray) 
{
	pressure = 0;
	density = 0;
	// F = std::array<double, 2> {0, 0};

	for (int i : neighbours) {
		Particle& p = particleArray[i];
		std::array<double, 2> dVector = {p.location[0] - location[0], p.location[1] - location[1]};
		double d = sqrt(dVector[0] * dVector[0] + dVector[1] * dVector[1]);
		if (d < attractionRadius && d > 0) {
			density += p.mass * KernalFunction(d);
		}
	}

	if (density == 0) {
		density = rhoZero * 0.001;
	}

	pressure = (pow((density /rhoZero), gamma)-1) * AttractionCoefficient;
}
std::array<double,2> Particle::calculatePressureForce(double r, Particle P) 
{
    double pF = mass * (pressure / pow(density, 2) + P.pressure / pow(P.density, 2)) * KernalDerivative(r);
    //cout << "pressure Force = " << pF << endl;
    return {pF , pF };
}

std::array<double, 2> Particle::calculateViscousForce(double r, Particle P) 
{
	double k =  KernalDerivative(r);
    return {(velocity[0] - P.velocity[0]) *k, (velocity[1] - P.velocity[1]) *k};
}

std::array<double, 2> Particle::calculateBodyForces(double r) 
{
	//if (location[1] > 700 and location[0] > 700)
	//{
		//return {-100, -100};
	//}
    return {0,1000};
}

void Particle::calculateNetForce(std::vector<Particle>& particleArray) 
{
	std::array<double, 2> F = {0.0, 0.0};
	for (int i : neighbours) {
		Particle& p = particleArray[i];
		std::vector<double> dVector(2);
		for (int j = 0; j < 2; ++j)
			dVector[j] = p.location[j] - location[j];
		double d = std::sqrt(dVector[0] * dVector[0] + dVector[1] * dVector[1]);
		if (d < attractionRadius && d > 0) {
			std::array<double, 2> pressureForce = calculatePressureForce(d, p);
			for (int j = 0; j < 2; ++j)
				F[j] += pressureForce[j] * dVector[j] / d;
			std::array<double, 2> viscousForce = calculateViscousForce(d, p);
			for (int j = 0; j < 2; ++j)
				F[j] += viscousForce[j];
		}
	}
	std::array<double, 2> bodyForce = calculateBodyForces(10);
	for (int j = 0; j < 2; ++j)
		F[j] += bodyForce[j];
	force = F;
	//cout << "Index = " << idx << " Net force = " << force[0] << " : " << force[1] << endl;
}

void Particle::move() 
{
	std::array<double, 2> acceleration;
	std::array<double, 2> newVelocity;
	std::array<double, 2> newLocation;
	for (int i = 0; i < 2; ++i) {
		acceleration[i] = force[i] / density;
		newVelocity[i] = inertiaFactor * velocity[i] + acceleration[i] * dt;
		newLocation[i] = location[i] + newVelocity[i] * dt;
	}

	// Update the particle's velocity and location
	velocity = newVelocity;
	location = newLocation;

	// Call detectCorners() if it's a member function of Particle class
	detectCorners();
}

//--- End of functions for Particle Class ----

//----------------------------------------------------------------------------

// --- Functions for the World Class

World::World() 
{
	particleCount = 2000;
    dt = 0.001;
}
void World::drawParticle(sf::RenderWindow& window, Particle& p) 
{
    sf::CircleShape circle(p.radius);
    circle.setFillColor(p.color);
    circle.setPosition(p.location[0] - p.radius, p.location[1] - p.radius); // Adjust position to center
    window.draw(circle);
}

void World::generateParticleArray(int count) 
{
    int margin = 300;
    int side = static_cast<int>(std::sqrt(count));
    double x = margin;
    double y = margin;
    double delta = (width - 2 * margin) / side;

    for (int i = 0; i < count; ++i) {
        x += delta;
        if (x > margin + side * delta) {
            y += delta;
            x = margin;
        }
        Particle p(i, {x, y}); // Assuming Particle constructor takes id and location
        particleArray.push_back(p);
    }
    
    //Generate domain grid template to map the particles
    Particle pp = Particle(0, {0,0});
    int x_size = static_cast<int>(width/pp.attractionRadius);
    int y_size = static_cast<int>(height/pp.attractionRadius);
    for(size_t i = 0; i < x_size; i++)
    {
		std::vector<vector<int> > row = {};
		for(size_t i = 0; i < x_size; i++)
		{
			std::vector<int> a = {};
			row.push_back(a);
		}
		domainGrid.push_back(row);
	}
}

void World::mapParticleToArray(Particle& p) 
{  
    int ulim = domainGrid.size();
    
    int x = static_cast<int>(p.location[0] / p.attractionRadius);
    int y = static_cast<int>(p.location[1] / p.attractionRadius);

    if (x >= ulim)
        x = ulim - 1;
    if (y >= ulim)
        y = ulim - 1;
    if (x <= 0)
        x = 0;
    if (y <= 0)
        y = 0;

    domainGrid[x][y].push_back(p.idx);
}

void World::clearGrid() 
{
    for (size_t i = 0; i < domainGrid.size(); ++i) {
        for (size_t j = 0; j < domainGrid[i].size(); ++j) {
            domainGrid[i][j].clear();
        }
    }
}


void World::runStep(sf::RenderWindow& window) 
{
    // Detect neighbours and calculate properties for each particle
    for (Particle& p : particleArray) 
    {
        p.detectNeighbours(domainGrid);
        p.calculateProperties(particleArray);
        //cout << p.idx << "  " << p.density << endl;
    }

    // Calculate net force for each particle
    for (Particle& p : particleArray) 
    {
        p.calculateNetForce(particleArray);
    }

    // Clear the grid
    clearGrid();

    // Move particles, map to array, and draw
    for (Particle& p : particleArray) 
    {
        p.move();
        mapParticleToArray(p);
        drawParticle(window, p);
    }
}

void World::addDietoFluid()
{
	for(size_t i = 0; i < particleArray.size(); i++)
	{
		if(particleArray[i].location[0] < 650 && particleArray[i].location[0] > 350)
		{
			if(particleArray[i].location[1] < 700 && particleArray[i].location[1] > 500)
			{
				particleArray[i].mass *= 0.8;
				particleArray[i].color = sf::Color::Red;
				//particleArray[i].velocity = {1000, 1000};
			}
		}
	}
}

void World::addParticles(double x, double y)
{
	for(int i = 0; i < 50; i++)
	{
		double random1 = rand() % 200;
		double random2 = rand() % 200;
		Particle p = Particle(particleArray.size(), {x+random1,y+random2});
		particleArray.push_back(p);
	}	
}


void World::run()
{
	width = height = 1000;
	sf::RenderWindow window(sf::VideoMode(width, height), "My SFML Application");
	window.clear(sf::Color::Black);

    window.setFramerateLimit(240.f);
    
	generateParticleArray(particleCount);
	
	for(size_t i = 0; i < particleCount; i++)
	{
		//cout << i << "  " << particleArray[i].location[0] << "   " << particleArray[i].location[1] << endl;
		mapParticleToArray(particleArray[i]);
	}
	int frame = 0;

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
                
            if (event.type == sf::Event::KeyReleased && event.key.code == sf::Keyboard::Space) 
            {
                cout << "Spacebar pressed" << endl;
                addDietoFluid();
            }
            if (event.type == sf::Event::MouseButtonPressed)
			{
				if (event.mouseButton.button == sf::Mouse::Left)
				{
					std::cout << "the right button was pressed" << std::endl;
					double x = event.mouseButton.x;
					double y = event.mouseButton.y;
					addParticles(x, y);	
				}  
			}             
        }
        window.clear(sf::Color::Black);
        runStep(window);
        window.display();
        //cout << "Frame : " << frame << endl;
        frame++;
    }
}

// --- End of Functions for the World Class ---
int main()
{
	World w = World();
	w.run();
	return 0;
	
}

