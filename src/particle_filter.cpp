/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;
	normal_distribution<double> nx(x, std[0]);
	normal_distribution<double> ny(y, std[1]);
	normal_distribution<double> ntheta(theta, std[2]);
	num_particles = 10;
	for (int i = 0; i < num_particles; i++) {
		Particle particleTemp;
		particleTemp.x = nx(gen);
		particleTemp.y = ny(gen);
		particleTemp.theta = ntheta(gen);
		particleTemp.id = i;
		particleTemp.weight = 1.0;
		particles.push_back(particleTemp);
	}
	is_initialized = true;

	//cout << "Init Done" << endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> nx(0, std_pos[0]);
	normal_distribution<double> ny(0, std_pos[1]);
	normal_distribution<double> ntheta(0, std_pos[2]);

	//cout << "Prediction start" << endl;
	for (int i = 0; i < num_particles; i++) {
		if (fabs(yaw_rate) < 0.001) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		else {
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}
		
		particles[i].x += nx(gen);
		particles[i].y += ny(gen);
		particles[i].theta += ntheta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//cout << "Data Association " << endl;
	for (int i = 0; i < observations.size(); i++) {
		int minDist = numeric_limits<int>::max();
		int index = -1;
		for (int j = 0; j < predicted.size(); j++) {
			int currDist = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if (currDist < minDist) {
				minDist = currDist;
				index = predicted[j].id;
			}
		}
		observations[i].id = index;
		//cout << "ObservationId:" << observations[i].id << endl;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < num_particles; i++) {

		//1. Find Landmarks in range
		vector<LandmarkObs> landmarksInRange;

		for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			
			//Boxing out landmarks out of range
			if (fabs(map_landmarks.landmark_list[j].x_f - particles[i].x) <= sensor_range && fabs(map_landmarks.landmark_list[j].y_f - particles[i].y) <= sensor_range) {
				landmarksInRange.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f });
			}
		}

		//2. Transform meas from Car to Map co-ordinates(as seen from Particle)

		vector<LandmarkObs> tObs;
		for (int j = 0; j < observations.size(); j++) {
			double tObsX = cos(particles[i].theta) * observations[j].x -
				sin(particles[i].theta) * observations[j].y + particles[i].x;
			double tObsY = sin(particles[i].theta) * observations[j].x +
				cos(particles[i].theta) * observations[j].y + particles[i].y;
			//ID is actually irrelevant here
			tObs.push_back(LandmarkObs{ observations[j].id,tObsX,tObsY });
		}

		//3. Tag the landmarks in Particles range 
		dataAssociation(landmarksInRange, tObs);

		particles[i].weight = 1.0;

		for (int j = 0; j < tObs.size(); j++) {
			double tObsX = tObs[j].x;
			double tObsY = tObs[j].y;
			double landmarkX, landmarkY;

			for (int k = 0; k < landmarksInRange.size(); k++) {
				if (tObs[j].id == landmarksInRange[k].id) {
					landmarkX = landmarksInRange[k].x;
					landmarkY = landmarksInRange[k].y;
					//break;
				}
			}


			double normalizer = 1 / (2 * M_PI *std_landmark[0] * std_landmark[1]);
			double term1 = pow(tObsX - landmarkX, 2) / (2 * pow(std_landmark[0], 2));
			double term2 = pow(tObsY - landmarkY, 2) / (2 * pow(std_landmark[1], 2));
			double weight = normalizer * exp(-(term1 + term2));
			particles[i].weight *= weight;
		}

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//cout << "Resampling" << endl;
	vector<double> weights;
	double maxWeight = numeric_limits<double>::min();

	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
		
		if (particles[i].weight > maxWeight) {
			maxWeight = particles[i].weight;
		}
	}

	//cout << "MaxWeight is " << maxWeight << endl;

	default_random_engine gen;
	uniform_real_distribution<double> randDist(0, 2*maxWeight);
	uniform_int_distribution<int> randIndex(0, num_particles-1);
	int index = randIndex(gen);
	double beta = 0;
	vector<Particle> resampled;

	for (int i = 0; i < num_particles; i++) {
		beta += randDist(gen);
		//cout << "Beta is " << beta << endl;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampled.push_back(particles[index]);
	}

	particles = resampled;	
	//cout << "Resampling done" << endl;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
