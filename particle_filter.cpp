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
#include <random>

#include "particle_filter.h"

using namespace std;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;
	
    // This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(0.0, std[0]);
	normal_distribution<double> dist_y(0.0, std[1]);
	normal_distribution<double> dist_theta(0.0, std[2]);
	
	for (int i=0;  i<num_particles; ++i) {
		Particle p;
		p.id = i;
      	
		p.x = x + dist_x(gen) ;
        p.y = y + dist_y(gen) ;
		p.theta = theta + dist_theta(gen) ;
        p.weight =1./num_particles;
		
		particles.push_back(p);
		weights.push_back(p.weight);		
		
    }
		
	cout << "Initialized!" << endl;
	is_initialized = true;

}




void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);
		
    for (int i=0; i<num_particles; i++)	{	
		if (fabs(yaw_rate) > 0.001) {			
			particles[i].x += velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta)) +dist_x(gen);
			particles[i].y += velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta + yaw_rate*delta_t)) +dist_y(gen);
			particles[i].theta += yaw_rate*delta_t +dist_theta(gen);
		}
		else {
			particles[i].x += velocity*delta_t*cos(particles[i].theta) + dist_x(gen);
			particles[i].y += velocity*delta_t*sin(particles[i].theta) + dist_y(gen);
		}
	}
}
 
	

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	
	
	double distance_new;
	for (int i=0; i<observations.size();++i) {
		double distance_old = numeric_limits<double>::max();;
		for (int j=0; j<predicted.size(); ++j) {
			distance_new = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if (distance_new < distance_old) {
				distance_old = distance_new;
				observations[i].id = predicted[j].id;
			}
		}
	}
}
			


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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
	

	// transform observations from vehicle to map coordinate system
	for (int i=0; i<num_particles;++i)
	{

		// transform observations from vehicle coordinate system to global coordinate system
		vector<LandmarkObs> transformedObservations;
		
		for (int j=0; j < observations.size(); ++j) {
			LandmarkObs obs_global;
			obs_global.x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
			obs_global.y = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta);
			obs_global.id = observations[j].id;
			transformedObservations.push_back(obs_global);
		}
	

	
		// look for landmarks in range
		vector<LandmarkObs>  landmarksInRange;
	
		for (int j=0; j < map_landmarks.landmark_list.size(); j++){
			
			double landmark_x = map_landmarks.landmark_list[j].x_f;
			double landmark_y  = map_landmarks.landmark_list[j].y_f;
			int landmark_id  = map_landmarks.landmark_list[j].id_i;			

			if (dist(particles[i].x, particles[i].y, landmark_x, landmark_y) <= sensor_range) {
				LandmarkObs pred;
				pred.id = landmark_id;
				pred.x = landmark_x;
				pred.y = landmark_y;
				landmarksInRange.push_back(pred);			
			} 
			
		}
	
		 
		// association
		dataAssociation(landmarksInRange, transformedObservations);

		// calculating final weight by using a  multivariate gaussian, final weight is product of probabilities
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;
		double prob = 0.0;
		

		for (int j=0; j<transformedObservations.size();++j) {
			associations.push_back(transformedObservations[j].id);
			sense_x.push_back(transformedObservations[j].x);
			sense_y.push_back(transformedObservations[j].y);
		
			for (int k=0; k<landmarksInRange.size(); ++k){
				if (landmarksInRange[k].id == transformedObservations[j].id) {
					double diff_x = transformedObservations[j].x-landmarksInRange[k].x;
					double diff_y = transformedObservations[j].y-landmarksInRange[k].y;
					prob = 1/(2.*M_PI*std_landmark[0]*std_landmark[1])*exp(-(diff_x*diff_x/(2.*std_landmark[0]*std_landmark[0]) + diff_y*diff_y/(2.*std_landmark[1]*std_landmark[1]) ) );
				}	
			}				
		}
		particles[i].weight = prob;		
		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
	
	}
	
	for (int i = 0; i < num_particles; ++i) {
		weights[i] = particles[i].weight; // no need for normalization
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	

	vector<Particle> new_particles;	
	discrete_distribution<> dist_weight(weights.begin(), weights.end());
	
	
	for (int i = 0; i < num_particles; ++i) {
		int number = dist_weight(gen);
		new_particles.push_back(particles[number]);		
		}
		
	particles = new_particles;
}	
	

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
