/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // Set the number of particles
  
  std::default_random_engine gen;
  
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2]; 
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  //Add random Gaussian noise to each particle.
  for (int i = 0; i < num_particles; i++) {
    Particle p;
	p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen); 
	p.weight = 1.0;
	particles.push_back(p);
    weights.push_back(p.weight);
  }
  /*ParticleFilter() : num_particles(0), is_initialized(false) {}*/
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   //Create Gaussian distribution for each parameter
   std::default_random_engine gen;
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);
   //Add measurements to each particle
   for(int i = 0; i < num_particles; ++i) {
	   double p_x = particles[i].x;
	   double p_y = particles[i].y;
	   double p_theta = particles[i].theta;
	   
	   double pred_x;
	   double pred_y;
	   double pred_theta;
	   
	   if (fabs(yaw_rate) < 0.0001) {
	    pred_x = p_x + velocity * cos(p_theta) * delta_t;
	    pred_y = p_y + velocity * sin(p_theta) * delta_t;
	    pred_theta = p_theta;
	  } else {
	    pred_x = p_x + (velocity/yaw_rate) * (sin(p_theta + (yaw_rate * delta_t)) - sin(p_theta));
	    pred_y = p_y + (velocity/yaw_rate) * (cos(p_theta) - cos(p_theta + (yaw_rate * delta_t)));
	    pred_theta = p_theta + (yaw_rate * delta_t);
	  }
	  

	  //add random Gaussian noise.
	  particles[i].x = pred_x + dist_x(gen);
	  particles[i].y = pred_y + dist_y(gen);
	  particles[i].theta = pred_theta + dist_theta(gen);
	   
	   
	   
   }
   
   

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

   	for (unsigned int i = 0; i < observations.size(); i++) {
	//
		double lowest_dist = std::numeric_limits<double>::max();
		int closest_landmark_id = -1;
		double obs_x = observations[i].x;
		double obs_y = observations[i].y;

		for (unsigned int j = 0; j < predicted.size(); j++) {
		  double pred_x = predicted[j].x;
		  double pred_y = predicted[j].y;
		  int pred_id = predicted[j].id;
		  double current_dist = dist(obs_x, obs_y, pred_x, pred_y);

		  if (current_dist < lowest_dist) {
		    lowest_dist = current_dist;
		    closest_landmark_id = pred_id;
		  }
		}
		observations[i].id = closest_landmark_id;

	}

   
 
   
   

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   for (int i = 0; i < num_particles; i++)
    {
       
        double px = particles[i].x;
        double py = particles[i].y;
        double ptheta = particles[i].theta;

        vector<LandmarkObs> tf_observations;

        for (unsigned int j = 0; j < observations.size(); j++)
        {
            // Transform observation to map coords, px = particles, ox = observations
            LandmarkObs tf_observation;

            double ox = observations[j].x;
            double oy = observations[j].y;

            double xm = px + (std::cos(ptheta) * ox) - (std::sin(ptheta) * oy);
            double ym = py + (std::sin(ptheta) * ox) + (std::cos(ptheta) * oy);

            tf_observation.x = xm;
            tf_observation.y = ym;
            tf_observation.id = j;

            tf_observations.push_back(tf_observation);
        }

        // Filter landmarks to only that which is within the sensor range
        vector<LandmarkObs> landmarks_best;

        for (unsigned int k = 0; k < map_landmarks.landmark_list.size(); k++)
        {
            // Calculat Distance to Landmark
            double distance = dist(px, py, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);

            // Filter the distance within the sensor range
            if (distance <= sensor_range)
            {
                LandmarkObs landmark_best;
                landmark_best.x = map_landmarks.landmark_list[k].x_f;
                landmark_best.y = map_landmarks.landmark_list[k].y_f;
                landmark_best.id = map_landmarks.landmark_list[k].id_i;

                landmarks_best.push_back(landmark_best);
            }
        }

        // Use dataAssociation to update the observations with the ID of the nearest landmark
        dataAssociation(landmarks_best, tf_observations);

        // Calculate the weight of each particle using multivariate gaussian probability density function
        double var_x = pow(std_landmark[0], 2);
        double var_y = pow(std_landmark[1], 2);
        double gaussian_norm = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);

      
        particles[i].weight = 1.0;

        double probability = 1.0;
        //loop for finding matched id
        for (unsigned int j = 0; j < tf_observations.size(); j++)
        {
            double tf_x = tf_observations[j].x;
            double tf_y = tf_observations[j].y;
            double tf_id = tf_observations[j].id;

            for (unsigned int k = 0; k < landmarks_best.size(); k++)
            {
                double ldmk_x = landmarks_best[k].x;
                double ldmk_y = landmarks_best[k].y;
                double ldmk_id = landmarks_best[k].id;

                // Find matching landmarks and observations
                if (tf_id == ldmk_id)
                {
                    double mu_x = ldmk_x;
                    double mu_y = ldmk_y;

                    // Calculate the normal probability 
                    double exponent_x = ((tf_x - mu_x) * (tf_x - mu_x)) / (2.0 * var_x);
                    double exponent_y = ((tf_y - mu_y) * (tf_y - mu_y)) / (2.0 * var_y);

                    probability = probability * gaussian_norm * exp(-(exponent_x + exponent_y));
                    break;
                }
            }
        }

        particles[i].weight = probability;
        weights[i] = probability;
    }

    // Normalise the probability to 1
    double weight_normalizer = std::accumulate(weights.begin(), weights.end(), 0.0f);

    for (int p = 0; p < num_particles; ++p)
    {
        particles[p].weight = particles[p].weight / weight_normalizer;
        weights[p] = particles[p].weight;
    }


}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   
   	vector<Particle> resampled_particles;

	
	std::default_random_engine gen;
	
	
	std::discrete_distribution<int> resample_d(weights.begin(), weights.end());
	std::vector<Particle> new_particles;

	for(int i = 0; i < num_particles; i++) {
		new_particles.push_back(particles[resample_d(gen)]);
	}

	particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}