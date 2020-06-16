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

vector<Particle> particles;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 500;  // TODO: Set the number of particles
  
  // create ramdom number generator engine
  std::default_random_engine gen;

 
  
 
   for(int i=0; i<num_particles; i++){

      // create gaussian distribution for x, y and theta
      normal_distribution<double> dist_x(x,std[0]);
      normal_distribution<double> dist_y(y, std[1]);
      normal_distribution<double> dist_theta(theta,std[2]);
    
     // sample from gaussian distribution
     double sample_x = dist_x(gen);
     double sample_y = dist_y(gen);
     double sample_theta = dist_x(gen);
    
     // intialize particle from gaussian distribution 
     Particle particle;
     particle.id = i;
     particle.x = sample_x;
     particle.y = sample_y;
     particle.theta = sample_theta;
    
    // Append each particle to the list of particles
    particles.push_back(particle);
  
    is_initialized = true;
  }
  

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
    
   // create ramdom number generator engine
   std::default_random_engine gen;
    
   for(int i=0; i<num_particles; i++) {
          // theta is zero
         if(particles[i].theta==0) {
            particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
            particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
            
            
         } 
         // theta is zero
         else {
            particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
            particles[i].y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta+yaw_rate*delta_t)-cos(particles[i].theta));
            particles[i].theta = particles[i].theta + yaw_rate*delta_t;
        }
          // Add sensor noise to the particle state
            // create gaussian distribution with zero mean
     	    normal_distribution<double> dist_x(0 ,std_pos[0]);
            normal_distribution<double> dist_y(0, std_pos[1]);
            normal_distribution<double> dist_theta(0,std_pos[2]);
            
            // sample from gaussian distribution
            double noise_x = dist_x(gen);
            double noise_y = dist_y(gen);
            double noise_theta = dist_theta(gen);
           // Add noise
            particles[i].x += noise_x;
            particles[i].y += noise_y;
            particles[i].theta += noise_theta;
         
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
    double shortest_dist = 1000.0;
    int shortest_dist_id = -1;

   for(int i=0; i<observations.size(); i++){
      for(int j=0; j<predicted.size(); j++){
          double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
          if(distance < shortest_dist){
             shortest_dist = distance;
             shortest_dist_id =predicted[j].id;
            }
      {
      
      observations[i].id = predicted[shortest_dist_id].id;
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
    for(int i=0; i<num_particles; i++) {
      
       double x_part, y_part, theta_part, x_obs, y_obs;
       int id_part;
       
       id_part = particles[i].id;
       x_part = particles[i].x;
       y_part = particles[i].y;
       theta_part = particles[i].theta; 
      
       // Map landmarks predicted to be within the particle sensor range
       vector<LandmarkObs> predicted_landmarks;
       for(int k=0; k<map_landmarks.landmark_lis.size(); k++) {
            
            double x_ldm, y_ldm;
            int id_ldm;
            
            id_ldm = map_landmarks.landmark_list[k].id;
            x_ldm = map_landmarks.landmark_list[k].x_f;
            y_ldm = map_landmarks.landmark_list[k].x_f;
            
            // Find the landmarks in the map that are within the particle sensor range
            double distance = dist(x_part, y_part, x_ldm, x_ldm);

            if(distance <= sensor_range){
             // Add landmark to the predicted landmarks
             map_landmarks.single_landmark_s.id = id_ldm;
             map_landmarks.single_landmark_s.x_f = x_ldm;
             map_landmarks.single_landmark_s.y_f = y_ldm;
             predicted_landmarks.push_back(single_landmark_s);
             }
        }
       
       
       // Transform observed landmarks from vehicle frame to map frame
       vector<LandmarkObs> trans_obs_list;
       for(int k=0; k<observations.size(); k++) {
            
            double x_obs, y_obs; 
            LandmarkObs trans_obs;

            x_obs = observations[k].x;
            y_obs = observations[k].y;
           
            trans_obs_list.x = x_part + cos(theta_part)*x_obs - sin(theta_part)*y_obs;
            trans_obs_list.y = y_part + cos(theta_part)*x_obs - sin(theta_part)*y_obs;
            trans_obs_list.id = observations[k].id;
            trans_obs_list.push_back(trans_obs);
       }
       
       
       dataAssociation(predicted_landmarks, trans_obs_list);
       
       for(int k=0; k<trans_obs_list.size(); k++){
            
            double  x = trans_obs_list[k].x; 
            double y = trans_obs_list[k].y;
            double sig_x = std_landmark[0];
            double sig_y = std_landmark[1];
            double mu_x;
            double mu_y;
            
            for (int i=0; i<predicted_landmarks.size(); i++){
                if(trans_obs_list[k].id == predicted_landmarks[i].id){
                 mu_x= predicted_landmarks[i].x;
                 mu_y = predicted_landmarks[i].y;
                 break;
                }
             }   
            
       }
      
      particles[i].weight= measurement_prob(x,y,mu_x,mu_y,ig_x,sig_y);
       
   }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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
  particle.associations= associations;
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
