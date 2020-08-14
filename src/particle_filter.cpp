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
using std::uniform_int_distribution;
using std::uniform_real_distribution;


static std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
    
  if(is_initialized){
    return;
   }

  num_particles = 100;  // TODO: Set the number of particles
  
  // create ramdom number generator engine
  

 
   // create gaussian distribution for x, y and theta
      normal_distribution<double> dist_x(x,std[0]);
      normal_distribution<double> dist_y(y, std[1]);
      normal_distribution<double> dist_theta(theta,std[2]);
 
   for(unsigned int i=0; i<num_particles; i++){

       
     // sample from gaussian distribution
     double sample_x = dist_x(gen);
     double sample_y = dist_y(gen);
     double sample_theta = dist_theta(gen);
    
     // intialize particle from gaussian distribution 
     Particle particle;
     particle.id = i;
     particle.x = sample_x;
     particle.y = sample_y;
     particle.theta = sample_theta;
     particle.weight = 1.0;
   
    
    // Append each particle to the list of particles
    particles.push_back(particle);
    
  }
  
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
    
   // create ramdom number generator engine
   
    // create gaussian distribution with zero mean
    normal_distribution<double> dist_x(0 ,std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0,std_pos[2]);
   for(unsigned int i=0; i<num_particles; i++) {
          // theta is zero
         if(fabs(yaw_rate) <0.00001) {
            particles[i].x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
            particles[i].y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
            
            
         } 
         // theta is not zero
         else {
            particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
            particles[i].y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate*delta_t));
            particles[i].theta = particles[i].theta + yaw_rate*delta_t;
        }
          // Add sensor noise to the particle state

            
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
    

   for(unsigned int i=0; i<observations.size(); i++){

      double shortest_dist = std::numeric_limits<double>::max();
      int shortest_dist_id = -1;

      for(unsigned int j=0; j<predicted.size(); j++){
          double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
          if(distance < shortest_dist){
             shortest_dist = distance;
             shortest_dist_id =predicted[j].id;
             
            }
         observations[i].id = shortest_dist_id;
      }
      
      
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
    
    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
    double gauss_norm =  1.0/(2.0*M_PI*sig_x*sig_y);

    for(unsigned int i=0; i<num_particles; i++) {
      
     
       //double id_part = particles[i].id;
       double x_part = particles[i].x;
       double y_part = particles[i].y;
       double  theta_part = particles[i].theta; 
       
       
       // Map landmarks predicted to be within the particle sensor range
       vector<LandmarkObs> predicted_landmarks;

       for(unsigned int k=0; k<map_landmarks.landmark_list.size(); k++) {
            
            double x_ldm, y_ldm;
            int id_ldm;
            
            id_ldm = map_landmarks.landmark_list[k].id_i;
            x_ldm = map_landmarks.landmark_list[k].x_f;
            y_ldm = map_landmarks.landmark_list[k].y_f;
            
            // Find the landmarks in the map that are within the particle sensor range
            //double dx = x_ldm-x_part;
            //double dy = y_ldm-y_part;
            //double distance =dx*dx+dy*dy;
            
            if(fabs(x_ldm-x_part) <= sensor_range && fabs(y_ldm-y_part) <= sensor_range ){
             // Add landmark to the predicted landmarks
            
             predicted_landmarks.push_back(LandmarkObs{id_ldm, x_ldm, y_ldm});
             }
        }
       
       
       // Transform observed landmarks from vehicle frame to map frame
       vector<LandmarkObs> trans_obs_list;
       for(unsigned int k=0; k<observations.size(); k++) {
            
            double x_obs, y_obs; 
            LandmarkObs trans_obs;

            x_obs = observations[k].x;
            y_obs = observations[k].y;
           
            trans_obs.x = x_part + cos(theta_part)*x_obs - sin(theta_part)*y_obs;
            trans_obs.y = y_part + sin(theta_part)*x_obs + cos(theta_part)*y_obs;
            trans_obs.id = observations[k].id;
            trans_obs_list.push_back(trans_obs);
       }
       
       
       dataAssociation(predicted_landmarks, trans_obs_list);
       particles[i].weight =1.0;
       
       for(unsigned int k=0; k<trans_obs_list.size(); k++){
            
            double  x = trans_obs_list[k].x; 
            double y = trans_obs_list[k].y;
            
            double mu_x;
            double mu_y;
            
            for (unsigned int i=0; i<predicted_landmarks.size(); i++){
                if(trans_obs_list[k].id == predicted_landmarks[i].id){
                 mu_x= predicted_landmarks[i].x;
                 mu_y = predicted_landmarks[i].y;
                }
             } 

           double obs_w = gauss_norm * exp( -( pow(mu_x-x,2)/(2*pow(sig_x, 2)) + (pow(mu_y-y,2)/(2*pow(sig_y, 2))) ) );
           std::cout << "particles[i].weight" << particles[i].weight << std::endl;
           particles[i].weight*=obs_w;
       }
    
   }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    // Get weight and max weight
    vector<double> weights;
    double maxWeight = std::numeric_limits<double>::min();
    for(unsigned int i=0; i<num_particles; i++) {
        weights.push_back(particles[i].weight);
        if(particles[i].weight > maxWeight) {
			maxWeight = particles[i].weight;
		}
    }

    //maxWeight = *max_element(weights.begin(), weights.end());
    std::cout << "maxWeight" << maxWeight << std::endl;
    uniform_real_distribution<double> distDouble(0.0,maxWeight);
    uniform_int_distribution<int> distInt(0,num_particles-1);
    int index = distInt(gen);
    double beta = 0.0;
    vector<Particle> resampledParticles;
    for(unsigned int i=0; i<num_particles; i++) {
        beta+=distDouble(gen)*2.0;
        while(weights[index] < beta){
             beta-=weights[index];
             index=(index+1)%num_particles;
        }
        resampledParticles.push_back(particles[index]);
    } 
   particles = resampledParticles;
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
