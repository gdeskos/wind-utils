//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#ifndef MESHMOTIONFIELDS_H
#define MESHMOTIONFIELDS_H

#include "PreProcessingTask.h"

namespace sierra {
namespace nalu {

/**
 * Generate mesh motion field   
 */
class MeshMotionFields: public PreProcessingTask
{
public:
    
    MeshMotionFields(CFDMesh&, const YAML::Node&);

    virtual ~MeshMotionFields() {}

    //! Declare velocity fields and register them for output
    void initialize();

    //! Initialize the velocity fields by linear interpolation
    void run();

private:
    MeshMotionFields() = delete;
    MeshMotionFields(const MeshMotionFields&) = delete;

    //! Parse the YAML file and initialize parameters
    void load(const YAML::Node&);

    //! Parts of the fluid mesh where velocity is initialized
    stk::mesh::PartVector fluid_parts_;

    stk::mesh::PartVector bc_parts_;

    std::vector<double> xs_; 
    std::vector<double> ys_;
    std::vector<double> etas_;
    std::vector<double> phis_;

    double length_; 
    double width_;
    double height_;

    //! X dimension for reading the file
    int Nx_;

    //! Y dimension for reading the file
    int Ny_;

    //! Dimensionality of the mesh
    int ndim_;
 
    //! Exponent for vertical dampening
    double n_exp_;
     
    //! Number of time steps
    size_t numSteps_{0};

    //! Input waves_filename_
    std::string waves_filename_;
    
    std::string output_db_{"waves_history.exo"};

};

}
}

#endif /* MESHMOTIONFIELDS_H */
