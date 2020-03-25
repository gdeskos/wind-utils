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


#include "MeshMotionFields.h"
#include "core/YamlUtils.h"
#include "core/KokkosWrappers.h"
#include "core/PerfUtils.h"
#include "core/ParallelInfo.h"

#include "stk_mesh/base/TopologyDimensions.hpp"
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_mesh/base/Field.hpp"

#include <fstream>
#include <vector>
#include <string>

namespace sierra {
namespace nalu {

REGISTER_DERIVED_CLASS(PreProcessingTask, MeshMotionFields, "create_meshmotion_fields");

MeshMotionFields::MeshMotionFields(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void MeshMotionFields::load(const YAML::Node& node)
{

    auto fluid_partnames = node["fluid_parts"].as<std::vector<std::string>>();
    fluid_parts_.resize(fluid_partnames.size());
    

    auto& meta = mesh_.meta();
    for(size_t i=0; i<fluid_partnames.size(); i++) {
        auto* part = meta.get_part(fluid_partnames[i]);
        if (nullptr == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     fluid_partnames[i]);
        } else {
            fluid_parts_[i] = part;
        }
    } 
    
    Nx_ = node["Nx"].as<int>();
    Ny_ = node["Ny"].as<int>();
    n_exp_=node["damping_exponent"].as<int>(); 
    waves_filename_ = node["waves_file"].as<std::string>();
    output_db_      = node["time_history_db"].as<std::string>();
    numSteps_       = node["num_timesteps"].as<int>();
}

void MeshMotionFields::initialize()
{

    const std::string timerName = "HOSWaves::initialize";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    
    VectorFieldType& mesh_displacement = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");
    VectorFieldType& mesh_velocity = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_velocity");
    for (auto part: fluid_parts_) {
        stk::mesh::put_field_on_mesh(
            mesh_displacement, *part, meta.spatial_dimension(), nullptr);
        stk::mesh::put_field_on_mesh(
            mesh_velocity, *part, meta.spatial_dimension(), nullptr);
    }
    // Initialize the vectors that read into the .dat files
    xs_.resize(Nx_,Ny_); ys_.resize(Nx_,Ny_);
    etas_.resize(Nx_,Ny_); phis_.resize(Nx_,Ny_);

}

void MeshMotionFields::run()
{
    const std::string timerName = "MeshMotion::run";
    auto timeMon = get_stopwatch(timerName);

    auto& pinfo = get_mpi();
    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();
    
    const stk::mesh::Selector sel = stk::mesh::selectUnion(fluid_parts_);
    auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);
 
    // Computing the box dimensions
    auto bbox = mesh_.calc_bounding_box(sel,false);
    length_ = bbox.get_x_max() - bbox.get_x_min();
    width_ = bbox.get_y_max() - bbox.get_y_min();
    height_ = bbox.get_z_max() - bbox.get_z_min();
    
    VectorFieldType* coords = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* mesh_disp = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_diplacement");
    VectorFieldType* mesh_vel = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_velocity");
    
    std::set<std::string> outfields{"mesh_displacement","mesh_velocity"};
    double time, xs,ys,etas, phis;
    mesh_.write_timesteps(
        output_db_, numSteps_, outfields,
        [&](int itime) { 
                std::string file_tmp=waves_filename_+"_"+std::to_string(itime)+".dat";
                std::cerr<<file_tmp<<std::endl;
                std::ifstream waves(file_tmp, std::ios::in);
                if (!waves.is_open())
                    throw std::runtime_error(
                    "HOSWaves:: Error opening file: " + waves_filename_);
                waves >> time; 
                std::cerr<<time<<std::endl;
                for (int i=0; i<Nx_;i++){
                    for (int j=0; j<Ny_; j++){
                        waves>>xs>>ys>>etas>>phis;
                        //std::cerr<<xs<<ys<<etas<<phis<<std::endl;
                        //xs_[i][j]=xs;
                        //ys_[i][j]=ys;
                        //etas_[i][j]=etas;
                        //phis_[i][j]=phisz;
                    }
                }
                
                for (auto b: bkts) {
                    for (auto node: *b) {
                        double* xyz = stk::mesh::field_data(*coords,node);
                        double* disp = stk::mesh::field_data(*mesh_disp, node);
                        double* vel = stk::mesh::field_data(*mesh_vel, node); 
                        disp[0]=-time;
                        disp[1]=-2*time;
                        disp[2]=-3*time;
                        //vel[0]=time;
                        //vel[1]=2*time;
                        //vel[2]=3*time;
                    }
                }
    //                const int i = std::floor(x / dx);
    //                const int j = std::floor(y / dy);
            return time;
        });
    pinfo.info() << numSteps_ << " timesteps written successfully" << std::endl;
}

} // nalu
} // sierra
