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

#include "preprocessing/HOSWaves.h"

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

REGISTER_DERIVED_CLASS(PreProcessingTask, HOSWaves, "write_waves_fields");

HOSWaves::HOSWaves(
    CFDMesh& mesh,
    const YAML::Node& node
) : PreProcessingTask(mesh)
{
    load(node);
}

void HOSWaves::load(const YAML::Node& node)
{
    auto partnames = node["fluid_parts"].as<std::vector<std::string>>();
    partVec_.resize(partnames.size());

    auto& meta = mesh_.meta();
    for(size_t i=0; i < partnames.size(); i++) {
        auto* part = meta.get_part(partnames[i]);
        if (nullptr == part) {
            throw std::runtime_error("Missing fluid part in mesh database: " +
                                     partnames[i]);
        } else {
            partVec_[i] = part;
        }
    }

    waves_filename_ = node["waves_file"].as<std::string>();
    output_db_ = node["time_history_db"].as<std::string>();
    numSteps_ = node["num_timesteps"].as<int>();
}

void HOSWaves::initialize()
{
    const std::string timerName = "HOSWaves::initialize";
    auto timeMon = get_stopwatch(timerName);

    auto& meta = mesh_.meta();
    VectorFieldType& mesh_displacement = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");
    VectorFieldType& mesh_velocity = meta.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_velocity");
    for (auto part: partVec_) {
        stk::mesh::put_field_on_mesh(
            mesh_displacement, *part, meta.spatial_dimension(), nullptr);
        stk::mesh::put_field_on_mesh(
            mesh_velocity, *part, meta.spatial_dimension(), nullptr);
    }
}

void HOSWaves::run()
{
    const std::string timerName = "HOSWaves::run";
    auto timeMon = get_stopwatch(timerName);

    auto& pinfo = get_mpi();
    auto& meta = mesh_.meta();
    auto& bulk = mesh_.bulk();

    const stk::mesh::Selector sel = stk::mesh::selectUnion(partVec_);
    auto& bkts = bulk.get_buckets(stk::topology::NODE_RANK, sel);


    auto* mesh_displacement = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");
    auto* mesh_velocity = meta.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_velocity");
    
    pinfo.info() << "Writing time-history file = " << output_db_ << std::endl;
    std::set<std::string> outfields{"mesh_displacement","mesh_velocity"};
    double time, x,y,eta,phiS;
    mesh_.write_timesteps(
        output_db_, numSteps_, outfields,
        [&](int itime) { 
                std::string file_tmp=waves_filename_+"_"+std::to_string(itime+1)+".dat";
                std::cerr<<file_tmp<<std::endl;
                std::ifstream waves(file_tmp, std::ios::in);
                if (!waves.is_open())
                    throw std::runtime_error(
                    "HOSWaves:: Error opening file: " + waves_filename_);
                waves >> time; 
                std::cerr<<time<<std::endl;
                for (auto b: bkts) {
                    for (auto node: *b) {
                    double* disp = stk::mesh::field_data(*mesh_displacement, node);
                    double* vel = stk::mesh::field_data(*mesh_velocity, node);
                    waves >> x >> y >> eta >> phiS;
                    disp[0] = 0.;
                    disp[1] = 0.;
                    disp[2] = eta;
                    vel[0]=0.;
                    vel[1]=0.;
                    vel[2]=0.;
                    }
                }
            return time;
        });
    pinfo.info() << numSteps_ << " timesteps written successfully" << std::endl;
}

}  // nalu
}  // sierra
