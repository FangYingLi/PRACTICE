//
//  stacking.hpp
//  gmx_stacking
//
//  Created by Yiming Tang on 09/04/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#ifndef stacking_hpp
#define stacking_hpp

#include <stdio.h>
#include <gromacs/trajectoryanalysis.h>
#include <string>
#include <vector>

using namespace std;
using namespace gmx;

struct coordinate
{
    float x;
    float y;
    float z;
};

class stacking: public TrajectoryAnalysisModule
{
public:
    stacking();
    
    virtual void initOptions(IOptionsContainer          *options,
                             TrajectoryAnalysisSettings *settings);
    
    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                              const TopologyInformation        &top);
    
    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);
    
    virtual void finishAnalysis(int nframes);
    
    virtual void writeOutput();
    
    //double calDistance(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2);
    //double calAngle(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2);

    float calDistance(coordinate atom11, coordinate atom21);
    float calDistance(coordinate atom11, coordinate atom12, coordinate atom21, coordinate atom22);
    float calDistance(coordinate atom11, coordinate atom12, coordinate atom13,
            coordinate atom21, coordinate atom22, coordinate atom23);
    float calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom14,
                       coordinate atom21, coordinate atom22, coordinate atom23, coordinate atom24);
    float calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom14, coordinate atom15,
                       coordinate atom21, coordinate atom22, coordinate atom23, coordinate atom24, coordinate atom25);
    float calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom14, coordinate atom15, coordinate atom16,
                       coordinate atom21, coordinate atom22, coordinate atom23, coordinate atom24, coordinate atom25, coordinate atom26);



    // If vectorAngle == False, the "benzene" is treated as a plane and the angle is the plane-plane angle.
    float calAngle(coordinate atom11, coordinate atom12, coordinate atom13,
                    coordinate atom21, coordinate atom22, coordinate atom33);

    // If vectorAngle == True, The "main chain" is treated as a vector and the angle is between to vectors;
    float calAngle(coordinate atom11, coordinate atom12,
                    coordinate atom21, coordinate atom22);



private:
    
    AnalysisData    data_probability_;
    
    std::string     fn_density_;
    
    SelectionList   ref_;
    SelectionList   sel_;
	SelectionList   hyd_;
    
    double          cutoff_;
    
    // AnalysisDataAverageModulePointer avem_probability_;

    unsigned long **probability_;
    
    //t_topology     *top_;
    //t_atoms         atoms_;
    
    
    std::string     fnEnergySurface_;
    std::string     fnEnergySurfaceRaw_;
    std::string     fnPropability_;
    
    double          maxDistance_;
    double          stepDistance_;
    double          maxAngle_;
    double          stepAngle_;
    
    double          temperature_;

    // If vector is true: main-chain-like calculation for angles; otherwise benzene-like.
    bool            vector_;
	bool            hbond_;

    bool            not_same_residue_;
    bool            only_same_residue_;

    t_topology     *top_;
    t_atoms         atoms_;
    
};

#endif /* stacking_hpp */
