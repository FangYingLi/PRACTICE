//
//  stacking.cpp
//  gmx_stacking
//
//  Created by Yiming Tang on 09/04/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#include "stacking.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <gromacs/trajectoryanalysis.h>
#include <gromacs/fileio/matio.h>
#include <gromacs/fileio/gmxfio.h>



double stacking::calDistance(coordinate atom11, coordinate atom21)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = atom11.x;
    y_1 = atom11.y;
    z_1 = atom11.z;

    x_2 = atom21.x;
    y_2 = atom21.y;
    z_2 = atom21.y;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));

}

double stacking::calDistance(coordinate atom11, coordinate atom12, coordinate atom21, coordinate atom22)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x) / 2;
    y_1 = (atom11.y + atom12.y) / 2;
    z_1 = (atom11.z + atom12.z) / 2;

    x_2 = (atom21.x + atom22.x) / 2;
    y_2 = (atom21.y + atom22.y) / 2;
    z_2 = (atom21.z + atom22.z) / 2;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));

}

double stacking::calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom21,
                             coordinate atom22, coordinate atom23)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x + atom13.x) / 3;
    y_1 = (atom11.y + atom12.y + atom13.y) / 3;
    z_1 = (atom11.z + atom12.z + atom13.z) / 3;

    x_2 = (atom21.x + atom22.x + atom23.x) / 3;
    y_2 = (atom21.y + atom22.y + atom23.y) / 3;
    z_2 = (atom21.z + atom22.z + atom23.z) / 3;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
}

double stacking::calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom14,
                             coordinate atom21, coordinate atom22, coordinate atom23, coordinate atom24)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x + atom13.x + atom14.x) / 4;
    y_1 = (atom11.y + atom12.y + atom13.y + atom14.y) / 4;
    z_1 = (atom11.z + atom12.z + atom13.z + atom14.z) / 4;

    x_2 = (atom21.x + atom22.x + atom23.x + atom24.x) / 4;
    y_2 = (atom21.y + atom22.y + atom23.y + atom24.y) / 4;
    z_2 = (atom21.z + atom22.z + atom23.z + atom24.z) / 4;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
}

double stacking::calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom14,
                             coordinate atom15, coordinate atom21, coordinate atom22, coordinate atom23,
                             coordinate atom24, coordinate atom25)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x + atom13.x + atom14.x + atom15.x) / 5;
    y_1 = (atom11.y + atom12.y + atom13.y + atom14.y + atom15.y) / 5;
    z_1 = (atom11.z + atom12.z + atom13.z + atom14.z + atom15.z) / 5;

    x_2 = (atom21.x + atom22.x + atom23.x + atom24.x + atom25.x) / 5;
    y_2 = (atom21.y + atom22.y + atom23.y + atom24.y + atom25.y) / 5;
    z_2 = (atom21.z + atom22.z + atom23.z + atom24.z + atom25.z) / 5;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
}

double stacking::calDistance(coordinate atom11, coordinate atom12, coordinate atom13, coordinate atom14,
                             coordinate atom15, coordinate atom16, coordinate atom21, coordinate atom22,
                             coordinate atom23, coordinate atom24, coordinate atom25, coordinate atom26)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    x_1 = (atom11.x + atom12.x + atom13.x + atom14.x + atom15.x + atom16.x) / 6;
    y_1 = (atom11.y + atom12.y + atom13.y + atom14.y + atom15.y + atom16.y) / 6;
    z_1 = (atom11.z + atom12.z + atom13.z + atom14.z + atom15.z + atom16.z) / 6;

    x_2 = (atom21.x + atom22.x + atom23.x + atom24.x + atom25.x + atom26.x) / 6;
    y_2 = (atom21.y + atom22.y + atom23.y + atom24.y + atom25.y + atom26.y) / 6;
    z_2 = (atom21.z + atom22.z + atom23.z + atom24.z + atom25.z + atom26.z) / 6;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
}

/*
double stacking::calDistance(int atomNumber, coordinate atom1, ...)
{
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;

    coordinate *ptr;
    ptr = &atom1;

    for(int i=0; i<atomNumber; i++ )
    {
        x_1 += ptr->x;
        y_1 += ptr->y;
        z_1 += ptr->z;
        ptr++;
    }

    x_1 = x_1 / atomNumber;
    y_1 = y_1 / atomNumber;
    z_1 = z_1 / atomNumber;


    for(int i=0; i<atomNumber; i++ )
    {
        x_2 += ptr->x;
        y_2 += ptr->y;
        z_2 += ptr->z;
        ptr++;
    }

    x_2 = x_2 / atomNumber;
    y_2 = y_2 / atomNumber;
    z_2 = z_2 / atomNumber;

    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));


}
*/
/*
double stacking::calDistance(ConstArrayRef<rvec> Group1, ConstArrayRef<rvec> Group2)
{
    ConstArrayRef<rvec>::iterator iter_group1;
    ConstArrayRef<rvec>::iterator iter_group2;
    
    //cout << Group1.size() << '\t';
    //cout << Group2.size() << endl;;
    
    double x_1 = 0.0, y_1 = 0.0, z_1 = 0.0;
    double x_2 = 0.0, y_2 = 0.0, z_2 = 0.0;
    
    for(iter_group1 = Group1.begin(); iter_group1 != Group1.end(); ++iter_group1)
    {
        x_1 += *iter_group1[0];
        y_1 += *iter_group1[1];
        z_1 += *iter_group1[2];
    }
    
    x_1 = x_1 / Group1.size();
    y_1 = y_1 / Group1.size();
    z_1 = z_1 / Group1.size();
    
    for(iter_group2 = Group2.begin(); iter_group2 != Group2.end(); ++iter_group2)
    {
        x_2 += *iter_group2[0];
        y_2 += *iter_group2[1];
        z_2 += *iter_group2[2];
    }
    
    x_2 = x_2 / Group2.size();
    y_2 = y_2 / Group2.size();
    z_2 = z_2 / Group2.size();
    
    
    //cout << x_1 << '\t' << y_1 << '\t' << z_1 << endl;
    //cout << x_2 << '\t' << y_2 << '\t' << z_2 << endl;
    //cout << sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0)) << endl << endl;
    return sqrt( pow(x_1 - x_2, 2.0) + pow(y_1 - y_2, 2.0) + pow(z_1 - z_2, 2.0));
    
}
*/

double stacking::calAngle(coordinate atom11, coordinate atom12, coordinate atom13,
                          coordinate atom21, coordinate atom22, coordinate atom23) {
    double a1 = ((atom12.y - atom11.y) * (atom13.z - atom11.z) - (atom12.z - atom11.z) * (atom13.y - atom11.y));
    double b1 = ((atom12.z - atom11.z) * (atom13.x - atom11.x) - (atom12.x - atom11.x) * (atom13.z - atom11.z));
    double c1 = ((atom12.x - atom11.x) * (atom13.y - atom11.y) - (atom12.y - atom11.y) * (atom13.x - atom11.x));

    double a2 = ((atom22.y - atom21.y) * (atom23.z - atom21.z) - (atom22.z - atom21.z) * (atom23.y - atom21.y));
    double b2 = ((atom22.z - atom21.z) * (atom23.x - atom21.x) - (atom22.x - atom21.x) * (atom23.z - atom21.z));
    double c2 = ((atom22.x - atom21.x) * (atom23.y - atom21.y) - (atom22.y - atom21.y) * (atom23.x - atom21.x));


    double angle = acos((a1 * a2 + b1 * b2 + c1 * c2) / sqrt(a1 * a1 + b1 * b1 + c1 * c1) /
                        sqrt(a2 * a2 + b2 * b2 + c2 * c2)) / (M_PI) * 180;

    double myAngle = min(abs(angle - 0), abs(180 - angle));

    return myAngle;
}

double stacking::calAngle(coordinate atom11, coordinate atom12,
                coordinate atom21, coordinate atom22) {
    double a1 = atom12.x - atom11.x;
    double b1 = atom12.y - atom11.y;
    double c1 = atom12.z - atom11.z;

    double a2 = atom22.x - atom21.x;
    double b2 = atom22.y - atom21.y;
    double c2 = atom22.z - atom21.z;

    double angle = acos((a1 * a2 + b1 * b2 + c1 * c2) / sqrt(a1 * a1 + b1 * b1 + c1 * c1) /
                        sqrt(a2 * a2 + b2 * b2 + c2 * c2)) / (M_PI) * 180;

    return angle;
}

stacking::stacking(): cutoff_(1.5), maxDistance_(1.2), maxAngle_(90.0), stepDistance_(0.01), stepAngle_(1.0)
{
    registerAnalysisDataset(&data_probability_, "probability");
}

void stacking::initOptions(gmx::IOptionsContainer *options, gmx::TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
    {
        "Stacking Free Energy Surface (SFES) Calculator by Yiming Tang @ Fudan. \n",
        "This Tool calculates the Joint Probabilit Density (JPD) of two molecules",
        "as a function of angle (0-90) - distance (0-1.2). The tool ouputs the free",
        "energy surface defined by E = - R * T * log(H), where H is the probability.",
        "This tool takes two selection groups named ref and sel, each should contains",
        "groups each of which contains one benzene. Angles and Distances will be calculated",
        "on each group pair between ref and sel."
    };
    
    //settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    
    settings->setHelpText(desc);
    
    options->addOption(FileNameOption("energy").filetype(eftUnknown).legacyType(efXPM).outputFile()
                       .store(&fnEnergySurface_).defaultBasename("energy")
                       .description("Energy Surface xpm map").required());
    
    options->addOption(FileNameOption("raw-energy").filetype(eftGenericData).outputFile()
                       .store(&fnEnergySurfaceRaw_)
                       .description("Ras Energy Surface Data for further analysis"));
    
    options->addOption(FileNameOption("raw-probability").filetype(eftGenericData).outputFile()
                       .store(&fnPropability_)
                       .description("Raw Probability Data for further analysis"));
    
    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                       .description("Benzenes whose minimum distance are beyond this cutoff is not considered contacted"));
    
    options->addOption(DoubleOption("maxD").store(&maxDistance_).defaultValue(1.2)
                       .description("Maxmimum Centroid Distance to output."));
    
    options->addOption(DoubleOption("stepD").store(&stepDistance_).defaultValue(0.01)
                       .description("Step of Centroid Distance to output."));
    
    options->addOption(DoubleOption("maxA").store(&maxAngle_).defaultValue(90)
                       .description("Maxmimum Centroid Angle to output, default is 90 for plane and 180 for vector."));
    
    options->addOption(DoubleOption("stepA").store(&stepAngle_).defaultValue(1)
                       .description("Step of Centroid Angle to output."));
    
    options->addOption(SelectionOption("ref")
                       .storeVector(&ref_).required().multiValue()
                       .description("Reference Groups of benzenes to calculate angle and distance"));
    
    options->addOption(SelectionOption("sel")
                       .storeVector(&sel_).required().multiValue()
                       .description("Selection Groups of benzenes to calculate angle and distance"));
    
    options->addOption(DoubleOption("temperature").store(&temperature_).defaultValue(298)
                       .description("Temperature in K for energy calculation"));

    options->addOption(BooleanOption("vector")
                       .store(&vector_).required().defaultValue(false)
                       .description("Whether selections are processed as vectors when calculating angles."));
}

void stacking::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings, const gmx::TopologyInformation &top)
{
    //atoms_ = top.topology()->atoms;
    //top_   = top.topology();
    
 
    if(vector_ && maxAngle_ == 90)
    {
        maxAngle_ = 180;
    }

   // initial data set columns and rows
    
    data_probability_.setDataSetCount((int)(maxDistance_ / stepDistance_));
    
    for (int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        data_probability_.setColumnCount(i, (int)(maxAngle_ / stepAngle_));
    }
    
    avem_probability_.reset(new AnalysisDataAverageModule());
    data_probability_.addModule(avem_probability_);
    
}

void stacking::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, gmx::TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle dhProbability = pdata->dataHandle(data_probability_);
    dhProbability.startFrame(frnr, fr.time);
    
    // Constructed a new temp all-zero matrix
    
    int **tempProbability = new int*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        tempProbability[i] = new int[(int)(maxAngle_ / stepAngle_)];
    }
    for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
        {
            tempProbability[i][j] = 0;
        }
    }
    
    // Begin analysis
    
    for (size_t gRef = 0 ; gRef < ref_.size(); ++gRef)
    {
        for (size_t gSel = 0; gSel < sel_.size(); ++gSel)
        {
            //cout << "ref_size = " << ref_.size() << endl;
            //cout << "sel_size = " << sel_.size() << endl;
            
            
            const Selection &ref = pdata->parallelSelection(ref_[gRef]);
            const Selection &sel = pdata->parallelSelection(sel_[gSel]);
            
            if (ref.atomIndices()[0] == sel.atomIndices()[0])
            {
                continue;
            }
            

            
            ////// Calculating distance //////////////
            


            coordinate *molecule1;
            coordinate *molecule2;
            double tempDistance;


            int atomCount = (int) min(ref.coordinates().size(), sel.coordinates().size());

            if(atomCount <= 6)
            {
                molecule1 = new coordinate[atomCount];
                for (int i=0; i<atomCount; i++)
                {
                    molecule1[i].x = ref.coordinates()[i][0];
                    molecule1[i].y = ref.coordinates()[i][1];
                    molecule1[i].z = ref.coordinates()[i][2];
                }

                molecule2 = new coordinate[atomCount];
                for (int i=0; i<atomCount; i++)
                {
                    molecule2[i].x = sel.coordinates()[i][0];
                    molecule2[i].y = sel.coordinates()[i][1];
                    molecule2[i].z = sel.coordinates()[i][2];
                }

            }
            else
            {
                molecule1 = new coordinate[6];
                for (int i=0; i<3; i++)
                {
                    molecule1[i].x = ref.coordinates()[i][0];
                    molecule1[i].y = ref.coordinates()[i][1];
                    molecule1[i].z = ref.coordinates()[i][2];

                    molecule1[3+i].x = ref.coordinates()[atomCount-4+i][0];
                    molecule1[3+i].y = ref.coordinates()[atomCount-4+i][1];
                    molecule1[3+i].z = ref.coordinates()[atomCount-4+i][2];
                }

                molecule2 = new coordinate[6];
                for (int i=0; i<3; i++)
                {
                    molecule2[i].x = sel.coordinates()[i][0];
                    molecule2[i].y = sel.coordinates()[i][1];
                    molecule2[i].z = sel.coordinates()[i][2];

                    molecule2[3+i].x = sel.coordinates()[atomCount-4+i][0];
                    molecule2[3+i].y = sel.coordinates()[atomCount-4+i][1];
                    molecule2[3+i].z = sel.coordinates()[atomCount-4+i][2];
                }
                atomCount = 6;
            }


            switch(atomCount) {
                case 1:
                    tempDistance = calDistance(molecule1[0], molecule2[0]);
                    break;
                case 2:
                    tempDistance = calDistance(molecule1[0], molecule1[1],
                                               molecule2[0], molecule2[1]);
                    break;
                case 3:
                    tempDistance = calDistance(molecule1[0], molecule1[1], molecule1[2],
                                               molecule2[0], molecule2[1], molecule2[2]);
                    break;
                case 4:
                    tempDistance = calDistance(molecule1[0], molecule1[1], molecule1[2], molecule1[3],
                                               molecule2[0], molecule2[1], molecule2[2], molecule2[3]);
                    break;
                case 5:
                    tempDistance = calDistance(molecule1[0], molecule1[1], molecule1[2],
                                               molecule1[3], molecule1[4],
                                               molecule2[0], molecule2[1], molecule2[2],
                                               molecule2[3], molecule2[4]);
                    break;
                case 6:
                    tempDistance = calDistance(molecule1[0], molecule1[1], molecule1[2],
                                               molecule1[3], molecule1[4], molecule1[5],
                                               molecule2[0], molecule2[1], molecule2[2],
                                               molecule2[3], molecule2[4], molecule2[5]);
                    break;
                default:
                    cerr << "Wrong Number of atomCount" << endl;
                    exit(1);


            }

            /// Calculating Angle ////////////////////

            atomCount = (int) min(ref.coordinates().size(), sel.coordinates().size());

            double tempAngle;

            if(vector_) {
                if (atomCount == 2 || atomCount == 3) {
                    coordinate p11 = {ref.coordinates()[0][0], ref.coordinates()[0][1], ref.coordinates()[0][2]};
                    coordinate p12 = {ref.coordinates()[atomCount - 1][0], ref.coordinates()[atomCount - 1][1],
                                      ref.coordinates()[atomCount - 1][2]};

                    coordinate p21 = {sel.coordinates()[0][0], sel.coordinates()[0][1], sel.coordinates()[0][2]};
                    coordinate p22 = {sel.coordinates()[atomCount - 1][0], sel.coordinates()[atomCount - 1][1],
                                      sel.coordinates()[atomCount - 1][2]};

                    tempAngle = calAngle(p11, p12, p21, p22);

                } else {
                    if (atomCount > 3) {
                        coordinate p11 = {(ref.coordinates()[0][0] + ref.coordinates()[1][0]) / 2,
                                          (ref.coordinates()[0][1] + ref.coordinates()[1][1]) / 2,
                                          (ref.coordinates()[0][2] + ref.coordinates()[1][2]) / 2,};

                        coordinate p12 = {
                                (ref.coordinates()[atomCount - 1][0] + ref.coordinates()[atomCount - 2][0]) / 2,
                                (ref.coordinates()[atomCount - 1][1] + ref.coordinates()[atomCount - 2][1]) / 2,
                                (ref.coordinates()[atomCount - 1][2] + ref.coordinates()[atomCount - 2][2]) / 2,};

                        coordinate p21 = {(sel.coordinates()[0][0] + sel.coordinates()[1][0]) / 2,
                                          (sel.coordinates()[0][1] + sel.coordinates()[1][1]) / 2,
                                          (sel.coordinates()[0][2] + sel.coordinates()[1][2]) / 2,};

                        coordinate p22 = {
                                (sel.coordinates()[atomCount - 1][0] + sel.coordinates()[atomCount - 2][0]) / 2,
                                (sel.coordinates()[atomCount - 1][1] + sel.coordinates()[atomCount - 2][1]) / 2,
                                (sel.coordinates()[atomCount - 1][2] + sel.coordinates()[atomCount - 2][2]) / 2,};

                        tempAngle = calAngle(p11, p12, p21, p22);
                    } else {
                        cerr << "You have groups that contains less than three atoms!" << endl;
                        exit(1);
                    }
                }
            } else {

                coordinate p11 = {ref.coordinates()[0][0], ref.coordinates()[0][1], ref.coordinates()[0][2]};
                coordinate p12 = {ref.coordinates()[1][0], ref.coordinates()[1][1], ref.coordinates()[1][2]};
                coordinate p13 = {ref.coordinates()[2][0], ref.coordinates()[2][1], ref.coordinates()[2][2]};

                coordinate p21 = {sel.coordinates()[0][0], sel.coordinates()[0][1], sel.coordinates()[0][2]};
                coordinate p22 = {sel.coordinates()[1][0], sel.coordinates()[1][1], sel.coordinates()[1][2]};
                coordinate p23 = {sel.coordinates()[2][0], sel.coordinates()[2][1], sel.coordinates()[2][2]};

                tempAngle = calAngle(p11, p12, p13, p21, p22, p23);

            }


            
            if (tempDistance <= maxDistance_ && tempAngle <= maxAngle_ )
            {
                tempProbability[(int)floor(tempDistance / stepDistance_)][(int)floor(tempAngle/stepAngle_)] += 1;
            }
        }
    }
    
    for(int i = 0; i < maxDistance_ / stepDistance_; i++)
    {
        for(int j = 0; j < maxAngle_ / stepAngle_; j++)
        {
            dhProbability.selectDataSet(i);
            dhProbability.setPoint(j, tempProbability[i][j]);
            if (tempProbability[i][j] != 0)
            {
            //cout << tempProbability[i][j];
            }
        }
    }
    
    dhProbability.finishFrame();
    
    
}

void stacking::finishAnalysis(int /*nframes*/)
{
    
}

void stacking::writeOutput()
{
    
    real DistanceVector[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_) ; i++)
    {
        DistanceVector[i] = real(i * stepDistance_);
    }
    
    real AngleVector[(int)(maxAngle_ / stepAngle_)];
    
    for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
    {
        AngleVector[j] = real(j * stepAngle_);
    }
    
    
    if       (fnEnergySurface_.empty())         {fnEnergySurface_ = "energy.xpm";}
    else if  (fnEnergySurface_.compare(".xpm")) {}
    else                                        {fnEnergySurface_ += ".xpm";}
    
    // Construct matrix probability
    real **matProbability = new real*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matProbability[i] = new real[(int)(maxAngle_ / stepAngle_)];
    }
    
    real **matEnergy = new real*[(int)(maxDistance_ / stepDistance_)];
    for(int i = 0; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        matEnergy[i] = new real[(int)(maxAngle_ / stepAngle_)];
    }
    
    int scale = (int)(sel_.size()) * (int)(ref_.size());
    
    for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
    {
        for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
        {
            matProbability[i][j] = avem_probability_->average(i, j) ;
            //cout << i << '\t' << j << '\t' << avem_probability_->average(i, j) << '\t' << avem_probability_->average(i, j) / (float)scale << '\t' << matProbability[i][j] << endl;
            matEnergy[i][j]      = real(- 8.31447 * temperature_ * log(matProbability[i][j] / (float)scale) * 0.0002389);
        }
    }
    
    t_rgb rlo, rhi;
    rlo.r = 0.0; rlo.g = 0.0; rlo.b = 1.0;
    rhi.r = 1.0; rhi.g = 0.0; rhi.b = 0.0;
    int nlevels  = 400;
    
    FILE *fpEnergySurface;
    fpEnergySurface = fopen(fnEnergySurface_.c_str(), "w");
    
    
    
    write_xpm(fpEnergySurface, 0, "Free Energy Surface"
              , "Contact Probability", "Distance", "angle"
              , (int)(maxDistance_ / stepDistance_), (int)(maxAngle_ / stepAngle_), DistanceVector, AngleVector
              , matEnergy, 0, 20, rlo, rhi, &nlevels);
    
    fclose(fpEnergySurface);
    
    if (!fnEnergySurfaceRaw_.empty())
    {
        if (fnEnergySurfaceRaw_.compare(".dat"))    {}
        else                                        {fnEnergySurfaceRaw_ += ".dat";}
        
        FILE *fpEnergySurfaceRaw;
        fpEnergySurfaceRaw = fopen(fnEnergySurfaceRaw_.c_str(), "w" );
        
        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
            {
                fprintf(fpEnergySurfaceRaw, "%f ", matEnergy[i][j]);
            }
            fprintf(fpEnergySurfaceRaw, "\n");
        }
        fclose(fpEnergySurfaceRaw);
        
        
    }
    
    if (!fnPropability_.empty())
    {
        if (fnPropability_.compare(".dat"))    {}
        else                                        {fnPropability_ += ".dat";}
        
        FILE *fpProbability;
        fpProbability = fopen(fnPropability_.c_str(), "w" );
        
        for(int i = 0 ; i < (int)(maxDistance_ / stepDistance_); i++)
        {
            for(int j = 0; j < (int)(maxAngle_ / stepAngle_); j++)
            {
                fprintf(fpProbability, "%f ", matProbability[i][j]);
            }
            fprintf(fpProbability, "\n");
        }
        fclose(fpProbability);
        
    }
    
}


















