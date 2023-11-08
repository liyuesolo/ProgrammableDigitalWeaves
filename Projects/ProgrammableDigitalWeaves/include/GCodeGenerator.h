#ifndef GCODE_GENERATOR_H
#define GCODE_GENERATOR_H

#include <fstream>

#include "EoLRodSim.h"
template<class T, int dim>
class EoLRodSim;

enum PrinterType
{
    PrusaI3, UltiMaker
};

enum ExtrusionMode
{
    Absolute, Relative
};

template<class T, int dim>
class GCodeGenerator
{
public:
    using TV = Vector<T, dim>;
    using TV3 = Vector<T, 3>;
    using TV2 = Vector<T, 2>;
    using Range = Vector<T, 2>;
    using Offset = Vector<int, dim + 1>;

// generic printing stuff
private:
    const EoLRodSim<T, dim>& sim;
    std::string gcode_file;
	T nozzle_diameter;
    T filament_diameter;
    T extruder_temperature;
    T bed_temperature;
    T first_layer_height;
    T layer_height;
    T feed_rate_move;
    T feed_rate_print;
    PrinterType printer;
	ExtrusionMode extrusion_mode;

    T current_E;
    TV current_position = TV::Zero();

    std::ofstream gcode;

// specific parameters
private:
    T tunnel_height = 0.2; //0.2mm
public:
    GCodeGenerator(const EoLRodSim<T, dim>& _sim) : sim(_sim), gcode_file("./rod.gcode") {}
    GCodeGenerator(const EoLRodSim<T, dim>& _sim, const std::string& filename);
    ~GCodeGenerator() {}

public:

    void buildGridScene(int n_row, int n_col, int type);

    void generateShelterScene();

    void generateGCodeShelter();
    void generateGCodeSingleStrand();

    void generateGCodeClosedGrid(int n_row, int n_col, int type);
    
    void circlePatchGCode(int n_row, int n_col, int type, bool add_bar);

    void slidingBlocksGCode(int n_row, int n_col, int type, bool add_bar);

    void activeTexticleGCode(bool fused = false);
    void activeTexticleGCode2(bool fused = false);

    void crossingTest();

    void generateGCodeFromRodsCurveGripperHardCoded();

    void generateGCodeFromRodsGridGripperHardCoded();
    void generateGCodeFromRodsFixedGridGripperHardCoded();

    void generateGCodeFromRodsShelterHardCoded();
    void generateGCodeFromRodsGridHardCoded(int n_row, int n_col, int type);

    void generateGCodeFromRodsNoTunnel();

    void writeCircle(const TV& center, T r, const TV& start, const Vector<T, 2>& range,
         const std::vector<TV>& lifting_points, int sub_div, T rod_diameter);

    void addRecBar(const TV& border_a, const TV& border_b, T width, T rod_diameter);

    void writeLine(const TV& from, const TV& to, T rod_radius, T speed = 800.0);
    void moveTo(const TV& to, T speed = 2000.0, bool do_retract = true);

    void addSingleTunnel(const TV& from, const TV& to, T height, T rod_diameter = 0.2);
    void addSingleTunnel(const TV& left, const TV& center, const TV& right, T rod_diameter);

    void addSingleTunnelOnCrossing(int crossing_id, const TV3& heights, 
        int direction, std::function<void(TV&)> scaleAndShift);
    
    void addSingleTunnelOnCrossingWithFixedRange(int crossing_id, const TV3& heights, 
        int direction, std::function<void(TV&)> scaleAndShift,
        const Range& range, T extend_right = 0.3, T speed_first_half = 100, T speed_second_half = 300);

    void addTunnelsAlongRod(int rod_idx, const TV3& heights,
        std::function<void(TV&)> scaleAndShift, int cross_n_seg,
        T extend_right = 0.3, T speed_first_half = 100, T speed_second_half = 300);

    void generateCodeSingleRod(int rod_idx, std::function<void(TV&)> scaleAndShift, 
        bool is_first_layer,
        T bd_height = 0.3, T inner_height = 0.3, T buffer_percentage = 0.3, T less = false, T extend = false);
    
    void generateCodeSingleRodMoveUpNozzle(int rod_idx, std::function<void(TV&)> scaleAndShift, 
        bool is_first_layer,
        T bd_height = 0.3, T inner_height = 0.3, T buffer_percentage = 0.3, T less = false, T extend = false);
private:

    T computeExtrusionAmount() const;
    T crossSectionAreaFromRod(T rod_radius) const;

    T extrusionWidth() const;
    T crossSectionArea(bool is_first_layer) const;
    void retract(T E);
    void extrude(T E);
    void writeHeader();
    void writeFooter();

    
};


#endif