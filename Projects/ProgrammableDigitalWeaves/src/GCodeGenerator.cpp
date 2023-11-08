#include "../include/GCodeGenerator.h"

template<class T, int dim>
GCodeGenerator<T, dim>::GCodeGenerator(const EoLRodSim<T, dim>& _sim, 
    const std::string& filename) 
    : sim(_sim), 
    gcode_file(filename),
    nozzle_diameter(0.4),
    filament_diameter(1.75),
    extruder_temperature(180),
    bed_temperature(60),
    first_layer_height(0.2), 
    layer_height(0.2),
    feed_rate_move(5200),
	feed_rate_print(300),
    printer(PrusaI3),
    extrusion_mode(Absolute),
    current_E(0.0)
{

}

template<class T, int dim>
void GCodeGenerator<T, dim>::buildGridScene(int n_row, int n_col, int type)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 40);
    };

    T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
    
    if constexpr (dim == 3)
    {
        writeHeader();
        T extend_percentage = 0.3;
        TV bottom_left, top_right;
        sim.computeBoundingBox(bottom_left, top_right);

        scaleAndShift(bottom_left); scaleAndShift(top_right);
        bottom_left[dim-1] = rod_radius_in_mm;
        top_right[dim-1] = rod_radius_in_mm;

        TV bottom_right = bottom_left;
        bottom_right[0] = top_right[0];
        TV top_left = top_right;
        top_left[0] = bottom_left[0];

        TV bottom_left_extend = bottom_left - (bottom_right - bottom_left) * 0.2;

        if (type == 0)
        {
            for (int row = 0; row < n_row; row++)
            {
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = 0.5 * rod_radius_in_mm; right[dim-1] = 0.5 * rod_radius_in_mm;
                from[dim-1] = 0.5 * rod_radius_in_mm; to[dim-1] = 0.5 * rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, 0.5 * rod_radius_in_mm, 600);
                    writeLine(from, to, 0.5 * rod_radius_in_mm, 1200);
                    writeLine(to, right, 0.5 * rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, 0.5 * rod_radius_in_mm, 600);
                    writeLine(to, from, 0.5 * rod_radius_in_mm, 1200);
                    writeLine(from, left, 0.5 * rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            for (int col = 0; col < n_col; col++)
            {
                // generateCodeSingleRod(col + n_row, scaleAndShift, true, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[col + n_row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                scaleAndShift(from); scaleAndShift(to);
                from[dim-1] = 3.0 * rod_radius_in_mm; to[dim-1] = 3.0 * rod_radius_in_mm;
                
                if (col % 2 == 0)
                {
                    
                    moveTo(right);
                    right[dim-1] = 0.2;
                    moveTo(right);
                    to[1] = top_right[1];
                    TV temp = 0.5 * (right + to);
                    temp[2] = 0.2;
                    writeLine(right, temp, rod_radius_in_mm, 200);
                    writeLine(temp, to, rod_radius_in_mm, 400);
                    writeLine(to, from, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (from + left);
                    temp[2] = 0.2;
                    writeLine(from, temp, rod_radius_in_mm, 200);
                    writeLine(temp, left, rod_radius_in_mm, 200);
                }
                else
                {
                    moveTo(left);
                    left[dim-1] = 0.2;
                    moveTo(left);
                    from[1] = bottom_left[1];
                    TV temp = 0.5 * (left + from);
                    temp[2] = 0.2;
                    writeLine(left, temp, rod_radius_in_mm, 200);
                    writeLine(temp, from, rod_radius_in_mm, 400);
                    writeLine(from, to, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (to + right);
                    temp[2] = 0.2;
                    writeLine(to, temp, rod_radius_in_mm, 200);
                    writeLine(temp, right, rod_radius_in_mm, 200);
                }
            }
            feed_rate_print = 800;
            for (int row = 0; row < n_row; row++)
                generateCodeSingleRodMoveUpNozzle(row, scaleAndShift, false, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0, false, true);
            
            TV top_right_extend = top_right + (top_right - top_left) * 0.2;

            moveTo(top_right_extend);
            top_left[dim - 1] = 1.5 * rod_radius_in_mm;
            bottom_left[dim - 1] = 1.5 * rod_radius_in_mm;
            bottom_right[dim - 1] = 1.5 * rod_radius_in_mm;
            top_right[dim - 1] = 1.5 * rod_radius_in_mm;
            writeLine(top_right_extend, top_left, rod_radius_in_mm, 1500);
            writeLine(top_left, bottom_left, rod_radius_in_mm, 1500);
            writeLine(bottom_left, bottom_right, rod_radius_in_mm, 1500);
            writeLine(bottom_right, top_right, rod_radius_in_mm, 1500);
        }
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeClosedGrid(int n_row, int n_col, int type)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 40);
    };

    T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
    
    if constexpr (dim == 3)
    {
        writeHeader();
        T extend_percentage = 0.2;
        TV bottom_left, top_right;
        sim.computeBoundingBox(bottom_left, top_right);

        scaleAndShift(bottom_left); scaleAndShift(top_right);
        bottom_left[dim-1] = rod_radius_in_mm;
        top_right[dim-1] = rod_radius_in_mm;

        TV bottom_right = bottom_left;
        bottom_right[0] = top_right[0];
        TV top_left = top_right;
        top_left[0] = bottom_left[0];

        TV bottom_left_extend = bottom_left - (bottom_right - bottom_left) * 0.2;

        
        
        if (type == 0)
        {
            for (int row = 1; row < n_row - 1; row++)
            {
                
                // generateCodeSingleRod(row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            for (int col = 0; col < n_col; col++)
            {
                // generateCodeSingleRod(col + n_row, scaleAndShift, true, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[col + n_row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                T rod_height = (col == 0 || col == n_col - 1) ? rod_radius_in_mm : 3.0 * rod_radius_in_mm;
                scaleAndShift(left); scaleAndShift(right);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                scaleAndShift(from); scaleAndShift(to);
                from[dim-1] = rod_height; to[dim-1] = rod_height;
                
                if (col % 2 == 0)
                {
                    
                    moveTo(right);
                    right[dim-1] = 0.2;
                    moveTo(right);
                    to[1] = top_right[1];
                    TV temp = 0.5 * (right + to);
                    temp[2] = 0.2;

                    writeLine(right, temp, rod_radius_in_mm, 200);
                    writeLine(temp, to, rod_radius_in_mm, 400);
                    writeLine(to, from, rod_height+ rod_radius_in_mm, 600);
                    temp = 0.5 * (from + left);
                    temp[2] = 0.2;
                    writeLine(from, temp, rod_radius_in_mm, 200);
                    writeLine(temp, left, rod_radius_in_mm, 200);
                }
                else
                {
                    moveTo(left);
                    left[dim-1] = 0.2;
                    moveTo(left);
                    from[1] = bottom_left[1];
                    TV temp = 0.5 * (left + from);
                    temp[2] = 0.2;
                    writeLine(left, temp, rod_radius_in_mm, 200);
                    writeLine(temp, from, rod_radius_in_mm, 400);
                    writeLine(from, to, rod_height + rod_radius_in_mm, 600);
                    temp = 0.5 * (to + right);
                    temp[2] = 0.2;
                    writeLine(to, temp, rod_radius_in_mm, 200);
                    writeLine(temp, right, rod_radius_in_mm, 200);
                }
            }

            for (int row : {0, n_row - 1})
            {
                // generateCodeSingleRod(row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            // T width = 0.034;
            T width = 0.015;
            tunnel_height = 0.2;
            for (int row = 1; row < n_row - 1; row++)
            {
                for (int col = 1; col < n_col - 1; col++)
                {
                    int crossing_id = row * n_col + col;
                    addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(width, width), 0.2, 100, 200);
                }
            }
            
            // moveTo(bottom_left_extend);
            // writeLine(bottom_left_extend, bottom_right, rod_radius_in_mm, 1500);
            // writeLine(bottom_right, top_right, rod_radius_in_mm, 1500);
            // writeLine(top_right, top_left, rod_radius_in_mm, 1500);
            // writeLine(top_left, bottom_left, rod_radius_in_mm, 1500);
        }
    }
    writeFooter();
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeSingleStrand()
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(60, 60);
    };

    if constexpr (dim == 3)
    {
        
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
        writeHeader();
        
        // sim.Rods[2]->fixed_by_crossing[0] = false;
        // sim.Rods[2]->fixed_by_crossing[2] = false;
        // for (int i : {0, 2, 4, 6})
        // {
        //     sim.Rods[0]->fixed_by_crossing[i] = false;
        //     sim.Rods[1]->fixed_by_crossing[i] = false;
        // }
        // for (int i = 0; i < sim.Rods.size(); i++)
        feed_rate_print = 400;
        for (int i : {0, 1, 2})
        {
            sim.Rods[i]->indices.push_back(sim.Rods[i]->indices[0]);
            sim.Rods[i]->indices.push_back(sim.Rods[i]->indices[1]);
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.0);   
            sim.Rods[i]->indices.pop_back();
            sim.Rods[i]->indices.pop_back();
        }
        for (int i : {5, 6, 7, 8})
        {
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.0);   
        }

        generateCodeSingleRod(4, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.0);   
        for (int i = 0; i < sim.Rods[3]->dof_node_location.size(); i++)
        {
            if (i>0)
                sim.Rods[3]->fixed_by_crossing[i] = false;
        }
        generateCodeSingleRod(3, scaleAndShift, true, rod_radius_in_mm, 4.0 * rod_radius_in_mm, 0.2, false, true);   

        TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
        T width = 0.5;

        addTunnelsAlongRod(0, heights, scaleAndShift, 4);
        addTunnelsAlongRod(1, heights, scaleAndShift, 2);
        addTunnelsAlongRod(2, heights, scaleAndShift, 2);

        addSingleTunnelOnCrossingWithFixedRange(24, heights, 1, scaleAndShift, Range(0.03, 0.05), 0.3, 150, 200);
        // for (int i : {0, 1, 2, 3, 4, 5} )
        //     addTunnelsAlongRod(i, heights, scaleAndShift, 10);

        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeShelter()
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(60, 60);
    };

    if constexpr (dim == 3)
    {
        
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
        writeHeader();
        
        // sim.Rods[2]->fixed_by_crossing[0] = false;
        // sim.Rods[2]->fixed_by_crossing[2] = false;
        // for (int i : {0, 2, 4, 6})
        // {
        //     sim.Rods[0]->fixed_by_crossing[i] = false;
        //     sim.Rods[1]->fixed_by_crossing[i] = false;
        // }
        // for (int i = 0; i < sim.Rods.size(); i++)
        feed_rate_print = 900;
        for (int i : {0, 1, 2, 3, 4, 5})
        {
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.0);   
        }
        for (int i : {6, 7, 8, 9, 10, 11} )
        {
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, 4.0 * rod_radius_in_mm, 0.0);   
        }

        TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
        T width = 0.5;
        for (int i : {0, 1, 2, 3, 4, 5} )
            addTunnelsAlongRod(i, heights, scaleAndShift, 10);

        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateShelterScene()
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(100, 100);
    };

    if constexpr (dim == 3)
    {
        
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
        writeHeader();
        feed_rate_print = 600;
        // for (int i = 0; i < sim.Rods.size(); i++)
        sim.Rods[2]->fixed_by_crossing[0] = false;
        sim.Rods[2]->fixed_by_crossing[2] = false;
        for (int i : {0, 2, 4, 6})
        {
            sim.Rods[0]->fixed_by_crossing[i] = false;
            sim.Rods[1]->fixed_by_crossing[i] = false;
        }
        for (int i : {0, 1, 2})
        {
            sim.Rods[i]->indices.push_back(sim.Rods[i]->indices[0]);
            sim.Rods[i]->indices.push_back(sim.Rods[i]->indices[1]);
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.0);   
            sim.Rods[i]->indices.pop_back();
            sim.Rods[i]->indices.pop_back();
        }

        for (int i : {5, 6, 7, 8})
        {
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.0);   
        }

        for (int i : {3, 4} )
        {
            generateCodeSingleRod(i, scaleAndShift, true, rod_radius_in_mm, 4.0 * rod_radius_in_mm, 0.2, false, true);   
        }

        TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 5.0 * rod_radius_in_mm);
        T width = 0.5;
        
        
        
        addTunnelsAlongRod(0, heights, scaleAndShift, 2);
        addTunnelsAlongRod(1, heights, scaleAndShift, 1);
        addTunnelsAlongRod(2, heights, scaleAndShift, 1);
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::circlePatchGCode(int n_row, int n_col, int type, bool add_bar)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 40);
    };

    if constexpr (dim == 3)
    {
        
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;

        T r = 0.008 * 1e3; 
        T d = 0.9 * 2.0 * r;

        writeHeader();
        if (type == 0)
        {   
            // TV bottom_left_center = TV(0, 0, rod_radius_in_mm);
            TV bottom_left_center = TV(40, 40, rod_radius_in_mm);

            // TV left = bottom_left_center;
            // TV start = left + TV(r, 0, 0);
            // writeCircle(left, r, start, TV2(0, M_PI * 2.0), {}, 20, rod_radius_in_mm);

            // TV right = left + TV(d, 0, 0);
            // start = right + TV(r, 0, 0);
            // TV ixn0, ixn1;
            // circleCircleIntersection(left, r, right, r, ixn0, ixn1);
            // writeCircle(right, r, start, TV2(0, M_PI * 2.0), {ixn0, ixn1}, 20, rod_radius_in_mm);

            // T theta = std::acos((ixn0 - left).normalized().dot(TV(1, 0, 0)));
            // T epsilon = 0.3;
            // TV tunnel_left = left + TV(r * std::cos(theta - epsilon), r * std::sin(theta - epsilon), 0);
            // TV tunnel_right = left + TV(r * std::cos(theta + epsilon + 0.1), r * std::sin(theta + epsilon + 0.1), 0);
            // ixn0[dim-1] = 4.0 * rod_radius_in_mm;
            // addSingleTunnel(tunnel_left, ixn0, tunnel_right, rod_radius_in_mm);


            for (int row = 0; row < n_row; row+=2)
            {
                for (int col = 0; col < n_col; col+=2)
                {
                    TV center = bottom_left_center + TV(40 + d * col, 40 + d * row, 0);
                    TV start = center + TV(r, 0, 0);
                    if (col == 0)
                    {
                        start = center - TV(0, r, 0);
                        writeCircle(center, r, start, TV2(-M_PI / 2.0, M_PI / 2.0), {}, 10, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0, M_PI * 1.5), {}, 10, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0), {}, 20, rod_radius_in_mm);
                }   
            }

            for (int row = 0; row < n_row; row++)
            {
                int col = row % 2 == 0 ? 1 : 0;
                int col_end = row % 2 == 0 ? n_col - 1 : n_col;
                for (; col < col_end; col+=2)
                {
                    TV center = bottom_left_center + TV(40 + d * col, 40 + d * row, 0);
                    std::vector<TV> ixns;
                    if (row % 2 == 0)
                    {
                        TV left_center = bottom_left_center + TV(40 + d * (col - 1), 40 + d * row, 0);
                        TV right_center = bottom_left_center + TV(40 + d * (col + 1), 40 + d * row, 0);
                        TV ixn_left0, ixn_left1, ixn_right0, ixn_right1;
                        circleCircleIntersection(left_center, r, center, r, ixn_left0, ixn_left1);
                        circleCircleIntersection(right_center, r, center, r, ixn_right0, ixn_right1);

                        ixns = { ixn_left0, ixn_left1, ixn_right1, ixn_right0 };
                    }
                    if (row % 2 == 1)
                    {
                        TV bottom_center = bottom_left_center + TV(40 + d * col, 40 + d * (row - 1), 0);
                        TV top_center = bottom_left_center + TV(40 + d * col, 40 + d * (row + 1), 0);
                        TV ixn_bottom0, ixn_bottom1, ixn_top0, ixn_top1;
                        circleCircleIntersection(bottom_center, r, center, r, ixn_bottom0, ixn_bottom1);
                        circleCircleIntersection(top_center, r, center, r, ixn_top0, ixn_top1);

                        ixns = { ixn_bottom0, ixn_bottom1, ixn_top0, ixn_top1};
                    }
                    TV start = center - TV(0, r, 0);
                    if (col == 0)
                    {
                        writeCircle(center, r, start, TV2(-M_PI / 2.0, M_PI / 2.0), ixns, 20, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0, M_PI * 1.5), ixns, 20, rod_radius_in_mm);
                    }
                    else
                    {
                        writeCircle(center, r, start, TV2(-M_PI / 2.0, M_PI * 1.5), ixns, 40, rod_radius_in_mm);
                    }

                    if (row % 2 == 0)
                    {
                        TV ixn0 = ixns[0];
                        TV left_center = bottom_left_center + TV(40 + d * (col - 1), 40 + d * row, 0);
                        TV right_center = bottom_left_center + TV(40 + d * (col + 1), 40 + d * row, 0);

                        T theta = std::acos((ixn0 - left_center).normalized().dot(TV(1, 0, 0)));
                        T epsilon_left = 0.25, epsilon_right = 0.25;

                        auto writeTunnel = [&](T _theta, const TV& _ixn, const TV& _center)
                        {
                            TV ixn = _ixn;
                            TV tunnel_left = _center + TV(r * std::cos(_theta - epsilon_left), r * std::sin(_theta - epsilon_left), 0);
                            TV tunnel_right = _center + TV(r * std::cos(_theta + epsilon_right), r * std::sin(_theta + epsilon_right), 0);
                            ixn[dim-1] = 5.0 * rod_radius_in_mm;
                            addSingleTunnel(tunnel_left, ixn, tunnel_right, rod_radius_in_mm);
                        };
                        
                        writeTunnel(theta, ixns[0], left_center);
                        writeTunnel(M_PI + M_PI - theta, ixns[1], left_center);
                        
                        writeTunnel(M_PI - theta, ixns[2], right_center);
                        writeTunnel(-M_PI + theta, ixns[3], right_center);
                    }
                    else if (row % 2 == 1)
                    {
                        TV ixn0 = ixns[2];

                        TV bottom_center = bottom_left_center + TV(40 + d * col, 40 + d * (row - 1), 0);
                        TV top_center = bottom_left_center + TV(40 + d * col, 40 + d * (row + 1), 0);

                        T theta = std::acos((ixn0 - center).normalized().dot(TV(1, 0, 0)));
                        T epsilon_left = 0.25, epsilon_right = 0.25;

                        auto writeTunnel = [&](T _theta, const TV& _ixn, const TV& _center)
                        {
                            TV ixn = _ixn;
                            TV tunnel_left = _center + TV(r * std::cos(_theta - epsilon_left), r * std::sin(_theta - epsilon_left), 0);
                            TV tunnel_right = _center + TV(r * std::cos(_theta + epsilon_right), r * std::sin(_theta + epsilon_right), 0);
                            ixn[dim-1] = 5.0 * rod_radius_in_mm;
                            addSingleTunnel(tunnel_left, ixn, tunnel_right, rod_radius_in_mm);
                        };
                        if (col != n_col - 1)
                        {
                            writeTunnel(theta, ixns[1], bottom_center);
                            writeTunnel(-theta, ixns[2], top_center);
                        }
                        if (col != 0)
                        {
                            writeTunnel(-M_PI + theta, ixns[3], top_center);
                            writeTunnel(M_PI - theta, ixns[0], bottom_center);
                        }
                            
                        
                    }
                    if (col == n_col - 1)
                    {
                        current_position[dim - 1] = 2.0;
                        moveTo(current_position);
                    }
                    
                }   
            }
            for (int row = 1; row < n_row - 1; row+=2)
            {
                for (int col = 1; col < n_col; col+=2)
                {
                    TV center = bottom_left_center + TV(40 + d * col, 40 + d * row, 0);
                    if (row == 1)
                    {
                        center[dim - 1] = 2.0;
                        moveTo(center);
                        center[dim - 1] = rod_radius_in_mm;
                    }
                    TV start = center + TV(r, 0, 0);

                    TV left_center = bottom_left_center + TV(40 + d * (col - 1), 40 + d * row, 0);
                    TV right_center = bottom_left_center + TV(40 + d * (col + 1), 40 + d * row, 0);
                    TV ixn_left0, ixn_left1, ixn_right0, ixn_right1;
                    circleCircleIntersection(left_center, r, center, r, ixn_left0, ixn_left1);
                    circleCircleIntersection(center, r, right_center, r, ixn_right0, ixn_right1);
                    TV bottom_center = bottom_left_center + TV(40 + d * col, 40 + d * (row - 1), 0);
                    TV top_center = bottom_left_center + TV(40 + d * col, 40 + d * (row + 1), 0);
                    TV ixn_bottom0, ixn_bottom1, ixn_top0, ixn_top1;
                    circleCircleIntersection(bottom_center, r, center, r, ixn_bottom0, ixn_bottom1);
                    circleCircleIntersection(center, r, top_center, r, ixn_top0, ixn_top1);

                    std::vector<TV> ixns = {ixn_right0, ixn_top1, ixn_top0, ixn_left0, ixn_left1, ixn_bottom0, ixn_bottom1, ixn_right1};

                    T theta = std::acos((ixns[0] - center).normalized().dot(TV(1, 0, 0)));
                    T epsilon_left = 0.25, epsilon_right = 0.25;

                    auto writeTunnel = [&](T _theta, const TV& _ixn, const TV& _center)
                    {
                        TV ixn = _ixn;
                        TV tunnel_left = _center + TV(r * std::cos(_theta - epsilon_left), r * std::sin(_theta - epsilon_left), 0);
                        TV tunnel_right = _center + TV(r * std::cos(_theta + epsilon_right), r * std::sin(_theta + epsilon_right), 0);
                        ixn[dim-1] = 5.0 * rod_radius_in_mm;
                        addSingleTunnel(tunnel_left, ixn, tunnel_right, rod_radius_in_mm);
                    };

                    writeCircle(center, r, start, TV2(0, M_PI * 2.0), ixns, 40, rod_radius_in_mm);
                    
                    writeTunnel(M_PI - theta, ixns[0], right_center);
                    epsilon_left = 0.1, epsilon_right = 0.25;
                    writeTunnel(-(M_PI * 0.5 - theta), ixns[1], top_center);

                    writeTunnel(-M_PI * 0.5 - theta, ixns[2], top_center);
                    writeTunnel(theta, ixns[3], left_center);
                    writeTunnel(-theta, ixns[4], left_center);
                    writeTunnel((M_PI * 0.5 + theta), ixns[5], bottom_center);
                    writeTunnel(M_PI * 0.5 - theta, ixns[6], bottom_center);
                    epsilon_left = 0.25, epsilon_right = 0.25;
                    writeTunnel(M_PI + theta, ixns[7], right_center);
                }
            }
            // TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
            // TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            // addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
            // TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
            // TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            // addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
        }
        // current one overlay the last one
        else if (type == 1)
        {
            d = 0.9 * 2.0 * r;
            if (add_bar)
            {
                TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
                TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
                bottom_right[dim - 1] = 2.0;
                moveTo(bottom_right);
                bottom_right[dim - 1] = rod_radius_in_mm;
                addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
                TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
                TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
                top_left[dim - 1] = 2.0;
                moveTo(top_left);
                top_left[dim - 1] = rod_radius_in_mm;
                addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
            }
            TV bottom_left_center = TV(40, 40, rod_radius_in_mm);
            bool intersect_bottom = true, intersect_left = false;
            for (int row = 0; row < n_row; row++)
            {
                if (row == 0)
                    intersect_bottom = false;
                else
                    intersect_bottom = true;
                for (int col = 0; col < n_col; col++)   
                {
                    if (col == 0)
                    {
                        intersect_left = false;
                    }
                    else
                        intersect_left = true;
                    std::vector<TV> ixns;
                    TV center = bottom_left_center + TV(d * col, d * row, 0);
                    TV start = center + TV(r, 0 , 0);
                    
                    if (intersect_left)
                    {
                        TV ixn0, ixn1;
                        TV center_left = bottom_left_center + TV(d * (col - 1), d * row, 0);
                        circleCircleIntersection(center_left, r, center, r, ixn0, ixn1);
                        ixns.push_back(ixn0);
                        ixns.push_back(ixn1);
                    }
                    if (intersect_bottom)
                    {
                        TV ixn0, ixn1;
                        TV center_bottom = bottom_left_center + TV(d * col, d * (row - 1), 0);
                        circleCircleIntersection(center_bottom, r, center, r, ixn0, ixn1);
                        ixns.push_back(ixn0);
                        if (col < n_col - 1)
                            ixns.push_back(ixn1);
                    }
                    if (col == 0 && add_bar)
                    {
                        start = center - TV(0, r, 0);
                        center[dim - 1] = 2.0;
                        moveTo(center);
                        center[dim - 1] = rod_radius_in_mm;
                        writeCircle(center, r, start, TV2(-M_PI / 2.0 - 0.2, M_PI / 2.0 + 0.2), {}, 15, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1 && add_bar)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0 - 0.2, M_PI * 1.5 + 0.2), ixns, 15, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.2), ixns, 30, rod_radius_in_mm);
                    // if ((col == 0 || col == n_col - 1)  && add_bar)
                    if (col == 0  && add_bar)
                        continue;
                    T epsilon_left = 0.25, epslion_right = 0.2;
                    for (int i = 0; i < ixns.size(); i++)
                    {
                        TV ixn_center = i < 2 && intersect_left ? center - TV(d, 0, 0) : center - TV(0, d, 0);
                        if (!intersect_left)
                        {
                            ixn_center = center - TV(0, d, 0);
                            if (i == 1)
                            {
                                epsilon_left = 0.2;
                                epslion_right = 0.3;
                            }
                        }
                        TV ixn0 = ixns[i];
                        T theta = std::acos((ixn0 - ixn_center).normalized().dot(TV(1, 0, 0)));
                        if (intersect_left)
                        {
                            if (i == 1)
                            {
                                theta = M_PI + (M_PI - theta);
                                epsilon_left = 0.1;
                                epslion_right = 0.3;
                            }
                            else if (i == 2)
                            {
                                epsilon_left = 0.2;
                                epslion_right = 0.3;
                            }
                            else if (i == 3)
                            {
                                epsilon_left = 0.2;
                                epslion_right = 0.3;
                            }
                        }
                        
                        TV tunnel_left = ixn_center + TV(r * std::cos(theta - epsilon_left), r * std::sin(theta - epsilon_left), 0);
                        TV tunnel_right = ixn_center + TV(r * std::cos(theta + epslion_right), r * std::sin(theta + epslion_right), 0);
                        ixn0[dim-1] = 5.0 * rod_radius_in_mm;
                        addSingleTunnel(tunnel_left, ixn0, tunnel_right, rod_radius_in_mm);
                    }
                    if (col == n_col - 1)
                    {
                        current_position[dim - 1] = 2.0;
                        moveTo(current_position);
                    }
                }
            }
            // TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
            // TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            // bottom_right[dim - 1] = 2.0;
            // moveTo(bottom_right);
            // bottom_right[dim - 1] = rod_radius_in_mm;
            // addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
            // TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
            // TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            // top_left[dim - 1] = 2.0;
            // moveTo(top_left);
            // top_left[dim - 1] = rod_radius_in_mm;
            // addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
        }
        // fused
        else if (type == 2)
        {
            d = 0.9 * 2.0 * r;
            TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
            TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            bottom_right[dim - 1] = 2.0;
            moveTo(bottom_right);
            bottom_right[dim - 1] = rod_radius_in_mm;
            addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
            TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
            TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            top_left[dim - 1] = 2.0;
            moveTo(top_left);
            top_left[dim - 1] = rod_radius_in_mm;
            addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
            TV bottom_left_center = TV(40, 40, rod_radius_in_mm);
            bool intersect_bottom = true, intersect_left = false;
            for (int row = 0; row < n_row; row++)
            {
                for (int col = 0; col < n_col; col++)   
                {
                    std::vector<TV> ixns;
                    TV center = bottom_left_center + TV(d * col, d * row, 0);
                    TV start = center + TV(r, 0 , 0);
                    if (col == 0 && add_bar)
                    {
                        start = center - TV(0, r, 0);
                        writeCircle(center, r, start, TV2(-M_PI / 2.0 - 0.2, M_PI / 2.0 + 0.2), {}, 15, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1 && add_bar)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0 - 0.2, M_PI * 1.5 + 0.2), {}, 15, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.2), ixns, 30, rod_radius_in_mm);
                }
            }
        }
        else if (type == 3)
        {
            TV bottom_left_center = TV(40, 40, rod_radius_in_mm);

            for (int row = 0; row < n_row; row++)
            {
                int col = (row % 2 == 0) ? 0 : 1;
                int col_end = (row % 2 == 0) ? n_col : n_col - 1;
                for (; col < col_end; col+=2)
                {
                    TV center = bottom_left_center + TV(0 + d * col, 0 + d * row, 0);
                    TV start = center + TV(r, 0, 0);
                    if (col == 0)
                    {
                        start = center - TV(0, r, 0);
                        writeCircle(center, r, start, TV2(-M_PI / 2.0 - 0.2, M_PI / 2.0 + 0.2), {}, 10, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0 - 0.2, M_PI * 1.5 + 0.2), {}, 10, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.2), {}, 20, rod_radius_in_mm);
                }   
            }

            for (int row = 0; row < n_row; row++)
            {
                int col = (row % 2 == 0) ? 1 : 0;
                int col_end = (row % 2 == 0) ? n_col - 1 : n_col;
                for (; col < col_end; col+=2)
                {
                    TV center = bottom_left_center + TV(d * col, d * row, 0);
                    TV start = center + TV(r, 0, 0);

                    TV left_center = bottom_left_center + TV(d * (col - 1),  d * row, 0);
                    TV right_center = bottom_left_center + TV(d * (col + 1), d * row, 0);
                    TV ixn_left0, ixn_left1, ixn_right0, ixn_right1;
                    TV bottom_center = bottom_left_center + TV( d * col, d * (row - 1), 0);
                    TV top_center = bottom_left_center + TV(d * col,d * (row + 1), 0);
                    TV ixn_bottom0, ixn_bottom1, ixn_top0, ixn_top1;
                    
                    circleCircleIntersection(bottom_center, r, center, r, ixn_bottom0, ixn_bottom1);
                    circleCircleIntersection(center, r, top_center, r, ixn_top0, ixn_top1);
                    circleCircleIntersection(left_center, r, center, r, ixn_left0, ixn_left1);
                    circleCircleIntersection(center, r, right_center, r, ixn_right0, ixn_right1);
                    
                    std::vector<TV> ixns = {ixn_right0, ixn_top1, ixn_top0, ixn_left0, ixn_left1, ixn_bottom0, ixn_bottom1, ixn_right1};

                    std::vector<TV> ixn_write;

                    if (col > 0)
                    {
                        ixn_write.push_back(ixn_left0);
                        ixn_write.push_back(ixn_left1);
                    }
                    if (col < n_col-1)
                    {
                        ixn_write.push_back(ixn_right0);
                        ixn_write.push_back(ixn_right1);
                    }
                    if (row > 0)
                    {
                        ixn_write.push_back(ixn_bottom0);
                        ixn_write.push_back(ixn_bottom1);
                    }
                    if (row < n_row-1)
                    {
                        ixn_write.push_back(ixn_top0);
                        ixn_write.push_back(ixn_top1);
                    }
                        
                    T theta = std::acos((ixns[0] - center).normalized().dot(TV(1, 0, 0)));
                    
                    T epsilon_left = 0.25, epsilon_right = 0.25;

                    auto writeTunnel = [&](T _theta, const TV& _ixn, const TV& _center)
                    {
                        TV ixn = _ixn;
                        TV tunnel_left = _center + TV(r * std::cos(_theta - epsilon_left), r * std::sin(_theta - epsilon_left), 0);
                        TV tunnel_right = _center + TV(r * std::cos(_theta + epsilon_right), r * std::sin(_theta + epsilon_right), 0);
                        ixn[dim-1] = 5.0 * rod_radius_in_mm;
                        addSingleTunnel(tunnel_left, ixn, tunnel_right, rod_radius_in_mm);
                    };
                    if (col == 0)
                    {
                        start = center - TV(0, r, 0);
                        writeCircle(center, r, start, TV2(-M_PI / 2.0, M_PI / 2.0), ixn_write, 20, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0, M_PI * 1.5), ixn_write, 20, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.1), ixn_write, 40, rod_radius_in_mm);

                    if (col < n_col - 1)
                    {
                        writeTunnel(M_PI - theta, ixns[0], right_center);
                    }
                    if (row < n_row - 1 && col < n_col - 1)
                    {
                        writeTunnel(-(M_PI * 0.5 - theta), ixns[1], top_center);    
                    }
                    if (row < n_row - 1 && col > 0)
                        writeTunnel(-M_PI * 0.5 - theta, ixns[2], top_center);
                    if (col > 0)
                    {
                        writeTunnel(theta, ixns[3], left_center);
                        writeTunnel(-theta, ixns[4], left_center);
                    }
                    if (row > 0 && col > 0)
                        writeTunnel((M_PI * 0.5 + theta), ixns[5], bottom_center);
                    if (row > 0 && col < n_col - 1)
                        writeTunnel(M_PI * 0.5 - theta, ixns[6], bottom_center);
                    if (col < n_col - 1)
                        writeTunnel(M_PI + theta, ixns[7], right_center);
                }
            }

        }
        else if (type == 4)
        {
            d = 0.8 * 2.0 * r;
            TV bottom_left_center = TV(40, 40, rod_radius_in_mm);
            TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
            TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            bottom_right[dim - 1] = 2.0;
            moveTo(bottom_right);
            bottom_right[dim - 1] = rod_radius_in_mm;
            addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
            TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
            TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            top_left[dim - 1] = 2.0;
            moveTo(top_left);
            top_left[dim - 1] = rod_radius_in_mm;
            addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
            

            for (int row = 0; row < n_row; row++)
            {
                int col = (row % 2 == 0) ? 0 : 1;
                int col_end = (row % 2 == 0) ? n_col : n_col - 1;
                for (; col < col_end; col+=2)
                {
                    TV center = bottom_left_center + TV(0 + d * col, 0 + d * row, 0);
                    TV start = center + TV(r, 0, 0);
                    if (col == 0)
                    {
                        start = center - TV(0, r, 0);
                        writeCircle(center, r, start, TV2(-M_PI / 2.0 - 0.2, M_PI / 2.0 + 0.2), {}, 10, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0 - 0.2, M_PI * 1.5 + 0.2), {}, 10, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.2), {}, 20, rod_radius_in_mm);
                }   
            }

            for (int row = 0; row < n_row; row++)
            {
                int col = (row % 2 == 0) ? 1 : 0;
                int col_end = (row % 2 == 0) ? n_col - 1 : n_col;
                for (; col < col_end; col+=2)
                {
                    TV center = bottom_left_center + TV(d * col, d * row, 0);
                    TV start = center + TV(r, 0, 0);

                    TV left_center = bottom_left_center + TV(d * (col - 1),  d * row, 0);
                    TV right_center = bottom_left_center + TV(d * (col + 1), d * row, 0);
                    TV ixn_left0, ixn_left1, ixn_right0, ixn_right1;
                    TV bottom_center = bottom_left_center + TV( d * col, d * (row - 1), 0);
                    TV top_center = bottom_left_center + TV(d * col,d * (row + 1), 0);
                    TV ixn_bottom0, ixn_bottom1, ixn_top0, ixn_top1;
                    
                    circleCircleIntersection(bottom_center, r, center, r, ixn_bottom0, ixn_bottom1);
                    circleCircleIntersection(center, r, top_center, r, ixn_top0, ixn_top1);
                    circleCircleIntersection(left_center, r, center, r, ixn_left0, ixn_left1);
                    circleCircleIntersection(center, r, right_center, r, ixn_right0, ixn_right1);
                    
                    std::vector<TV> ixns = {ixn_right0, ixn_top1, ixn_top0, ixn_left0, ixn_left1, ixn_bottom0, ixn_bottom1, ixn_right1};

                    // std::vector<TV> ixns = {ixn_top0, ixn_left0, ixn_left1, ixn_bottom0};

                    std::vector<TV> ixn_write;

                    

                    if (col > 0)
                    {
                        ixn_write.push_back(ixn_left0);
                        ixn_write.push_back(ixn_left1);
                    }
                    if (col < n_col-1)
                    {
                        // ixn_write.push_back(ixn_right0);
                        // ixn_write.push_back(ixn_right1);
                    }
                    if (row > 0)
                    {
                        ixn_write.push_back(ixn_bottom0);
                        // ixn_write.push_back(ixn_bottom1);
                    }
                    if (row < n_row-1)
                    {
                        ixn_write.push_back(ixn_top0);
                        // ixn_write.push_back(ixn_top1);
                    }
                        
                    T theta = std::acos((ixns[0] - center).normalized().dot(TV(1, 0, 0)));
                    
                    T epsilon_left = 0.25, epsilon_right = 0.25;
                    

                    auto writeTunnel = [&](T _theta, const TV& _ixn, const TV& _center)
                    {
                        TV ixn = _ixn;
                        TV tunnel_left = _center + TV(r * std::cos(_theta - epsilon_left), r * std::sin(_theta - epsilon_left), 0);
                        TV tunnel_right = _center + TV(r * std::cos(_theta + epsilon_right), r * std::sin(_theta + epsilon_right), 0);
                        ixn[dim-1] = 5.0 * rod_radius_in_mm;
                        addSingleTunnel(tunnel_left, ixn, tunnel_right, rod_radius_in_mm);
                    };
                    if (col == 0)
                    {
                        start = center - TV(0, r, 0);
                        writeCircle(center, r, start, TV2(-M_PI / 2.0, M_PI / 2.0), ixn_write, 20, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0, M_PI * 1.5), ixn_write, 20, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.1), ixn_write, 40, rod_radius_in_mm);

                    if (col < n_col - 1)
                    {
                        // writeTunnel(M_PI - theta, ixns[0], right_center);
                    }
                    if (row < n_row - 1 && col < n_col - 1)
                    {
                        // writeTunnel(-(M_PI * 0.5 - theta), ixns[1], top_center);    
                    }
                    if (row < n_row - 1 && col > 0)
                        writeTunnel(-M_PI * 0.5 - theta, ixns[2], top_center);
                    if (col > 0)
                    {
                        writeTunnel(theta, ixns[3], left_center);
                        writeTunnel(-theta, ixns[4], left_center);
                    }
                    if (row > 0 && col > 0)
                        writeTunnel((M_PI * 0.5 + theta), ixns[5], bottom_center);
                    // if (row > 0 && col < n_col - 1)
                        // writeTunnel(M_PI * 0.5 - theta, ixns[6], bottom_center);
                    // if (col < n_col - 1)
                        // writeTunnel(M_PI + theta, ixns[7], right_center);
                }
            }

        }
        else if (type == 5)
        {
            d = 0.8 * 2.0 * r;
            if (add_bar)
            {
                TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
                TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
                bottom_right[dim - 1] = 2.0;
                moveTo(bottom_right);
                bottom_right[dim - 1] = rod_radius_in_mm;
                addRecBar(bottom_right, top_right, 2, rod_radius_in_mm);
                TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
                TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
                top_left[dim - 1] = 2.0;
                moveTo(top_left);
                top_left[dim - 1] = rod_radius_in_mm;
                addRecBar(top_left, bottom_left, 2, rod_radius_in_mm);
            }
            TV bottom_left_center = TV(40, 40, rod_radius_in_mm);
            bool intersect_bottom = true, intersect_left = false;
            for (int row = 0; row < n_row; row++)
            {
                if (row == 0)
                    intersect_bottom = false;
                else
                    intersect_bottom = true;
                for (int col = 0; col < n_col; col++)   
                {
                    if (col == 0)
                    {
                        intersect_left = false;
                    }
                    else
                        intersect_left = true;
                    std::vector<TV> ixns;
                    TV center = bottom_left_center + TV(d * col, d * row, 0);
                    TV start = center + TV(r, 0 , 0);
                    if (col == 0)
                    {
                        start[dim - 1] = 2.0;
                        moveTo(start);
                        start[dim - 1] = rod_radius_in_mm;
                    }
                    if (intersect_left)
                    {
                        // TV ixn0, ixn1;
                        // TV center_left = bottom_left_center + TV(d * (col - 1), d * row, 0);
                        // circleCircleIntersection(center_left, r, center, r, ixn0, ixn1);
                        // ixns.push_back(ixn0);
                        // ixns.push_back(ixn1);
                    }
                    if (intersect_bottom)
                    {
                        TV ixn0, ixn1;
                        TV center_bottom = bottom_left_center + TV(d * col, d * (row - 1), 0);
                        circleCircleIntersection(center_bottom, r, center, r, ixn0, ixn1);
                        ixns.push_back(ixn0);
                        if (col < n_col - 1)
                            ixns.push_back(ixn1);
                    }
                    if (col == 0 && add_bar)
                    {
                        start = center - TV(0, r, 0);
                        center[dim - 1] = 2.0;
                        moveTo(center);
                        center[dim - 1] = rod_radius_in_mm;
                        writeCircle(center, r, start, TV2(-M_PI / 2.0 - 0.2, M_PI / 2.0 + 0.2), {}, 15, rod_radius_in_mm);
                    }
                    else if (col == n_col - 1 && add_bar)
                    {
                        start = center + TV(0, r, 0);
                        writeCircle(center, r, start, TV2(M_PI / 2.0 - 0.2, M_PI * 1.5 + 0.2), ixns, 15, rod_radius_in_mm);
                    }
                    else
                        writeCircle(center, r, start, TV2(0, M_PI * 2.0 + 0.2), ixns, 30, rod_radius_in_mm);
                    // if ((col == 0 || col == n_col - 1)  && add_bar)
                    if (col == 0  && add_bar)
                        continue;
                    T epsilon_left = 0.25, epslion_right = 0.2;
                    if (ixns.size())
                    {
                        TV ixn_center = center - TV(0, d, 0);
                        TV ixn0 = ixns[0];
                        T theta = std::acos((ixn0 - ixn_center).normalized().dot(TV(1, 0, 0)));
                        TV tunnel_left = ixn_center + TV(r * std::cos(theta - epsilon_left), r * std::sin(theta - epsilon_left), 0);
                        TV tunnel_right = ixn_center + TV(r * std::cos(theta + epslion_right), r * std::sin(theta + epslion_right), 0);
                        ixn0[dim-1] = 5.0 * rod_radius_in_mm;
                        addSingleTunnel(tunnel_left, ixn0, tunnel_right, rod_radius_in_mm);
                        if (ixns.size() > 1)
                        {
                            TV ixn0 = ixns[1];
                            T theta = std::acos((ixn0 - ixn_center).normalized().dot(TV(1, 0, 0)));
                            TV tunnel_left = ixn_center + TV(r * std::cos(theta - epsilon_left), r * std::sin(theta - epsilon_left), 0);
                            TV tunnel_right = ixn_center + TV(r * std::cos(theta + epslion_right), r * std::sin(theta + epslion_right), 0);
                            ixn0[dim-1] = 4.0 * rod_radius_in_mm;
                            addSingleTunnel(tunnel_left, ixn0, tunnel_right, rod_radius_in_mm);
                        }
                    }
                    if (col == n_col - 1)
                    {
                        current_position[dim - 1] = 2.0;
                        moveTo(current_position);
                    }
                }
            }
            // TV bottom_right = TV(40 + d * (n_col - 1), 40 -r, rod_radius_in_mm);
            // TV top_right = TV(40 + d * (n_col - 1), 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            // bottom_right[dim - 1] = 2.0;
            // moveTo(bottom_right);
            // bottom_right[dim - 1] = rod_radius_in_mm;
            // addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
            // TV bottom_left = TV(40, 40 -r, rod_radius_in_mm);
            // TV top_left = TV(40, 40 + d * (n_row - 1) + r, rod_radius_in_mm);
            // top_left[dim - 1] = 2.0;
            // moveTo(top_left);
            // top_left[dim - 1] = rod_radius_in_mm;
            // addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
        }
        
        
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::slidingBlocksGCode(int n_row, int n_col, int type, bool add_bar)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 40);
    };

    if constexpr (dim == 3)
    {
        TV bottom_left, top_right;
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
        sim.computeBoundingBox(bottom_left, top_right);

        scaleAndShift(bottom_left); scaleAndShift(top_right);
        bottom_left[dim-1] = rod_radius_in_mm;
        top_right[dim-1] = rod_radius_in_mm;

        TV bottom_right = bottom_left;
        bottom_right[0] = top_right[0];
        TV top_left = top_right;
        top_left[0] = bottom_left[0];   

        writeHeader();
        T extend_percentage = 0.05;
        T inner_height = 3.5 * rod_radius_in_mm;
        // x sliding
        if (type == 0)
        {
            int rod_cnt = 0; 
            int boundary_rod_cnt = 0;
            // outer contour
            

            if (add_bar)
            {
                for (int col = 0; col < n_col; col++)
                    generateCodeSingleRod(col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);   
                for (int col = 0; col < n_col; col++)
                {
                    if (col == 0)
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.1);
                    else
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
                    
                
                addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
                addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
                
                rod_cnt += 2 * n_col + 2 * n_row;
            }
            else
            {
                for (int col = 0; col < n_col; col++)
                {
                    generateCodeSingleRod(col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
                for (int row = 0; row < n_row; row++)
                {
                    generateCodeSingleRod(n_col * 2 + n_row + row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
                for (int col = 0; col < n_col; col++)
                {
                    generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
                for (int row = 0; row < n_row; row++)
                {
                    generateCodeSingleRod(n_col * 2 + row , scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
            }

            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_col - 1) * n_row)
            {
                generateCodeSingleRod(rod_cnt, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                rod_cnt++;
            }
            
            int temp = rod_cnt + 2 * (n_row - 1) * n_col;
            for (int col = 0; col < n_col - 1; col++)
            {
                for (int row = 0; row < n_row - 1; row++)
                {
                    int idx = row * (n_col - 1) + col;
                    if (row == n_row - 2)
                        generateCodeSingleRod(temp + idx * 4 + 1, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, false, true);
                    else
                        generateCodeSingleRod(temp + idx * 4 + 1, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
            }
            for (int col = 0; col < n_col - 1; col++)
            {
                for (int row = n_row - 2; row > -1; row--)
                {
                    int idx = row * (n_col - 1) + col;
                    if (row == 0)
                        generateCodeSingleRod(temp + idx * 4 + 3, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, false, true);
                    else
                        generateCodeSingleRod(temp + idx * 4 + 3, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
            }
            
            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_row - 1) * n_col)
                generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage, true);


            
            temp = rod_cnt;
            for (int col = 0; col < n_col - 1; col++)
            {
                for (int row = 0; row < n_row - 1; row++)
                {
                    int idx = row * (n_col - 1) + col;
                    if (row == 0)
                    {
                        auto rod = sim.Rods[temp + idx * 4 + 0];
                        TV x0; rod->x(rod->indices.front(), x0);
                        TV front, back;
                        rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);
                        x0 -= (back - front).normalized() * 0.3 * (front - back).norm();
                        scaleAndShift(x0);
                        x0[dim - 1] = 2.0;
                        moveTo(x0);
                    }
                    generateCodeSingleRod(temp + idx * 4 + 0, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage, true);
                    if (row == n_row - 2)
                    {
                        current_position[dim - 1] = 2.0;
                        moveTo(current_position, 300, true);
                    }
                }
            }

            for (int row = n_row - 2; row > -1; row--)
            {
                for (int col = n_col - 2; col > -1; col--)
                {
                    int idx = row * (n_col - 1) + col;
                    if (col == n_col - 2)
                    {
                        auto rod = sim.Rods[temp + idx * 4 + 2];
                        TV x0; rod->x(rod->indices.front(), x0);
                        TV front, back;
                        rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);
                        x0 -= (back - front).normalized() * 0.3 * (front - back).norm();
                        scaleAndShift(x0);
                        x0[dim - 1] = 2.0;
                        moveTo(x0);
                    }
                    generateCodeSingleRod(temp + idx * 4 + 2, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage, true);
                    if (col == 0)
                    {
                        current_position[dim - 1] = 2.0;
                        moveTo(current_position, 300, true);
                    }
                }
            }

            TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            T width;
            if (sim.unit == 0.05)
                width = 0.1;
            tunnel_height = 0.2;
            for (int row = 0; row < n_row - 1; row++)
            {
                for (int col = 0; col < n_col - 1; col++)
                {
                    int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
                    for (int corner = 0; corner < 4; corner++)
                    {
                        int crossing_id = base + corner * 3 + 1;
                        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.08, 0.12), 0.3, 150, 200);
                        crossing_id = base + corner * 3 + 2;
                        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 1, scaleAndShift, Range(0.08, 0.12), 0.3, 150, 200);
                    }
                }
            }
        }
        else if (type == 1)
        {
            int rod_cnt = 0; 
            int boundary_rod_cnt = 0;

            if (add_bar)
            {
                for (int col = 0; col < n_col; col++)
                    generateCodeSingleRod(col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);   
                for (int col = 0; col < n_col; col++)
                {
                    if (col == 0)
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.1);
                    else
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
                
                addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
                addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
                
                rod_cnt += 2 * n_col + 2 * n_row;
            }

            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_col - 1) * n_row)
            {
                generateCodeSingleRod(rod_cnt, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                rod_cnt++;
            }

            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_row - 1) * n_col)
                generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, true);

            for (int row = 0; row < n_row - 1; row++)
            {
                for (int col = 0; col < n_col - 1; col++)
                {
                    // if (row == 0 && col == 0)
                    //     generateCodeSingleRod(rod_cnt++, scaleAndShift, true, 0.5 * rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.0);
                    // else
                    auto rod = sim.Rods[rod_cnt];
                    TV x0; rod->x(rod->indices.front(), x0);
                    TV front, back;
                    rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);
                    x0 -= (back - front).normalized() * 0.3 * (front - back).norm();
                    scaleAndShift(x0);
                    x0[dim - 1] = 2.0;
                    moveTo(x0);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage, false, true);
                    current_position[dim - 1] = 2.0;
                    moveTo(current_position, 300, true);
                }
            }

            

            TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            T width;
            tunnel_height = 0.2;
            for (int row = 0; row < n_row - 1; row++)
            {
                for (int col = 0; col < n_col - 1; col++)
                {
                    int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
                    for (int corner = 0; corner < 4; corner++)
                    {
                        int crossing_id = base + corner * 3 + 1;
                        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.08, 0.12), 0.3, 150, 200);
                        crossing_id = base + corner * 3 + 2;
                        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.08, 0.12), 0.3, 150, 200);
                    }
                }
            }
        }
        // fused
        else if (type == 2)
        {
            for (auto& rod : sim.Rods)
                rod->fixed_by_crossing = std::vector<bool>(rod->dof_node_location.size(), true);

            int rod_cnt = 0; 
            int boundary_rod_cnt = 0;
            // outer contour
            
            if (add_bar)
            {
                for (int col = 0; col < n_col; col++)
                    generateCodeSingleRod(col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);   
                for (int col = 0; col < n_col; col++)
                {
                    if (col == 0)
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.1);
                    else
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
                
                addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
                addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
                
                rod_cnt += 2 * n_col + 2 * n_row;
            }
            else
            {
                for (int col = 0; col < n_col; col++)
                {
                    generateCodeSingleRod(col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
                for (int row = 0; row < n_row; row++)
                {
                    generateCodeSingleRod(n_col * 2 + n_row + row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
                for (int col = 0; col < n_col; col++)
                {
                    generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
                for (int row = 0; row < n_row; row++)
                {
                    generateCodeSingleRod(n_col * 2 + row , scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                    rod_cnt++;
                }
            }

            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_col - 1) * n_row)
            {
                generateCodeSingleRod(rod_cnt, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                rod_cnt++;
            }
            
            int temp = rod_cnt + 2 * (n_row - 1) * n_col;
            for (int col = 0; col < n_col - 1; col++)
            {
                for (int row = 0; row < n_row - 1; row++)
                {
                    int idx = row * (n_col - 1) + col;
                    if (row == n_row - 2)
                        generateCodeSingleRod(temp + idx * 4 + 1, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, false, true);
                    else
                        generateCodeSingleRod(temp + idx * 4 + 1, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
            }
            for (int col = 0; col < n_col - 1; col++)
            {
                for (int row = n_row - 2; row > -1; row--)
                {
                    int idx = row * (n_col - 1) + col;
                    if (row == 0)
                        generateCodeSingleRod(temp + idx * 4 + 3, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, false, true);
                    else
                        generateCodeSingleRod(temp + idx * 4 + 3, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
            }
            
            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_row - 1) * n_col)
                generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, true);


            
            temp = rod_cnt;
            for (int col = 0; col < n_col - 1; col++)
            {
                for (int row = 0; row < n_row - 1; row++)
                {
                    int idx = row * (n_col - 1) + col;
                    generateCodeSingleRod(temp + idx * 4 + 0, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, true);
                }
            }

            for (int row = n_row - 2; row > -1; row--)
            {
                for (int col = n_col - 2; col > -1; col--)
                {
                    int idx = row * (n_col - 1) + col;
                    generateCodeSingleRod(temp + idx * 4 + 2, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, true);           
                }
            }


            
        }
        // Fix corner 1X1Y
        else if (type == 3)
        {
            for (int row = 0; row < n_row - 1; row++)
            {
                for (int col = 0; col < n_col - 1; col++)
                {
                    int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
                    for (int corner = 0; corner < 4; corner++)
                    {
                        
                        auto crossing = sim.rod_crossings[base + corner * 3 + 1];
                        

                        if (corner == 0)
                        {
                            crossing->is_fixed = true;
                            sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[1] = true;
                        }

                        crossing = sim.rod_crossings[base + corner * 3 + 2];
                    
                        if (corner == 0)
                        {
                            sim.Rods[crossing->rods_involved[1]]->fixed_by_crossing[2] = true;
                        }

                    }
                }
            }
            int rod_cnt = 0; 
            int boundary_rod_cnt = 0;

            if (add_bar)
            {
                for (int col = 0; col < n_col; col++)
                    generateCodeSingleRod(col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);   
                for (int col = 0; col < n_col; col++)
                {
                    if (col == 0)
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.1);
                    else
                        generateCodeSingleRod(n_col + col, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                }
                
                addRecBar(top_left, bottom_left, 10, rod_radius_in_mm);
                addRecBar(bottom_right, top_right, 10, rod_radius_in_mm);
                
                rod_cnt += 2 * n_col + 2 * n_row;
            }

            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_col - 1) * n_row)
            {
                generateCodeSingleRod(rod_cnt, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage);
                rod_cnt++;
            }

            boundary_rod_cnt = rod_cnt;
            while (rod_cnt < boundary_rod_cnt + 2 * (n_row - 1) * n_col)
                generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, extend_percentage, true);

            for (int row = 0; row < n_row - 1; row++)
            {
                for (int col = 0; col < n_col - 1; col++)
                {
                    // if (row == 0 && col == 0)
                    //     generateCodeSingleRod(rod_cnt++, scaleAndShift, true, 0.5 * rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.0);
                    // else
                    auto rod = sim.Rods[rod_cnt];
                    TV x0; rod->x(rod->indices.front(), x0);
                    TV front, back;
                    rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);
                    x0 -= (back - front).normalized() * 0.3 * (front - back).norm();
                    scaleAndShift(x0);
                    x0[dim - 1] = 2.0;
                    moveTo(x0);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
                    generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage, false, true);
                    current_position[dim - 1] = 2.0;
                    moveTo(current_position, 300, true);
                }
            }

            TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            T width;
            tunnel_height = 0.2;
            for (int row = 0; row < n_row - 1; row++)
            {
                for (int col = 0; col < n_col - 1; col++)
                {
                    int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
                    for (int corner = 0; corner < 4; corner++)
                    {
                        if (corner == 0)
                            continue;
                        int crossing_id = base + corner * 3 + 1;
                        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.08, 0.12), 0.3, 150, 200);
                        crossing_id = base + corner * 3 + 2;
                        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.08, 0.12), 0.3, 150, 200);
                    }
                }
            }

            // int rod_cnt = 0; 
            // for (int row = 0; row < n_row; row++)
            // {
            //     for (int col = 0; col < n_col; col++)
            //     {
            //         TV bottom_left;
            //         sim.getCrossingPosition((row * n_col + col) * 4, bottom_left);
            //         scaleAndShift(bottom_left);
            //         bottom_left[dim-1] = rod_radius_in_mm; 
            //         // moveTo(bottom_left);
            //         for (int corner = 0; corner < 4; corner++)
            //         {
            //             int from_idx = (row * n_col + col) * 4 + corner;
            //             int to_idx = (row * n_col + col) * 4 + (corner + 1) % 4;
            //             TV from, to;
            //             sim.getCrossingPosition(from_idx, from);
            //             sim.getCrossingPosition(to_idx, to);

            //             TV left, right;
            //             left = from - (to - from) * extend_percentage;
            //             right = to + (to - from) * extend_percentage;

            //             scaleAndShift(left); scaleAndShift(right);
            //             left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                        
            //             moveTo(left);
            //             writeLine(left, right, rod_radius_in_mm);
            //             // for (int i = 0; i < 8; i++)
            //             // {
            //             //     writeLine(from + (to - from) * i/8.0, from + (to - from) * (i + 1) /8.0, rod_radius_in_mm);
            //             // }
            //         }
                    
            //     }
            // }
            // for (int col = 0; col < n_col; col++) rod_cnt+=2;
            // for (int row = 0; row < n_row; row++) rod_cnt+=2;
            // for (int row = 0; row < n_row - 1; row++)
            // {
            //     for (int col = 0; col < n_col - 1; col++)
            //     {
            //         if (row == 0) rod_cnt += 2;
            //         if (row == n_row - 2) rod_cnt += 2;
            //         if (row != n_row - 2) rod_cnt += 2;
            //         if (col == 0) rod_cnt += 2;
            //         if (col == n_col - 2) rod_cnt += 2;
            //         if (n_col - 2 != col) rod_cnt += 2;
            //     }
            // }
            // for (int row = 0; row < n_row - 1; row++)
            // {
            //     for (int col = 0; col < n_col - 1; col++)
            //     {
            //         // if (row == 0 && col == 0)
            //         //     generateCodeSingleRod(rod_cnt++, scaleAndShift, true, 0.5 * rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.0);
            //         // else
            //             generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
            //         generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
            //         generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
            //         generateCodeSingleRod(rod_cnt++, scaleAndShift, true, rod_radius_in_mm, inner_height, extend_percentage);
            //     }
            // }

            // TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            // T width;
            // tunnel_height = 0.17;
            // if (sim.unit == 0.05)
            //     width = 0.12;
            // for (int row = 0; row < n_row - 1; row++)
            // {
            //     for (int col = 0; col < n_col - 1; col++)
            //     {
            //         int base = n_row * n_col * 4 + (row * (n_col - 1) + col) * 12;
            //         for (int corner = 0; corner < 4; corner++)
            //         {
            //             if (corner == 0)
            //                 continue;
            //             int crossing_id = base + corner * 3 + 1;
            //             addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(width, width), 0.2, 150, 200);
            //             crossing_id = base + corner * 3 + 2;
            //             addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(width, width), 0.2, 150, 200);
            //         }
            //     }
            // }
        }
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeFromRodsGridHardCoded(int n_row, int n_col, int type)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 40);
    };

    T rod_radius_in_mm = sim.Rods[0]->a * 1e3 * 2.0;
    
    if constexpr (dim == 3)
    {
        writeHeader();
        T extend_percentage = 0.3;
        TV bottom_left, top_right;
        sim.computeBoundingBox(bottom_left, top_right);

        scaleAndShift(bottom_left); scaleAndShift(top_right);
        bottom_left[dim-1] = rod_radius_in_mm;
        top_right[dim-1] = rod_radius_in_mm;

        TV bottom_right = bottom_left;
        bottom_right[0] = top_right[0];
        TV top_left = top_right;
        top_left[0] = bottom_left[0];

        TV bottom_left_extend = bottom_left - (bottom_right - bottom_left) * 0.2;

        
        
        if (type == 0)
        {
            for (int row = 0; row < n_row; row++)
            {
                // generateCodeSingleRod(row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            for (int col = 0; col < n_col; col++)
            {
                // generateCodeSingleRod(col + n_row, scaleAndShift, true, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[col + n_row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                scaleAndShift(from); scaleAndShift(to);
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (col % 2 == 0)
                {
                    
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
            }
            
            moveTo(bottom_left_extend);
            writeLine(bottom_left_extend, bottom_right, rod_radius_in_mm, 1500);
            writeLine(bottom_right, top_right, rod_radius_in_mm, 1500);
            writeLine(top_right, top_left, rod_radius_in_mm, 1500);
            writeLine(top_left, bottom_left, rod_radius_in_mm, 1500);
        }
        else if (type == 1)
        {
            for (int row = 0; row < n_row; row++)
            {
                // generateCodeSingleRod(row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            for (int col = 0; col < n_col; col++)
            {
                // generateCodeSingleRod(col + n_row, scaleAndShift, true, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[col + n_row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                scaleAndShift(from); scaleAndShift(to);
                from[dim-1] = 3.0 * rod_radius_in_mm; to[dim-1] = 3.0 * rod_radius_in_mm;
                
                if (col % 2 == 0)
                {
                    
                    moveTo(right);
                    right[dim-1] = 0.2;
                    moveTo(right);
                    to[1] = top_right[1];
                    TV temp = 0.5 * (right + to);
                    temp[2] = 0.2;

                    writeLine(right, temp, rod_radius_in_mm, 200);
                    writeLine(temp, to, rod_radius_in_mm, 400);
                    writeLine(to, from, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (from + left);
                    temp[2] = 0.2;
                    writeLine(from, temp, rod_radius_in_mm, 200);
                    writeLine(temp, left, rod_radius_in_mm, 200);
                }
                else
                {
                    moveTo(left);
                    left[dim-1] = 0.2;
                    moveTo(left);
                    from[1] = bottom_left[1];
                    TV temp = 0.5 * (left + from);
                    temp[2] = 0.2;
                    writeLine(left, temp, rod_radius_in_mm, 200);
                    writeLine(temp, from, rod_radius_in_mm, 400);
                    writeLine(from, to, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (to + right);
                    temp[2] = 0.2;
                    writeLine(to, temp, rod_radius_in_mm, 200);
                    writeLine(temp, right, rod_radius_in_mm, 200);
                }
            }

            TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            // T width = 0.034;
            T width = 0.015;
            tunnel_height = 0.2;
            for (int row = 0; row < n_row; row++)
            {
                for (int col = 0; col < n_col; col++)
                {
                    int crossing_id = row * n_col + col;
                    addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(width, width), 0.2, 100, 200);
                }
            }
            
            TV top_right_extend = top_right + (top_right - top_left) * 0.2;

            moveTo(top_right_extend);
            top_left[dim - 1] = 1.5 * rod_radius_in_mm;
            bottom_left[dim - 1] = 1.5 * rod_radius_in_mm;
            bottom_right[dim - 1] = 1.5 * rod_radius_in_mm;
            top_right[dim - 1] = 1.5 * rod_radius_in_mm;
            writeLine(top_right_extend, top_left, rod_radius_in_mm, 1500);
            writeLine(top_left, bottom_left, rod_radius_in_mm, 1500);
            writeLine(bottom_left, bottom_right, rod_radius_in_mm, 1500);
            writeLine(bottom_right, top_right, rod_radius_in_mm, 1500);
        }

        else if (type == 2)
        {
            TV top_right_extend = top_right + (top_right - top_left) * 0.2;

            moveTo(top_right_extend);
            writeLine(top_right_extend, top_left, rod_radius_in_mm, 1500);
            writeLine(top_left, bottom_left, rod_radius_in_mm, 1500);
            writeLine(bottom_left, bottom_right, rod_radius_in_mm, 1500);
            writeLine(bottom_right, top_right, rod_radius_in_mm, 1500);

            for (int row = 0; row < n_row; row++)
            {
                // generateCodeSingleRod(row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            for (int col = 0; col < n_col; col++)
            {
                // generateCodeSingleRod(col + n_row, scaleAndShift, true, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[col + n_row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                scaleAndShift(from); scaleAndShift(to);
                from[dim-1] = 3.0 * rod_radius_in_mm; to[dim-1] = 3.0 * rod_radius_in_mm;
                
                if (col % 2 == 0)
                {
                    
                    moveTo(right);
                    right[dim-1] = 0.2;
                    moveTo(right);
                    to[1] = top_right[1];
                    TV temp = 0.5 * (right + to);
                    temp[2] = 0.2;

                    writeLine(right, temp, rod_radius_in_mm, 200);
                    writeLine(temp, to, rod_radius_in_mm, 400);
                    writeLine(to, from, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (from + left);
                    temp[2] = 0.2;
                    writeLine(from, temp, rod_radius_in_mm, 200);
                    writeLine(temp, left, rod_radius_in_mm, 200);
                }
                else
                {
                    moveTo(left);
                    left[dim-1] = 0.2;
                    moveTo(left);
                    from[1] = bottom_left[1];
                    TV temp = 0.5 * (left + from);
                    temp[2] = 0.2;
                    writeLine(left, temp, rod_radius_in_mm, 200);
                    writeLine(temp, from, rod_radius_in_mm, 400);
                    writeLine(from, to, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (to + right);
                    temp[2] = 0.2;
                    writeLine(to, temp, rod_radius_in_mm, 200);
                    writeLine(temp, right, rod_radius_in_mm, 200);
                }
            }

        }

        else if (type == 3)
        {
            TV top_right_extend = top_right + (top_right - top_left) * 0.2;

            moveTo(top_right_extend);
            writeLine(top_right_extend, top_left, rod_radius_in_mm, 1500);
            writeLine(top_left, bottom_left, rod_radius_in_mm, 1500);
            writeLine(bottom_left, bottom_right, rod_radius_in_mm, 1500);
            writeLine(bottom_right, top_right, rod_radius_in_mm, 1500);

            for (int row = 0; row < n_row; row++)
            {
                // generateCodeSingleRod(row, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[row];
                TV from, to, left, right;
                rod->x(rod->indices.front(), from); rod->x(rod->indices.back(), to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                scaleAndShift(from); scaleAndShift(to);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                from[dim-1] = rod_radius_in_mm; to[dim-1] = rod_radius_in_mm;
                
                if (row % 2 == 0)
                {
                    moveTo(left, 3000, true);
                    writeLine(left, from, rod_radius_in_mm, 600);
                    writeLine(from, to, rod_radius_in_mm, 1200);
                    writeLine(to, right, rod_radius_in_mm, 600);
                    TV extend = right;
                    extend += (right - left) * 0.1;
                    moveTo(extend);
                }
                else
                {
                    moveTo(right, 3000, true);
                    writeLine(right, to, rod_radius_in_mm, 600);
                    writeLine(to, from, rod_radius_in_mm, 1200);
                    writeLine(from, left, rod_radius_in_mm, 600);
                    TV extend = left;
                    extend -= (right - left) * 0.1;
                    moveTo(extend);
                }
            }

            for (int col = 0; col < n_col; col++)
            {
                // generateCodeSingleRod(col + n_row, scaleAndShift, true, rod_radius_in_mm, 3.0 * rod_radius_in_mm, 0.3, false, true);
                auto rod = sim.Rods[col + n_row];
                TV from, to, left, right;
                rod->x(rod->indices[rod->dof_node_location[0]], from); rod->x(rod->indices[rod->dof_node_location.back()], to);
                left = from - (to - from) * extend_percentage;
                right = to + (to - from) * extend_percentage;

                scaleAndShift(left); scaleAndShift(right);
                left[dim-1] = rod_radius_in_mm; right[dim-1] = rod_radius_in_mm;
                scaleAndShift(from); scaleAndShift(to);
                from[dim-1] = 3.0 * rod_radius_in_mm; to[dim-1] = 3.0 * rod_radius_in_mm;
                
                if (col % 2 == 0)
                {
                    
                    moveTo(right);
                    right[dim-1] = 0.2;
                    moveTo(right);
                    to[1] = top_right[1];
                    TV temp = 0.5 * (right + to);
                    temp[2] = 0.2;

                    writeLine(right, temp, rod_radius_in_mm, 200);
                    writeLine(temp, to, rod_radius_in_mm, 400);
                    writeLine(to, from, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (from + left);
                    temp[2] = 0.2;
                    writeLine(from, temp, rod_radius_in_mm, 200);
                    writeLine(temp, left, rod_radius_in_mm, 200);
                }
                else
                {
                    moveTo(left);
                    left[dim-1] = 0.2;
                    moveTo(left);
                    from[1] = bottom_left[1];
                    TV temp = 0.5 * (left + from);
                    temp[2] = 0.2;
                    writeLine(left, temp, rod_radius_in_mm, 200);
                    writeLine(temp, from, rod_radius_in_mm, 400);
                    writeLine(from, to, 4.0 * rod_radius_in_mm, 600);
                    temp = 0.5 * (to + right);
                    temp[2] = 0.2;
                    writeLine(to, temp, rod_radius_in_mm, 200);
                    writeLine(temp, right, rod_radius_in_mm, 200);
                }
            }

        }

       
        

        writeFooter();
    }
}


template<class T, int dim>
void GCodeGenerator<T, dim>::activeTexticleGCode(bool fused)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(40, 80);
    };


    if constexpr (dim == 3)
    {
        writeHeader();
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3;

        if (fused)
        {
            // step one bottom layers
            for (int rod_idx  : {0, 1, 2, 3, 4, 5, 6, 8, 9})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
            
            for (int rod_idx  : {7})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
                    4.0 * rod_radius_in_mm);

            
            TV3 heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);

            std::vector<int> crossings = {7, 12, 17, 22};

            for (int crossing_id : crossings)
            {
                addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.03, 0.03));
            }

            
        }
        else
        {
            // step one bottom layers
            for (int rod_idx  : {0, 1, 2, 3, 4, 5, 9})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
            
            for (int rod_idx  : {6, 7, 8})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
                    4.0 * rod_radius_in_mm);

            
            TV3 heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);

            std::vector<int> crossings = {6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23};

            for (int crossing_id : crossings)
            {
                addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.03, 0.03));
            }
        }
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::activeTexticleGCode2(bool fused)
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(40, 80);
    };


    if constexpr (dim == 3)
    {
        writeHeader();
        T rod_radius_in_mm = sim.Rods[0]->a * 1e3;

        // if (fused)
        // {
        //     for (int rod_idx  : {5})
        //         generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
        //             rod_radius_in_mm);

        //     // step one bottom layers
        //     for (int rod_idx  : {0, 1, 2, 3, 4, 6, 8, 9})
        //         generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
            
        //     for (int rod_idx  : {7})
        //         generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
        //             4.0 * rod_radius_in_mm);

        //     TV3 heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);

        //     std::vector<int> crossings = {7, 12, 17, 22};

        //     for (int crossing_id : crossings)
        //     {
        //         addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.03, 0.03));
        //     }

        //     crossings = {10};
        //     for (int crossing_id : crossings)
        //     {
        //         addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 1, scaleAndShift, Range(0.03, 0.03));
        //     }
            
        // }
        if (fused)
        {
            // step one bottom layers
            for (int rod_idx  : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
               
        }
        else
        {

            for (int rod_idx  : {5})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
                    rod_radius_in_mm);

            // step one bottom layers
            for (int rod_idx  : {0, 1, 2, 3, 4})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 4.0 * rod_radius_in_mm);

            for (int rod_idx  : {9})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
            
            for (int rod_idx  : {6, 7, 8})
                generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
                    4.0 * rod_radius_in_mm);

            
            TV3 heights = TV3(first_layer_height, first_layer_height, 7.0 * first_layer_height);

            std::vector<int> crossings = {6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23};

            for (int crossing_id : crossings)
            {
                addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.04, 0.04));
            }

            crossings = {5, 10, 15};
            for (int crossing_id : crossings)
            {
                addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 1, scaleAndShift, Range(0.04, 0.04));
            }
        }
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::crossingTest()
{
    if constexpr (dim == 3)
    {
        layer_height = 0.3;
        writeHeader();
        TV from = TV(20, 40, layer_height);
        TV to = TV(60, 40, layer_height);
        TV extend = TV(80, 40, layer_height);

        TV left(38, 40, layer_height);
        TV right(42, 40, layer_height);

        moveTo(from);
        writeLine(from, to, layer_height, 300);
        moveTo(extend);
        extend[dim - 1] += 4.0;
        moveTo(extend);
        extend[dim - 1] -= 4.0;

        left[dim - 1] += 2.0;
        moveTo(left);
        left[dim - 1] -= 2.0;
        // addSingleTunnel(left, right, 2.0);
        
        right[dim - 1] += 4.0;
        moveTo(right);
        right[dim - 1] -= 4.0;
        
        for (int i = 0; i < 9; i++)
        {
            from[dim - 1] += layer_height;
            left[dim - 1] += layer_height;
            
            from[dim - 1] += 4.0;
            moveTo(from);
            from[dim - 1] -= 4.0;
            moveTo(from);

            writeLine(from, left, layer_height, 300);
            left[dim - 1] += 2.0;
            moveTo(left);
            left[dim - 1] -= 2.0;
        }

        for (int i = 0; i < 9; i++)
        {
            to[dim - 1] += layer_height;
            right[dim - 1] += layer_height;
            
            right[dim - 1] += 4.0;
            moveTo(right);
            right[dim - 1] -= 4.0;
            moveTo(right);
            
            writeLine(right, to, layer_height, 300);
            to[dim - 1] += 2.0;
            moveTo(to);
            to[dim - 1] -= 2.0;
        }

        

        from = TV(40, 20, layer_height);
        to = TV(40, 60, layer_height);
        extend = TV(40, 80, layer_height);

        TV middle = TV(40, 40, layer_height * 4);
        TV middle0 = TV(40, 30, layer_height);
        TV middle1 = TV(40, 50, layer_height);
        
        for (int i = 0; i < 9; i++)
        {
            from[dim - 1] += layer_height;
            to[dim - 1] += layer_height;
            middle[dim - 1] += layer_height;
            middle0[dim - 1] += layer_height;
            middle1[dim - 1] += layer_height;

            from[dim - 1] += 4.0;
            moveTo(from);
            from[dim - 1] -= 4.0;
            moveTo(from);

            writeLine(from, middle0, layer_height, 300);
            writeLine(middle0, middle, layer_height, 300);
            writeLine(middle, middle1, layer_height, 300);
            writeLine(middle1, to, layer_height, 300);

            to[dim - 1] += 2.0;
            moveTo(to);
            to[dim - 1] -= 2.0;
        }
    
        left[dim - 1] += 8.0;
        moveTo(left);
        left[dim - 1] -= 8.0;
        right[0] += 5;
        addSingleTunnel(left, right, 3.0);
        
        writeFooter();
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::addSingleTunnel(const TV& from, const TV& to, T height, T rod_diameter)
{
    if constexpr (dim == 3)
    {
        TV3 mid_point = TV3::Zero();
        mid_point.template segment<dim>(0) = 0.5 * (from + to);
        mid_point[2] += height;    
        moveTo(from);
        
        writeLine(from, mid_point, 0.2, 150);
        writeLine(mid_point, to, 0.2, 200);
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::addSingleTunnel(const TV& left, const TV& center, const TV& right, T rod_diameter)
{
    if constexpr (dim == 3)
    {
        TV _left = left, _right = right;
        _left[dim - 1] += 2.0;
        moveTo(_left);
        _left[dim - 1] -= 2.0;
        moveTo(_left, 100);
        
        writeLine(_left, center, tunnel_height, 400);
        writeLine(center, _right, tunnel_height, 600);
        _right[dim - 1] += 2.0;
        moveTo(_right);

        // moveTo(left, 600);
        
        // writeLine(left, center, 0.2, 600);
        // writeLine(center, right, 0.2, 600);
    }
}


template<class T, int dim>
void GCodeGenerator<T, dim>::generateCodeSingleRodMoveUpNozzle(int rod_idx, std::function<void(TV&)> scaleAndShift, 
        bool is_first_layer, T bd_height, T inner_height, T buffer_percentage, T less, T extend)
{
    auto rod = sim.Rods[rod_idx];

    T rod_radius_in_mm;
    rod_radius_in_mm = rod->a * 1e3 * 0.5;
    
    TV x0; rod->x(rod->indices.front(), x0);
    TV front, back;
    rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);

    TV extension = back + (back - front).normalized() * 0.6 * (front - back).norm();
    scaleAndShift(extension);

    //move slightly out of domain in case it doesn't stick at the beginning
    x0 -= (back - front).normalized() * buffer_percentage * (front - back).norm();
    scaleAndShift(x0);
    TV front_scaled = front;
    scaleAndShift(front_scaled);
    front_scaled[dim-1] = bd_height;

    // 2.0 is to avoid nozzle touching existing rods
    // x0[dim - 1] = 1.5;
    // retract(current_E - 0.5);
    // moveTo(x0);

    // // 0.2 is used for better sticking at the beginning
    // x0[dim - 1] = 0.2;
    // moveTo(x0, 200, false);
    // retract(current_E + 0.5);

    
    x0[dim - 1] = bd_height;
    moveTo(x0);

    if (buffer_percentage > 1e-6)
        writeLine(x0, front_scaled, rod_radius_in_mm);


    // writeLine(x0, front_scaled, rod_radius_in_mm);
    int running_cnt =0;

    std::vector<bool> is_fused;
    rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
    {
        is_fused.push_back(rod->isFixedNodeForPrinting(node_i, rod_idx));
        if (rod_idx == rod->numSeg() - 1)
            is_fused.push_back(rod->isFixedNodeForPrinting(node_j, rod_idx));        
        // std::cout << "is_fused " << std::endl;
        if (rod_idx == 0)
            std::cout << is_fused.back() << std::endl;
    }); 
    
    bool lift_head = false;

    std::vector<bool> fused_buffer = is_fused;
    for (int i = 0; i < is_fused.size(); i++)
    {
        if (!is_fused[i])
        {
            // lift_head = true;
            // break;       
            for (int j = i  ; j < i + 2; j++)
            {
                if (j >= 0 && j < rod->numSeg())
                {
                    fused_buffer[j] = false;
                }
            }
        }
    }



    int node_cnt = 0;
    rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
    {
        // std::cout << is_fused[node_cnt] << std::endl;
        TV xi, xj;
        rod->x(node_i, xi); rod->x(node_j, xj);
        // if (rod_idx == rod->numSeg() - 1 && buffer_percentage > 1e-6)
        //     xj += (back - front).normalized() * buffer_percentage * (front - back).norm();
        // else if (rod_idx == 0 && buffer_percentage > 1e-6)
        //     xi -= (back - front).normalized() * buffer_percentage * (front - back).norm();
        scaleAndShift(xi); scaleAndShift(xj);
        
        if (fused_buffer[node_cnt]) xi[dim - 1] = bd_height;
        else xi[dim - 1] = inner_height; 
        
        if (fused_buffer[node_cnt + 1]) xj[dim - 1] = bd_height;
        else xj[dim - 1] = inner_height; 
        node_cnt++;

        
        if (less)
            writeLine(xi, xj, 0.8 * rod_radius_in_mm, feed_rate_print);
        else
            writeLine(xi, xj, rod_radius_in_mm, feed_rate_print);
    });
    
    
    TV xn = back + (back - front).normalized() * buffer_percentage * (front - back).norm();
    // 
    scaleAndShift(xn);
    scaleAndShift(back);
    xn[dim - 1] = bd_height;
    // // moveTo(xn, 100);
    if (buffer_percentage > 1e-6)
        writeLine(back, xn, rod_radius_in_mm);
    // xn[dim - 1] = 2.0;
    // retract(current_E - 0.5);
    // moveTo(xn, 500, false);

    // // move nozzle along printing direction to avoid detaching of current print
    extension[dim - 1] = bd_height;
    if (extend)
        moveTo(extension);
    // retract(current_E + 0.5);
}


template<class T, int dim>
void GCodeGenerator<T, dim>::generateCodeSingleRod(int rod_idx, 
    std::function<void(TV&)> scaleAndShift, bool is_first_layer,
    T bd_height, T inner_height, T buffer_percentage, T less, T extend)
{
    auto rod = sim.Rods[rod_idx];

    T rod_radius_in_mm;
    rod_radius_in_mm = rod->a * 1e3 * 2.0;
    
    TV x0; rod->x(rod->indices.front(), x0);
    TV front, back;
    rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);

    TV extension = back + (back - front).normalized() * 0.6 * (front - back).norm();
    scaleAndShift(extension);

    //move slightly out of domain in case it doesn't stick at the beginning
    x0 -= (back - front).normalized() * buffer_percentage * (front - back).norm();
    scaleAndShift(x0);
    TV front_scaled = front;
    scaleAndShift(front_scaled);
    front_scaled[dim-1] = bd_height;

    // 2.0 is to avoid nozzle touching existing rods
    // x0[dim - 1] = 1.5;
    // retract(current_E - 0.5);
    // moveTo(x0);

    // // 0.2 is used for better sticking at the beginning
    // x0[dim - 1] = 0.2;
    // moveTo(x0, 200, false);
    // retract(current_E + 0.5);

    
    x0[dim - 1] = bd_height;
    moveTo(x0);

    if (buffer_percentage > 1e-6)
        writeLine(x0, front_scaled, rod_radius_in_mm);


    // writeLine(x0, front_scaled, rod_radius_in_mm);
    int running_cnt =0;

    std::vector<bool> is_fused;
    rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
    {
        is_fused.push_back(rod->isFixedNodeForPrinting(node_i, rod_idx));
        if (rod_idx == rod->numSeg() - 1)
            is_fused.push_back(rod->isFixedNodeForPrinting(node_j, rod_idx));        
        // std::cout << "is_fused " << std::endl;
        // std::cout << is_fused.back() << std::endl;
    }); 
    
    bool lift_head = false;

    std::vector<bool> fused_buffer = is_fused;
    for (int i = 0; i < is_fused.size(); i++)
    {
        if (!is_fused[i])
        {
            // lift_head = true;
            // break;       
            for (int j = i - 2 ; j < i + 2; j++)
            {
                if (j >= 0 && j < rod->numSeg())
                {
                    fused_buffer[j] = false;
                }
            }
        }
    }



    int node_cnt = 0;
    rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
    {
        // std::cout << is_fused[node_cnt] << std::endl;
        TV xi, xj;
        rod->x(node_i, xi); rod->x(node_j, xj);
        // if (rod_idx == rod->numSeg() - 1 && buffer_percentage > 1e-6)
        //     xj += (back - front).normalized() * buffer_percentage * (front - back).norm();
        // else if (rod_idx == 0 && buffer_percentage > 1e-6)
        //     xi -= (back - front).normalized() * buffer_percentage * (front - back).norm();
        scaleAndShift(xi); scaleAndShift(xj);
        
        if (fused_buffer[node_cnt]) xi[dim - 1] = bd_height;
        else xi[dim - 1] = inner_height; 
        
        if (fused_buffer[node_cnt + 1]) xj[dim - 1] = bd_height;
        else xj[dim - 1] = inner_height; 
        node_cnt++;

        // if (lift_head && rod_idx != 0 && rod_idx < rod->numSeg()-2)
        // {
        //      xj[dim - 1] = inner_height; 
        //      xi[dim - 1] = inner_height; 
        // }
        // else
        // {
        //     xj[dim - 1] = bd_height; 
        //     xi[dim - 1] = bd_height; 
        // }
        if (less)
            writeLine(xi, xj, 0.8 * rod_radius_in_mm, feed_rate_print);
        else
            writeLine(xi, xj, rod_radius_in_mm, feed_rate_print);
    });
    
    
    TV xn = back + (back - front).normalized() * buffer_percentage * (front - back).norm();
    // 
    scaleAndShift(xn);
    scaleAndShift(back);
    xn[dim - 1] = bd_height;
    // // moveTo(xn, 100);
    if (buffer_percentage > 1e-6)
        writeLine(back, xn, rod_radius_in_mm);
    // xn[dim - 1] = 2.0;
    // retract(current_E - 0.5);
    // moveTo(xn, 500, false);

    // // move nozzle along printing direction to avoid detaching of current print
    extension[dim - 1] = bd_height;
    if (extend)
        moveTo(extension);
    // retract(current_E + 0.5);
}


// template<class T, int dim>
// void GCodeGenerator<T, dim>::generateCodeSingleRod(int rod_idx, 
//     std::function<void(TV&)> scaleAndShift, bool is_first_layer,
//     T bd_height, T inner_height, T buffer_percentage)
// {
//     auto rod = sim.Rods[rod_idx];

//     T rod_radius_in_mm;
//     rod_radius_in_mm = rod->a * 1e3 * 2.0;
    
//     TV x0; rod->x(rod->indices.front(), x0);
//     TV front, back;
//     rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);

//     TV extend = back + (back - front).normalized() * buffer_percentage * (front - back).norm();
//     scaleAndShift(extend);

//     //move slightly out of domain in case it doesn't stick at the beginning
//     x0 -= (back - front).normalized() * buffer_percentage * (front - back).norm();
//     scaleAndShift(x0);
//     TV front_scaled = front;
//     scaleAndShift(front_scaled);
//     front_scaled[dim-1] = bd_height;

//     // 2.0 is to avoid nozzle touching existing rods
//     // x0[dim - 1] = 2.0;
//     // retract(current_E - 0.5);
//     // moveTo(x0, 2000, false);

//     // // 0.2 is used for better sticking at the beginning
//     // x0[dim - 1] = 0.2;
//     // moveTo(x0, 50, false);
//     // retract(current_E + 0.5);
//     // if (buffer_percentage > 1e-6)
//     //     writeLine(x0, front_scaled, rod_radius_in_mm);

//     x0[dim - 1] = rod_radius_in_mm;
//     moveTo(x0);

//     // writeLine(x0, front_scaled, rod_radius_in_mm);
//     int running_cnt =0;

//     std::vector<bool> is_fused;
//     rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
//     {
//         is_fused.push_back(rod->isFixedNodeForPrinting(node_i, rod_idx));
//         if (rod_idx == rod->numSeg() - 1)
//             is_fused.push_back(rod->isFixedNodeForPrinting(node_j, rod_idx));        
//         // std::cout << "is_fused " << std::endl;
//         // std::cout << is_fused.back() << std::endl;
//     }); 
    
//     std::vector<bool> fused_buffer = is_fused;
//     for (int i = 0; i < is_fused.size(); i++)
//     {
//         if (!is_fused[i])
//         {
            
//             for (int j = i - 5 ; j < i + 6; j++)
//             {
//                 if (j >= 0 && j < rod->numSeg())
//                 {
//                     fused_buffer[j] = false;
//                 }
//             }
//         }
//     }

//     int node_cnt = 0;
//     rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
//     {
//         // std::cout << is_fused[node_cnt] << std::endl;
//         TV xi, xj;
//         rod->x(node_i, xi); rod->x(node_j, xj);
//         // if (rod_idx == rod->numSeg() - 1)
//         //     xj += (back - front).normalized() * 0.05 * (front - back).norm();
//         scaleAndShift(xi); scaleAndShift(xj);
//         if (is_fused[node_cnt]) xi[dim - 1] = bd_height;
//         else xi[dim - 1] = inner_height; 
        
//         if (is_fused[node_cnt + 1]) xj[dim - 1] = bd_height;
//         else xj[dim - 1] = inner_height; 
//         node_cnt++;

//         // if (rod_idx > rod->numSeg() - 2 || rod_idx < 2)
//         // {
//         //     xi[dim - 1] = bd_height;
//         //     xj[dim - 1] = bd_height;
//         // }
//         // else
//         // {
//         //     xi[dim - 1] = inner_height;
//         //     xj[dim - 1] = inner_height;
//         // }
//         writeLine(xi, xj, rod_radius_in_mm);
//     });
    
    
//     TV xn = back + (back - front).normalized() * buffer_percentage * (front - back).norm();
//     // 
//     // scaleAndShift(xn);
//     // scaleAndShift(back);
//     // xn[dim - 1] += 0.2;
//     // // moveTo(xn, 100);
//     // if (buffer_percentage > 1e-6)
//     //     writeLine(back, xn, rod_radius_in_mm);
//     // xn[dim - 1] = 2.0;
//     // retract(current_E - 0.5);
//     // moveTo(xn, 100, false);

//     // // move nozzle along printing direction to avoid detaching of current print
//     // extend[dim - 1] = 2.0;
//     // moveTo(extend, 2000, false);
//     // retract(current_E + 0.5);
// }


template<class T, int dim>
void GCodeGenerator<T, dim>::addTunnelsAlongRod(int rod_idx, const TV3& heights,
        std::function<void(TV&)> scaleAndShift, int cross_n_seg,
        T extend_right, T speed_first_half, T speed_second_half)
{
    if constexpr (dim == 3)
    {
        auto rod = sim.Rods[rod_idx];
        int cnt = -1;

        for (int idx : rod->dof_node_location)
        {
            // std::cout << rod->fixed_by_crossing[cnt]  << std::endl;
            cnt++;
            if (rod->fixed_by_crossing[cnt])
                continue;
            int extend_seg = extend_right * cross_n_seg;
            int left_node, right_node;
            left_node = rod->closed ? (idx - cross_n_seg + rod->numSeg()) % (rod->numSeg()) : std::max(0, idx-cross_n_seg);
            right_node = rod->closed ? (idx + cross_n_seg + extend_seg) % (rod->numSeg()) : std::min(idx+cross_n_seg + extend_seg, rod->numSeg());

            TV left, right, center;
            rod->x(rod->indices[left_node], left);
            rod->x(rod->indices[idx], center);
            rod->x(rod->indices[right_node], right);

            scaleAndShift(left); scaleAndShift(right); scaleAndShift(center); 

            left[dim - 1] = heights[0];
            right[dim - 1] = heights[1];
            center[2] = heights[2];  

            left[dim - 1] += 2.0;
            moveTo(left);
            left[dim - 1] -= 2.0;
            moveTo(left, 100);
            
            writeLine(left, center, tunnel_height, speed_first_half);
            writeLine(center, right, tunnel_height, speed_second_half);
            right[dim - 1] += 2.0;
            moveTo(right);
        }
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::addSingleTunnelOnCrossingWithFixedRange(int crossing_id, const TV3& heights, 
        int direction, std::function<void(TV&)> scaleAndShift,
        const Range& range, T extend_right, T speed_first_half, T speed_second_half)
{
    if constexpr (dim == 3)
    {
        
        auto crossing = sim.rod_crossings[crossing_id];
        auto rod = sim.Rods[crossing->rods_involved[direction]];
        TV front, back;
        rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);
        T rod_length = (front - back).norm();
        TV rod_dir = (back - front).normalized();
        Range absolute_dx_in_mm = range * rod_length;
        TV x_crossing;
        rod->x(crossing->node_idx, x_crossing);
        TV left = x_crossing - rod_dir * absolute_dx_in_mm(0);
        TV right = x_crossing + rod_dir * absolute_dx_in_mm(1);

        scaleAndShift(left); scaleAndShift(right); 
        left[dim - 1] = heights[0];
        right[dim - 1] = heights[1];

        TV3 mid_point = TV3::Zero();
        mid_point.template segment<dim>(0) = 0.5 * (left + right);
        right += (right - left).normalized() * extend_right * (right - left).norm();
        mid_point[2] += heights[2];  
        left[dim - 1] += 2.0;
        moveTo(left);
        left[dim - 1] -= 2.0;
        moveTo(left, 100);
        
        writeLine(left, mid_point, tunnel_height, speed_first_half);
        writeLine(mid_point, right, tunnel_height, speed_second_half);
        right[dim - 1] += 2.0;
        moveTo(right);
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::addSingleTunnelOnCrossing(int crossing_id, const TV3& heights,
    int direction,
    std::function<void(TV&)> scaleAndShift)
{
    if constexpr (dim == 3)
    {
        auto crossing = sim.rod_crossings[crossing_id];
        Range range = crossing->sliding_ranges[direction];
        auto rod = sim.Rods[crossing->rods_involved[direction]];
        TV front, back;
        rod->x(rod->indices.front(), front); rod->x(rod->indices.back(), back);
        T rod_length = (front - back).norm();
        TV rod_dir = (back - front).normalized();
        Range absolute_dx_in_mm = range * rod_length;
        TV x_crossing;
        rod->x(crossing->node_idx, x_crossing);
        TV left = x_crossing - rod_dir * absolute_dx_in_mm(0);
        TV right = x_crossing + rod_dir * absolute_dx_in_mm(1);
        
        scaleAndShift(left); scaleAndShift(right); 
        left[dim - 1] = heights[0];
        right[dim - 1] = heights[1];

        TV3 mid_point = TV3::Zero();
        mid_point.template segment<dim>(0) = 0.5 * (left + right);
        right += (right - left).normalized() * 0.3 * (right - left).norm();
        mid_point[2] += heights[2];  
        left[dim - 1] += 2.0;
        moveTo(left);
        left[dim - 1] -= 2.0;
        moveTo(left);

        writeLine(left, mid_point, tunnel_height);
        writeLine(mid_point, right, tunnel_height);

        
        right[dim - 1] += 2.0;
        moveTo(right);
    }

}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeFromRodsShelterHardCoded()
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 50);
    };

    writeHeader();

    T rod_radius_in_mm = sim.Rods[0]->a * 1e3;
    for (int rod_idx  : {0, 1, 2, 3})
        generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
    
    for (int rod_idx  : {5, 6})
        generateCodeSingleRod(rod_idx, scaleAndShift, true,  0.8 * rod_radius_in_mm, 1.2 * rod_radius_in_mm);

    for (int rod_idx  : {4, 7})
        generateCodeSingleRod(rod_idx, scaleAndShift, true,  0.8 * rod_radius_in_mm, 4.0 * rod_radius_in_mm);
    
    TV3 heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);

    for (int crossing_id : {0, 3, 4, 7, 8, 11})
    {
        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.03, 0.03));
    }

    writeFooter();
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeFromRodsFixedGridGripperHardCoded()
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(40, 80);
    };

    writeHeader();

    T rod_radius_in_mm = sim.Rods[0]->a * 1e3;

    // step one bottom layers
    for (int rod_idx  : {0, 1, 2, 3, 4, 5, 6})
        generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
    
    for (int rod_idx  : {7, 8, 9, 11, 12, 13})
        generateCodeSingleRod(rod_idx, scaleAndShift, true, 1.2 * rod_radius_in_mm, 1.2 * rod_radius_in_mm);

    generateCodeSingleRod(10, scaleAndShift, true, rod_radius_in_mm, 
        4.0 * rod_radius_in_mm);

    // TV3 heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
    TV3 heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);
    for (int crossing_id : {3, 10})
    {
        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.03, 0.03));
    }

    heights = TV3(first_layer_height * 2.0, first_layer_height * 2.0, 20.0 * first_layer_height);
    // for (int crossing_id : {3, 10})
    // {
    //     addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 1, scaleAndShift, Range(0.06, 0.06));
    // }

    // addSingleTunnelOnCrossingWithFixedRange(3, heights, 1, scaleAndShift, Range(0.08, 0.08));
    // addSingleTunnelOnCrossingWithFixedRange(10, heights, 1, scaleAndShift, Range(0.04, 0.04));

    writeFooter();

}


template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeFromRodsGridGripperHardCoded()
{
    auto scaleAndShift = [](TV& x)->void
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(40, 80);
    };

    writeHeader();

    T rod_radius_in_mm = sim.Rods[0]->a * 1e3;

    // step one bottom layers
    for (int rod_idx  : {0, 1, 2, 3, 4, 5, 6, 7, 13})
        generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, rod_radius_in_mm);
    // step two second layers with tunnels in y
    // print heigher in the middle for easy removal of unfixed crossings
    for (int rod_idx  : {8, 9, 11, 12})
        generateCodeSingleRod(rod_idx, scaleAndShift, true, rod_radius_in_mm, 
            4.0 * rod_radius_in_mm);
    
    // add tunnel
    // T base_height = first_layer_height + layer_height;
    TV3 heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);
    for (int crossing_id : {8, 9, 11, 12, 15, 16, 18, 19})
    {
        // addSingleTunnelOnCrossing(crossing_id, heights1, 0, scaleAndShift);
        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.04, 0.04));
    }

    heights = TV3(first_layer_height, first_layer_height, 14.0 * first_layer_height);
    for (int crossing_id : {22, 25})
    {
        // addSingleTunnelOnCrossing(crossing_id, heights1, 0, scaleAndShift);
        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.1, 0.1));
    }

    generateCodeSingleRod(10, scaleAndShift, true, rod_radius_in_mm, 
        4.0 * rod_radius_in_mm);

    heights = TV3(rod_radius_in_mm, rod_radius_in_mm, 4.0 * rod_radius_in_mm);
    for (int crossing_id : {3, 10})
    {
        addSingleTunnelOnCrossingWithFixedRange(crossing_id, heights, 0, scaleAndShift, Range(0.03, 0.03));
    }

    
    writeFooter();
}

template<class T, int dim>
void GCodeGenerator<T, dim>::writeCircle(const TV& center, T r, const TV& start, const Vector<T, 2>& range,
    const std::vector<TV>& lifting_points, int sub_div, T rod_diameter)
{
    
    std::vector<T> thetas;
    for (int i = 0; i < sub_div + 1; i++)
    {
        thetas.push_back(range[0] + T(i) / T(sub_div) * (range[1] - range[0]));
    }

    moveTo(start);
    int cnt = 0;
    for (T theta : thetas)
    {
        TV next = center; 
        next[0] += r * std::cos(theta);
        next[1] += r * std::sin(theta);
        bool lift = false;
        for (auto pt : lifting_points)
        {
            if ((next.template head<2>() - pt.template head<2>()).norm() < 1.5)
            {
                next[dim-1] = 3.0 * rod_diameter;
                lift = true;
                break;
            }
            else
            {
                next[dim-1] = rod_diameter;
            }

        }
        if (lift)
            writeLine(current_position, next, 0.5 * rod_diameter, 600);
        else
        {
            if (cnt == 0)
                writeLine(current_position, next, rod_diameter, 200);
            else
                writeLine(current_position, next, rod_diameter, 600);
        }
        cnt ++;
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeFromRodsNoTunnel()
{
    auto scaleAndShift = [&](TV& x)
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 50);
    };

    writeHeader();
    for (auto& rod : sim.Rods)
    {
        TV x0; rod->x(rod->indices.front(), x0);
        scaleAndShift(x0);
        x0[dim - 1] = first_layer_height;
        moveTo(x0);
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            TV xi, xj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            
            scaleAndShift(xi); scaleAndShift(xj);
            
            T rod_radius_in_mm = rod->a * 1e3 * 2.0;
            writeLine(xi, xj, rod_radius_in_mm);
        });
    }
    writeFooter();
}

template<class T, int dim>
void GCodeGenerator<T, dim>::generateGCodeFromRodsCurveGripperHardCoded()
{
    auto scaleAndShift = [&](TV& x)
    {
        x *= 1e3;
        x.template segment<2>(0) += Vector<T, 2>(50, 50);
    };

    writeHeader();
    for (auto& rod : sim.Rods)
    {
        TV x0; rod->x(rod->indices.front(), x0);
        scaleAndShift(x0);
        x0[dim - 1] = first_layer_height;
        moveTo(x0);
        rod->iterateSegments([&](int node_i, int node_j, int rod_idx)
        {
            TV xi, xj;
            rod->x(node_i, xi); rod->x(node_j, xj);
            
            scaleAndShift(xi); scaleAndShift(xj);
            if (rod->rod_id == 0 || rod_idx == rod->numSeg() - 1)
            {
                xi[dim - 1] = first_layer_height;
                xj[dim - 1] = first_layer_height;
            }
            else
            {
                xi[dim - 1] = first_layer_height + layer_height * 4.0;
                xj[dim - 1] = first_layer_height + layer_height * 4.0;
            }
            T rod_radius_in_mm = rod->a * 1e3 * 2.0;
            writeLine(xi, xj, rod_radius_in_mm);
        });
    }
    int n0 = sim.Rods[0]->indices[1];
    int n1 = sim.Rods[0]->indices[sim.Rods[0]->indices.size()-2];
    TV left, right;
    sim.Rods[0]->x(n0, left);
    sim.Rods[0]->x(n1, right);
    scaleAndShift(left); scaleAndShift(right);
    left[dim - 1] = first_layer_height;
    right[dim - 1] = first_layer_height;
    addSingleTunnel(left, right, 2.0);
    writeFooter();
}

// ################################## UTILITIES ##################################

template<class T, int dim>
void GCodeGenerator<T, dim>::addRecBar(const TV& border_a, const TV& border_b, T width, T rod_diameter)
{
    if constexpr (dim == 3)
    {
        moveTo(border_a);
        TV line_vector = border_b - border_a;
        Eigen::Matrix3d rotation_clockwise = rotationMatrixFromEulerAngle(-M_PI / 2.0, 0.0, 0.0);
        TV ortho_dir = rotation_clockwise * (border_b - border_a);
        ortho_dir.normalize();
        // T dx = width / rod_diameter;
        // std::cout << dx << " " << width << " " << rod_diameter << std::endl;
        int cnt = 0;
        for (T delta = 0; delta < width; delta += rod_diameter)
        {
            TV from, to;
            if (cnt % 2 == 0)
            {
                from = current_position;
                to = current_position + line_vector;
                writeLine(from, to, rod_diameter);
                TV shift = to + ortho_dir * rod_diameter;
                moveTo(shift);
                // writeLine(to, shift, rod_diameter);
            }
            else
            {
                from = current_position;
                to = current_position - line_vector;
                writeLine(from, to, rod_diameter);
                TV shift = to + ortho_dir * rod_diameter;
                moveTo(shift);
                // writeLine(to, shift, rod_diameter);
            }
            cnt++;
            
        }
        
    }
    
}

template<class T, int dim>
void GCodeGenerator<T, dim>::writeLine(const TV& from, const TV& to, T rod_radius, T speed)
{
    
    T cross_section_area = crossSectionAreaFromRod(rod_radius);
	T amount = (to - from).norm();
    amount *= cross_section_area / (M_PI * filament_diameter * filament_diameter * 0.25);
    T current_amount = extrusion_mode == Absolute ? current_E + amount : amount;
    std::string cmd;
    if constexpr (dim == 3)
        cmd += "G1 F" + std::to_string(speed) + " X" + 
            std::to_string(to[0]) + " Y" + std::to_string(to[1]) +
            " Z" + std::to_string(to[2]) +
            " E" + std::to_string(current_amount) + "\n";
    else if constexpr (dim == 2)
        cmd += "G1 F" + std::to_string(speed) + " X" + 
        std::to_string(to[0]) + " Y" + std::to_string(to[1]) +
            " Z" + std::to_string(first_layer_height) +
            " E" + std::to_string(current_amount) + "\n";
    if (extrusion_mode == Absolute)
        current_E += amount;
    gcode << cmd;
    current_position = to;
}

template<class T, int dim>  
void GCodeGenerator<T, dim>::retract(T E)
{
    // gcode << "G1 E" << std::to_string(E) << " F2100.0" << std::endl;
    gcode << "G1 E" << std::to_string(E) << std::endl;
    current_E = E;
}

template<class T, int dim>  
void GCodeGenerator<T, dim>::extrude(T E)
{
    T current_amout = current_E + E;
    gcode << "G1 E" << std::to_string(current_amout) << " F2100.0" << std::endl;
    current_E = current_amout;
}

template<class T, int dim>  
void GCodeGenerator<T, dim>::moveTo(const TV& to, T speed, bool do_retract)
{
    if ((current_position - to).norm() < 1e-6)
        return;
    std::string cmd;
    if (extrusion_mode == Absolute)
    {
        if (do_retract)
            retract(current_E - 0.4);
        cmd += "G1 F" + std::to_string(speed) + " X" + 
            std::to_string(to[0]) + " Y" + std::to_string(to[1]) +
            " Z" + std::to_string(to[2]) + "\n";
        gcode << cmd;
        if (do_retract)
            retract(current_E + 0.4);
    }
    else if (extrusion_mode == Relative)
    {
        retract(-0.8);
        cmd += "G1 F" + std::to_string(speed) + " X" + 
            std::to_string(to[0]) + " Y" + std::to_string(to[1]) +
            " Z" + std::to_string(to[2]) + "\n";
        gcode << cmd;
        retract(0.8);
    }
    current_position = to;
}

template<class T, int dim>
T GCodeGenerator<T, dim>::extrusionWidth() const 
{
    return 1.2 * nozzle_diameter;
}

template<class T, int dim>
T GCodeGenerator<T, dim>::crossSectionArea(bool is_first_layer) const
{
	// Approximating cross sectino area, based on http://hydraraptor.blogspot.ch/2014/06/why-slicers-get-dimensions-wrong.html
	T extrusion_width = extrusionWidth();
	if (is_first_layer)
		return M_PI * first_layer_height * first_layer_height / 4 + first_layer_height * (extrusion_width - first_layer_height);
	else
		return M_PI * layer_height * layer_height / 4 + layer_height * (extrusion_width - layer_height);
}

template<class T, int dim>
T GCodeGenerator<T, dim>::crossSectionAreaFromRod(T rod_radius) const
{
	T extrusion_width = extrusionWidth();
	return M_PI * rod_radius * rod_radius / 4 + rod_radius * (extrusion_width - rod_radius);
}

template<class T, int dim>
void GCodeGenerator<T, dim>::writeHeader()
{
    gcode.open(gcode_file);
    if (printer == PrusaI3)
    {
        gcode << "G21 ; set units to millimeters\n";
        gcode << "G90 ; use absolute positioning\n";
        gcode << "M104 S" << std::to_string(extruder_temperature) << " ;set extruder temp\n";
        gcode << "M140 S"<< std::to_string(bed_temperature) << " ; set bed temp\n";
        gcode << "M190 S"<< std::to_string(bed_temperature) << " ; wait for bed temp\n";
        gcode << "M109 S" << std::to_string(extruder_temperature) << " ;wait for extruder temp\n";
        
        gcode << "G28 W ; home all without mesh bed level\n";
        gcode << "G80 ; mesh bed leveling\n";
        gcode << "G1 Y-3.0 F1000.0 ; go outside print area\n";
        gcode << "G92 E0.0 ; reset extruder distance position\n";
        gcode << "G1 X100.0 E9.0 F1000.0 ; intro line\n";
        gcode << "G92 E0.0 ; reset extruder distance position\n";
        if (extrusion_mode == Absolute)
            gcode << "M82 ;absolute extrusion mode\n";
        else if (extrusion_mode == Relative)
            gcode << "M83 ;relative extrusion mode\n";
        else
        {
            std::cout << "unexpected extrusion mode" << std::endl;
            std::exit(0);
        }
        
    }
    else
    {
        std::cout << "unrecognized printer type" << std::endl;
        std::exit(0);
    }
}

template<class T, int dim>
void GCodeGenerator<T, dim>::writeFooter()
{
    if (printer == PrusaI3)
    {
        gcode << "M107\n";
        gcode << "G4 ; wait\n";
        gcode << "M221 S100 ; reset flow\n";
        gcode << "M104 S0 ; turn off extruder\n";
        gcode << "M140 S0 ; turn off heatbed\n";
        gcode << "M107 ; turn off fan\n";
        gcode << "G1 Z33.6 ; Move print head up\n";
        gcode << "G1 X0 Y210; home X axis and push Y forward\n";
        gcode << "M84 ; disable motors\n";
        gcode << "M73 P100 R0\n";
        gcode << "M73 Q100 S0\n";
    }
    else
    {
        std::cout << "unrecognized printer type" << std::endl;
        gcode.close();
        std::exit(0);
    }
    gcode.close();
    
}

template class GCodeGenerator<double, 3>;
template class GCodeGenerator<double, 2>;  