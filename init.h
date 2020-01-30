//#include "pbc.h"

class Init
{   
        public:
        //double a{0.9};
        void linear(Segment molecules[n_molecules][n_segments])
        {
            for (int mol=0; mol < n_molecules; ++mol)
            {   
                double start[3];
                start[0] = uniform() * box_size[0] - 0.755 * n_segments;//;
                start[1] = uniform() * box_size[1] - 0.755 * n_segments;//;
                start[2] = uniform() * box_size[2] - 0.755 * n_segments;//;
                
                for (int seg=0; seg < n_segments; ++seg)
                {
                        
                        molecules[mol][seg].C[0] = start[0] + 1.51*seg;
                        molecules[mol][seg].C[0] = pbc(molecules[mol][seg].C[0], box_size[0]);
                        molecules[mol][seg].C[1] = start[1];
                        molecules[mol][seg].C[1] = pbc(molecules[mol][seg].C[1], box_size[1]);
                        molecules[mol][seg].C[2] = start[2];
                        molecules[mol][seg].C[2] = pbc(molecules[mol][seg].C[2], box_size[2]);
                }
                
            }
        }
        
        void debug(Segment molecules[n_molecules][n_segments])
        {
            for (int mol=0; mol < n_molecules; ++mol)
            {   
                molecules[mol][0].C[0] = uniform() * box_size[0] - 0.755 * n_segments;
                molecules[mol][0].C[1] = uniform() * box_size[1] - 0.755 * n_segments;
                molecules[mol][0].C[2] = uniform() * box_size[2] - 0.755 * n_segments;
                //std::cout << "warning debug start \n";
                //molecules[mol][0].C[0]  = box_size[0]/2 - 0.755 * n_segments;
                //molecules[mol][0].C[1]  = box_size[1]/2 - 0.755 * n_segments;
                //molecules[mol][0].C[2]  = box_size[2]/2 - 0.755 * n_segments;              
                molecules[mol][0].C[0] = pbc(molecules[mol][0].C[0], box_size[0]);
                molecules[mol][0].C[1] = pbc(molecules[mol][0].C[1], box_size[1]);
                molecules[mol][0].C[2] = pbc(molecules[mol][0].C[2], box_size[2]);
                
                for (int seg=1; seg < n_segments; ++seg)
                {
                        if (seg%2==0)
                        {
                            molecules[mol][seg].C[0] = molecules[mol][seg-1].C[0] + 1.51;                            
                            molecules[mol][seg].C[1] = molecules[mol][seg-1].C[1];                           
                        }                       
                        else
                        {
                            molecules[mol][seg].C[1] = molecules[mol][seg-1].C[1] + 1.51;
                            molecules[mol][seg].C[0] = molecules[mol][seg-1].C[0];
                        }

                        molecules[mol][seg].C[2] = molecules[mol][seg-1].C[2];
                        molecules[mol][seg].C[0] = pbc(molecules[mol][seg].C[0], box_size[0]);
                        molecules[mol][seg].C[1] = pbc(molecules[mol][seg].C[1], box_size[1]);
                        molecules[mol][seg].C[2] = pbc(molecules[mol][seg].C[2], box_size[2]);
                 }
            }                
        }

};