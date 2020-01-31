//#include "pbc.h"

class ChangePos
{   
    private:
        mutable int mol;
        int t_seg;
        double vec[3];
        double center[3];
        mutable Segment old_mol[n_segments];
        mutable double old_com[3];
        double rot_angle;


        
        
    public:
        mutable int style;
        Trig TrigChange;
        void save(Segment molecules[n_molecules][n_segments], double com[n_molecules][3])
        {   
            for (int dim=0; dim<3; ++dim)
            {   
                old_com[dim] = com[mol][dim];
                for (int seg=0; seg<n_segments; ++seg)
                {
                    old_mol[seg].C[dim] = molecules[mol][seg].C[dim];
                }
            }
        }
        
        void reset(Segment molecules[n_molecules][n_segments], double com[n_molecules][3])
        {   
            for (int dim=0; dim<3; ++dim)
            {   
                com[mol][dim] = old_com[dim];
                for (int seg=0; seg<n_segments; ++seg)
                {
                    molecules[mol][seg].C[dim] = old_mol[seg].C[dim];
                }
            }    
        }
        
        void choose_mol()
        {
            mol = floor(n_molecules * uniform());
        }
        
        void propose(Segment molecules[n_molecules][n_segments], double com[n_molecules][3])
        {   
            if (n_molecules!=1)
            {
                if (n_segments >= 4)
                {
                    style = floor(4 * uniform());
                }
                else if (n_segments >= 3)
                {
                    style = floor(3 * uniform());
                }
                else
                {
                    style = floor(2 * uniform());
                }
            }
            else
            {
                style  = 2 + floor(2 * uniform());
            }
            
            if (style == 0)
            {   
                for (int dim=0; dim<3; ++dim) 
                {
                    vec[dim] = uniform() - 0.5;
                }

                for (int dim=0; dim<3; ++dim)
                {               
                    for (int seg=0; seg<n_segments; ++seg)
                    {
                        molecules[mol][seg].C[dim] += vec[dim];
                        molecules[mol][seg].C[dim] = pbc(molecules[mol][seg].C[dim], box_size[dim]);
                    }
                }
            }
            else
            {
                for (int seg=0; seg<n_segments; ++seg)
                {   
                    for (int dim=0; dim<3; ++dim)
                    // shift molecule to origin get rid of the pbc
                    {
                        molecules[mol][seg].C[dim] = ipbc(molecules[mol][seg].C[dim], box_size[dim]/2, com[mol][dim]);
                        molecules[mol][seg].C[dim] -= com[mol][dim];
                    }
                }
                switch (style)
                {
                    case 1: // rotate;
                    {
                        int iaxis = floor(3. * uniform());
                        for (int dim=0; dim<3; ++dim)
                        {
                            vec[dim] = 0;                          
                            if (dim == iaxis)
                            {
                                vec[dim] = 1.;
                            }
                         }   
                        
                        rot_angle = uniform() - 0.5;
                        for (int seg=0; seg<n_segments; ++seg)
                        {   
                            TrigChange.prepare_rodriguez(rot_angle);
                            TrigChange.rodriguez(vec, molecules[mol][seg].C);
                        }
                        break;
                    }
                                            
                    case 2:
                    {   
                        t_seg = 1 + floor((n_segments-2) * uniform());
                        double bound_vec_0[3];
                        double bound_vec_1[3];
                        
                        for (int dim=0; dim<3; ++dim) 
                        {
                            bound_vec_0[dim] = molecules[mol][t_seg-1].C[dim] - molecules[mol][t_seg].C[dim];                           
                            bound_vec_1[dim] = molecules[mol][t_seg+1].C[dim] - molecules[mol][t_seg].C[dim];
                        }
                        
                        vec[0] = bound_vec_0[1]*bound_vec_1[2] - bound_vec_0[2]*bound_vec_1[1];
                        vec[1] = bound_vec_0[2]*bound_vec_1[0] - bound_vec_0[0]*bound_vec_1[2];
                        vec[2] = bound_vec_0[0]*bound_vec_1[1] - bound_vec_0[1]*bound_vec_1[0];
                        if ((vec[0] + vec[1] +vec[2]) ==0)
                        {
                            for (int dim=0; dim<3; ++dim)
                            {
                                vec[dim] = uniform();
                                std::cout << "undefined rotation\n";
                            }
                        }
                        
                        rot_angle = uniform() - 0.5;
                        TrigChange.prepare_rodriguez(rot_angle);
                        
                        if (uniform() < 0.5)
                        {   //std::cout << "one side \n";
                            for (int seg=t_seg+1; seg<n_segments; ++seg)
                            {   
                                for (int dim = 0; dim<3; ++dim) 
                                {   
                                    molecules[mol][seg].C[dim] -= molecules[mol][t_seg].C[dim];
                                }
                                TrigChange.rodriguez(vec, molecules[mol][seg].C);
                                for (int dim = 0; dim<3; ++dim) 
                                {   
                                    molecules[mol][seg].C[dim] += molecules[mol][t_seg].C[dim];
                                }
                            }                        
                        }
                        else
                        {   
                            for (unsigned seg = t_seg; seg-- != 0;) 
                            {
                                for (int dim = 0; dim<3; ++dim) 
                                {
                                    molecules[mol][seg].C[dim] -= molecules[mol][t_seg].C[dim];
                                }
                                TrigChange.rodriguez(vec, molecules[mol][seg].C);
                                for (int dim = 0; dim<3; ++dim) 
                                {
                                    molecules[mol][seg].C[dim] += molecules[mol][t_seg].C[dim];                     
                                }
                            }
                        }
                        break;
                    
                    }
                
                    case 3:
                    {   
                        t_seg = 1 + floor((n_segments-3) * uniform());
                        for (int dim = 0; dim<3; ++dim)
                        {   
                            vec[dim] = molecules[mol][t_seg].C[dim] - molecules[mol][t_seg+1].C[dim];
                            center[dim] = molecules[mol][t_seg].C[dim] + molecules[mol][t_seg+1].C[dim];
                            center[dim] /= 2.;
                        }
                        rot_angle = M_PI*(uniform() - 0.5);
                        TrigChange.prepare_rodriguez(rot_angle);
                        
                        if (uniform() < 0.5)
                        {
                            for (int seg=t_seg+1; seg<n_segments; ++seg)
                            {   
                                for (int dim = 0; dim<3; ++dim) 
                                {
                                    molecules[mol][seg].C[dim] -= center[dim];
                                }
                                TrigChange.rodriguez(vec, molecules[mol][seg].C);
                                for (int dim = 0; dim<3; ++dim) 
                                {
                                    molecules[mol][seg].C[dim] += center[dim];
                                    
                                }
                            }                        
                        }
                        else                        
                        {                       
                            for (unsigned seg = t_seg; seg-- != 0;) 
                            {   
                                for (int dim = 0; dim<3; ++dim) 
                                {
                                    molecules[mol][seg].C[dim] -= center[dim];
                                }
                                
                                TrigChange.rodriguez(vec, molecules[mol][seg].C);
                                for (int dim = 0; dim<3; ++dim) 
                                {
                                    molecules[mol][seg].C[dim] += center[dim];
                                    
                                }
                            }
                        }
                        break;
                    }                                
                }
                
                for (int seg=0; seg<n_segments; ++seg)
                {   
                    // shift molecule back to right place reapply pbc
                    for (int dim=0; dim<3; ++dim)
                    {
                        molecules[mol][seg].C[dim] += com[mol][dim];
                        molecules[mol][seg].C[dim] = pbc(molecules[mol][seg].C[dim], box_size[dim]);
                    }
                }
                
                
            }
            
            double new_com[3];
            TrigChange.one_com(new_com, molecules, mol);
            for (int dim=0; dim<3; ++dim)
            {
                com[mol][dim] = new_com[dim];
            }
            
        }
        
        bool check(Segment molecules[n_molecules][n_segments])
        {
        bool negative{false};
            for (int mol=0; mol<n_molecules; ++mol)
            {
                for (int seg=0; seg<n_segments; ++seg)
                {
            
                    for (int dim=0; dim<3; ++dim)
                        {
                            if (molecules[mol][seg].C[dim] <0)
                                negative = true;
                                //std::cout << style <<" "<< molecules[mol][seg].C[dim] << "\n";
                        }
                }
            }
            return negative;
        }
};        
