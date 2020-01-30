
class Potential
{   
    
public:
    mutable Trig TrigPot;
    mutable double new_pot{0.};
    mutable double old_pot{0.};
    double local_com[3];
  
    void reset()
    {
        new_pot = 0.;
    }
    
    void intra(Segment molecule[n_molecules][n_segments], double com[n_molecules][3])
    {
        double intra_pot{0};
        double local;

        for (int mol=0; mol < n_molecules; ++mol)
        {   
            local_com[0] = com[mol][0];
            local_com[1] = com[mol][1];
            local_com[2] = com[mol][2];
            for (int seg=0; seg < n_segments; ++seg)
            {   
                if (seg != 0 && seg != (n_segments - 1))
                {               
                local = TrigPot.get_angle(molecule[mol][seg-1].C,  molecule[mol][seg+1].C, molecule[mol][seg].C, local_com);
                
                if (local >= 1)
                {
                    local -= pow(10,-14);
                }
                
                if (local <= -1)
                {
                    local += pow(10,-14);
                }
                
                intra_pot += spring_pot(local);
                }
                
                if ((seg < (n_segments - 3)))
                {
                    local = TrigPot.get_tors(molecule[mol][seg].C, molecule[mol][seg+1].C, 
                                                   molecule[mol][seg+2].C, molecule[mol][seg+3].C, local_com);                   
                    intra_pot += torsion_pot(local);
                    
                    if ((seg < (n_segments - 4)))
                    {
                        for (int seg_j = seg+4; seg_j < n_segments; ++seg_j)
                        {   
                            local = TrigPot.get_dist(molecule[mol][seg].C, molecule[mol][seg_j].C);
                            
                            if (seg==0)
                            {
                                if (seg_j == (n_segments -1))
                                {
                                    intra_pot += lennard_jones(local,2);
                                }
                                else
                                {
                                    intra_pot += lennard_jones(local,1);
                                }
                            }
                            else if (seg_j == (n_segments -1))
                            {
                                intra_pot += lennard_jones(local,1);
                            }
                            else
                            {
                                intra_pot += lennard_jones(local,0);
                            }                            
                        }
                    }
                }
                
            
            }
        new_pot += intra_pot;
        }
    }
    
    void inter(Segment molecule[n_molecules][n_segments])
    {
    double inter_pot{0.};
    double sq_dist{0.};
    //if (TrigPot.half_box[0]!=7)
    //    std::cout << TrigPot.half_box[0] << "half box \n";
    for (int i=0; i<n_molecules-1; ++i)
    {
        for (int j=i+1; j<n_molecules; ++j)
        {
            for (int seg_i = 1; seg_i<n_segments-1; ++seg_i)
            {   
                
                sq_dist = TrigPot.get_dist(molecule[i][0].C, molecule[j][seg_i].C);
                inter_pot += lennard_jones(sq_dist,1);
                
                sq_dist = TrigPot.get_dist(molecule[i][n_segments-1].C, molecule[j][seg_i].C);
                inter_pot += lennard_jones(sq_dist,1);
                
                sq_dist = TrigPot.get_dist(molecule[i][seg_i].C, molecule[j][0].C);
                inter_pot += lennard_jones(sq_dist,1);
                
                sq_dist = TrigPot.get_dist(molecule[i][n_segments-1].C, molecule[j][n_segments-1].C);
                inter_pot += lennard_jones(sq_dist,1);
                
                for (int seg_j = 1; seg_j<n_segments-1; ++seg_j)
                {
                    sq_dist = TrigPot.get_dist(molecule[i][seg_i].C, molecule[j][seg_j].C);
                    inter_pot += lennard_jones(sq_dist,0);
                    
                }
            }
            for (int seg_i=0; seg_i <n_segments; seg_i+=n_segments-1)
            {
                for (int seg_j=0; seg_j <n_segments; seg_j+=n_segments-1)
                    {
                    sq_dist = TrigPot.get_dist(molecule[i][seg_i].C, molecule[j][seg_j].C);
                    inter_pot += lennard_jones(sq_dist,2);
                    }
            }
    
        }
    }

    new_pot += inter_pot;
    }
        
        
private:
    
    double const lj_sig{15.4449}; //squared
    // CH2-CH2, CH2-CH3,CH3-CH3
    double const lj_eps[3]{188, 293, 456};   
    double const k_0{31250};
    double const b_angle_0{1.9897};    
    double const tors_param[4]{1009.99, 2018.95, 136.37, -3165.30};
        
    double lennard_jones(double sq_dist, int type)
    {   
        if (sq_dist > 38.61)
        {
            return 0.;
        }
        else
        {
            double att = pow(sq_dist/lj_sig, -3.);
            double rep = pow(att, 2.);
            return lj_eps[type] * (rep-att);
        }
        
    }

    double spring_pot(double b_angle)
    {
        float diff = b_angle - b_angle_0;
        diff *= diff;
        return k_0 * diff;
    }
    
    double torsion_pot(double cosangle)
    {   
        double tors_pot{0.};
        for (int i=0; i<4; ++i)
        {
            tors_pot += tors_param[i] * pow(cosangle, i);
        }
        return tors_pot;
    }
};