class NPT
{
public:
    mutable double pressure{101325 * 0.724297*pow(10,-7)};
    mutable double old_vol{0.};
    mutable double new_vol{0.};
    mutable double old_box[3];
    mutable double old_com[n_molecules][3];
    mutable double change{1.};
    Trig TrigNpt;
    
    double criterion(double new_pot, double old_pot)
    {   
        double prob = -new_pot;
        prob += old_pot;
        prob += pressure*(old_vol-new_vol);
        prob += n_molecules/beta * log(new_vol/old_vol);
        prob = exp(beta*prob);
        return prob;
    }
    
    void change_vol(Segment molecules[n_molecules][n_segments], double com[n_molecules][3], bool reset)
    {   
        double factor;
        factor = change - 1;
        for (int mol=0; mol<n_molecules; ++mol)
        {   
            for (int seg=0; seg<n_segments; ++seg)
            {   
                std::cout << "\t";
                for (int dim=0; dim<3; ++dim)
                {   
                    if (!reset)
                    {   
                        molecules[mol][seg].C[dim] = ipbc(molecules[mol][seg].C[dim], old_box[dim]/2, old_com[mol][dim]);
                        molecules[mol][seg].C[dim] += factor * old_com[mol][dim];
                        molecules[mol][seg].C[dim] = pbc(molecules[mol][seg].C[dim], box_size[dim]);
                        
                            
                    }
                    else
                    {   
                        molecules[mol][seg].C[dim] = ipbc(molecules[mol][seg].C[dim], box_size[dim]/2, com[mol][dim]);
                        molecules[mol][seg].C[dim] -= factor * old_com[mol][dim];
                        molecules[mol][seg].C[dim] = pbc(molecules[mol][seg].C[dim], old_box[dim]);
                    }
                }
            }                   
        }
    }
    
    void propose(double com[n_molecules][3])
    {

        change = 0.5 + uniform();
        old_vol = box_size[0] * box_size[1] * box_size[2];
        new_vol = old_vol * change;
        
        change = pow(change, 1./3.);
        
        for (int dim=0; dim<3; ++dim)
        {
            old_box[dim] = box_size[dim];
            for (int mol=0; mol<n_molecules; ++mol)
            {
                old_com[mol][dim] = com[mol][dim];
            }
        }
    }
    
    
};