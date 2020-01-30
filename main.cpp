#include "include.h"


class MC
{
public:
    int counter{0};
    mutable double volume{0.};
    double energy_diff;
    mutable double prob;
    mutable Segment molecules[n_molecules][n_segments];
    mutable double com[n_molecules][3];
    mutable Potential Pot;
    mutable ChangePos change_pos;
    mutable Observables Obs;
    mutable NPT npt;
    mutable Trig trig;
    Init init;
    
void prepare()   
{
    init.debug(molecules);
    trig.center_of_mass(com, molecules);
}

void run()
{   
    std::cout << molecules[0][0].C[0] << " first \n";
    std::string traj{"./traj/"};
    std::string id{"S_"};
    
    std::cout << beta <<"\n";
    id.append(std::to_string(n_segments));
    id.append("_N_");
    id.append(std::to_string(n_molecules));
    id.append("_T_");
    id.append(std::to_string(int( floor(1./beta) )));
    id.append("_P_");
    id.append(std::to_string(int( npt.pressure/(0.724297*pow(10,-7)))));
    id.append("_st_");
    id.append(std::to_string(n_steps));
    id.append(".");
    
    traj.append(id);
    traj.append("xyz");
    std::ofstream xyz;
  	xyz.open (traj);
    xyz.close();
    Obs.init();
    counter=0;
    volume=0;

    Pot.reset();
    Pot.inter(molecules);
    Pot.intra(molecules, com);
    Pot.old_pot = Pot.new_pot;
    std::ofstream bla;
  	bla.open ("pot.dat");
  	std::cout << 1/beta << " = T [K] \n";
    for (int step = 0; step < n_steps; ++step)
    {   
        if ((box_size[0] < 2.6*n_segments) || (box_size[1] < 2.6*n_segments) || (box_size[2] < 2.6*n_segments))
        {
            std::cout << "error box size to small \n";
            break; 
        }
        
        change_pos.choose_mol();
        change_pos.save(molecules, com);
        change_pos.propose(molecules, com);
        
        
        Pot.reset();
        Pot.inter(molecules);
        Pot.intra(molecules, com);
        
        energy_diff =  Pot.old_pot - Pot.new_pot;
        energy_diff *= beta;
        
        if (exp(energy_diff) > uniform())
        {   
            Pot.old_pot = Pot.new_pot;
        }
        else
        {   
            change_pos.reset(molecules, com);
        }

  
        std::cout << step << "\n";
        
        
        if ((npt_simu) && (step % n_molecules==0))
        {   
            
            npt.propose(com);
            for (int dim=0; dim<3; ++dim)
            {                  
                box_size[dim] *= npt.change;                
            }
            npt.change_vol(molecules, com, false);
            trig.center_of_mass(com, molecules);
            Pot.reset();
            Pot.intra(molecules, com);
            Pot.inter(molecules);
            prob = npt.criterion(Pot.new_pot, Pot.old_pot);
            
            if (prob > uniform())
            {                
                
                std::cout << float(step)/float(n_steps) << " " << Pot.old_pot<<" " << n_molecules/npt.new_vol <<"\n";
                Pot.old_pot = Pot.new_pot;
            }
            else
            {               
                npt.change_vol(molecules, com, true);
                for (int dim=0; dim<3; ++dim)
                {
                    box_size[dim] /= npt.change;
                }
                
                for (int mol=0; mol<n_molecules; ++mol)
                    for (int dim=0; dim<3; ++dim)
                        com[mol][dim] = npt.old_com[mol][dim];
            }
            
        }
        
        if ((step > n_steps*0.3))
        {   
            if (step%(10*n_molecules)==0)
            {   
                
                counter++;
                volume += box_size[0] * box_size[1] * box_size[2];
                
                Obs.get_gr(com);
                Obs.get_gee(molecules);
                Obs.get_rg(molecules, com);
                bla << Pot.old_pot << "\n";
            }
            if (step%(n_steps/2000)==0)
            {   
                std::cout << float(step)/float(n_steps) << "\n";
                write_xyz(molecules, com, traj);
            }
            
        }
        
        
    }
    bla.close();
    Obs.save_gr(id);
    Obs.save_gee(id);
    Obs.save_rg(id);
    std::string vol{"./vol/"};
    vol.append(id);
    vol.append("vol");
    std::ofstream volu;
  	volu.open (vol);
  	volu << volume/counter << "\n";
    volu.close();
    std::cout << molecules[0][0].C[0] << " lastmol \n";
}
    
};

int main()
{
    double betas[7]{1./700., 1./600., 1/500., 1/400., 1/300., 1./200., 1./100.};
    //double betas[1]{1/573.;};
    MC mc;
    mc.prepare();
    for (int b=0; b<7; ++b)
    {   
        beta = betas[b];
        std::cout << beta << "\n";
        mc.run();
    }

}