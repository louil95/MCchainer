class Observables

{

public:
    mutable Trig TrigObs;
    
    mutable double res_max;
    
    mutable double res{200};
       
    mutable double max_dist;
    
    mutable int count_calls{0};
    
    mutable long long int dist_hist[200];
    
    mutable long long int endend_hist[200];
    
    mutable long long int rg_hist[200];
     
    void init()
    {                 
        for (int i=0; i<res; ++i)
        {
            dist_hist[i] = 0;
            endend_hist[i] = 0;
            rg_hist[i] = 0;
            
        }
        count_calls=0;
        double min_box_len{box_size[0]};
        
        for (int dim=1; dim<3; ++dim)
        {
            if (min_box_len < box_size[dim])
            {
                min_box_len = box_size[dim];
            }
        }
        
        max_dist = min_box_len/2.;
        res_max = double(res)/max_dist;

    }
    
    void get_gr(double com[n_molecules][3])
    {   
        count_calls++;
        double distance;
        int index;
        
        for (int i=0; i < n_molecules-1; ++i)
        {
            for (int j=i+1; j<n_molecules; ++j)
            {   
                distance = TrigObs.get_dist(com[i], com[j]);
                distance = sqrt(distance);
                if (distance<max_dist)
                {   
                    distance = distance * res_max;
                    index = floor(distance);
                    dist_hist[index]++;
                }
            }
        }
        
    }
    
    void get_rg(Segment molecules[n_molecules][n_segments], double com[n_molecules][3])
    {   
        double distance;
        double local_com[3];
        int index;
        for (int mol=0; mol < n_molecules; ++mol)
        {   
            for (int dim=0; dim<3; ++dim)
            {
                local_com[dim] = com[mol][dim];
            }
            for (int seg=0; seg < n_segments; ++seg)
            {   
                distance = TrigObs.get_dist(local_com, molecules[mol][seg].C);
                distance = sqrt(distance);
                if (distance<max_dist)
                {   
                    distance = distance * res_max;
                    index = floor(distance);
                    rg_hist[index]++;
                }
            }
        }
     }
       
    
    void get_gee(Segment molecules[n_molecules][n_segments])
    {   
        double distance;
        int index;
        
        for (int mol=0; mol < n_molecules; ++mol)
        {

           distance = TrigObs.get_dist(molecules[mol][0].C, molecules[mol][n_segments-1].C);
           distance = sqrt(distance);
           if (distance<max_dist)
           {   
               distance = distance * res_max;
               index = floor(distance);
               endend_hist[index]++;
           }
        }
        
    }
    
    
    void save_gr(std::string id)
    {   
        double delta_vol;
        double radius;
        double g_r;
        double bin_width = max_dist/double(res);
        std::string gr_file{"./gr/"};
        gr_file.append(id);
        gr_file.append("gr");
        std::ofstream gr;
        gr.open (gr_file);
        for (int i=0; i < res;++i)
        {   
            radius = double(i)/double(res) * max_dist;
            delta_vol = 2./3. * M_PI;
            delta_vol *= pow(radius+bin_width, 3)- pow(radius, 3);
            g_r = double(dist_hist[i]);
            g_r *= box_size[0]*box_size[1]*box_size[2];
            g_r /= delta_vol;
            g_r /= n_molecules;
            g_r /= n_molecules;
            g_r /= count_calls;
            gr << radius << "\t" << g_r  << "\n";
            
        }
        gr.close();
    }
    
    void save_gee(std::string id)
    {
        std::string gee_file{"./gee/"};
        gee_file.append(id);
        gee_file.append("gee");
        std::ofstream gee;
        gee.open (gee_file);
        double g_ee, radius;
        for (int i=0; i < res;++i)
        {   
            radius = double(i)/double(res) * max_dist;
            g_ee = double(endend_hist[i]);
            g_ee /= double(count_calls);
            gee << radius << "\t" << g_ee  << "\n";
            
        }
        gee.close();
    }
    
    void save_rg(std::string id)
    {
        std::string rg_file{"./rg/"};
        rg_file.append(id);
        rg_file.append("rg");
        std::ofstream r_g;
        r_g.open (rg_file);
        double rg, radius;
        for (int i=0; i < res;++i)
        {   
            radius = double(i)/double(res) * max_dist;
            rg = double(rg_hist[i]);
            rg /= double(count_calls);
            rg /= n_molecules*n_segments;
            r_g << radius << "\t" << rg  << "\n";
            
        }
        r_g.close();
    }

};